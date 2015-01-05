#include <mitsuba/bidir/vertex.h>
#include <mitsuba/core/kdtree.h>
#include <mitsuba/bidir/edge.h>
#include <mitsuba/bidir/path.h>
#include <vector>
#include <utility>
#include <cmath>

#include <unistd.h>

#include "vcm.h"


MTS_NAMESPACE_BEGIN

typedef PathVertex* PathVertexPtr;
//    typedef VCMKDTree::IndexType    IndexType;
//    typedef VCMKDTree::SearchResult SearchResult;

enum EVCMVertexData {
	EdVCMData = 0,
	EdVCData = 1,
	EdVMData = 2
};

struct VCMVertex {
	VCMVertex() {
	}

	VCMVertex(const Point &p, const Intersection &intersection,
			const Spectrum &throughput, uint32_t pathLength,
			float dvcm, float dvc, float dvm,
			const BSDF *bsdf) :
			mPosition(p), mIntersect(intersection), mThroughput(throughput),
			mPathLength(pathLength), dVCM(dvcm),
			dVC(dvc), dVM(dvm), mBsdf(bsdf) {
	}

	Point mPosition;
	Intersection mIntersect;
	Spectrum mThroughput; // Path throughput (including emission)
	uint32_t mPathLength; // Number of segments between source and vertex

	// Stores all required local information, including incoming direction.
	// TODO: will it boom as it's a pointer?
	Float dVCM; // MIS quantity used for vertex connection and merging
	Float dVC;  // MIS quantity used for vertex connection
	Float dVM;  // MIS quantity used for vertex merging
	const BSDF *mBsdf;
};

// to make easier to carry those values around
struct SubPathState {
	Spectrum mThroughput;         // Path throughput
	bool mIsFiniteLight;
	bool mSpecularPath;
	uint32_t mPathLength;

	Float dVCM; // MIS quantity used for vertex connection and merging
	Float dVC;  // MIS quantity used for vertex connection
	Float dVM;  // MIS quantity used for vertex merging
};

struct VCMTreeEntry :
		public SimpleKDNode<Point, VCMVertex> {
public:
	/// Dummy constructor
	inline VCMTreeEntry() {
	}

	inline VCMTreeEntry(const VCMVertex &v) {
		position = v.mPosition;
		data = v;
	}

	/// Return a string representation (for debugging)
	std::string toString() const {
	}
};

class VCMIntegrator : public Integrator {
public:
	VCMIntegrator(const Properties &props) : Integrator(props) {
		/* Load the parameters / defaults */
		m_config.maxDepth = props.getInteger("maxDepth", 10);
		m_config.rrDepth = props.getInteger("rrDepth", 5);
		m_config.iterationCount = props.getInteger("iterationCount", 5);
		m_config.useVC = props.getBoolean("useVC", true);
		m_config.useVM = props.getBoolean("useVM", false);
		m_config.lightTraceOnly = props.getBoolean("lightTraceOnly", false);


		m_config.dump();

		if (m_config.rrDepth <= 0)
			Log(EError, "'rrDepth' must be set to a value greater than zero!");

		if (m_config.maxDepth <= 0 && m_config.maxDepth != -1)
			Log(EError, "'maxDepth' must be set to -1 (infinite) or a value greater than zero!");

	}

	/// Unserialize from a binary data stream
	VCMIntegrator(Stream *stream, InstanceManager *manager)
			: Integrator(stream, manager) {
		m_config = VCMConfiguration(stream);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Integrator::serialize(stream, manager);
		m_config.serialize(stream);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue,
			const RenderJob *job, int sceneResID, int sensorResID,
			int samplerResID) {
		Integrator::preprocess(scene, queue, job, sceneResID,
				sensorResID, samplerResID);

		return true;
	}

	void cancel() {
		Scheduler::getInstance()->cancel(m_process);
	}

	void configureSampler(const Scene *scene, Sampler *sampler) {
		/* Prepare the sampler for tile-based rendering */
		sampler->setFilmResolution(scene->getFilm()->getCropSize(), true);
	}

	bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
			int sceneResID, int sensorResID, int samplerResID) {

		Log(EInfo, "Start");

		mScene = scene;
		mSampler = scene->getSampler();
		// do we need this?
		//scene->initializeBidirectional();
		ref <Sensor> sensor = scene->getSensor();
		Film *film = sensor->getFilm();
		const Vector2i res = film->getSize();
		Float pathCount = res.x * res.y;
		mLightSubPathCount = float(res.x * res.y);

		mBitmaps.reserve(m_config.iterationCount);
		for (int i = 0; i < m_config.iterationCount; i++) {
			mBitmaps[i] = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, res);
			mBitmaps[i]->clear();
		}

		mBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getSize());
		mBitmap->clear();

		film->setBitmap(mBitmap);
		const Float radiusAlpha = 0.75f;

		for (int iterationNum = 0; iterationNum < m_config.iterationCount; iterationNum++) {

			Float radius = scene->getBSphere().radius * 0.003f;
			radius /= std::pow(Float(iterationNum + 1), 0.5f * (1 - radiusAlpha));
			const Float radiusSqr = radius * radius;

			Log(EInfo, "Scene BSphere radius: %f", scene->getBSphere().radius);
			Log(EInfo, "VCM Radius %f", radius);

			mEtaVCM = (M_PI * radiusSqr) * pathCount;
			mMisVMWeightFactor = 0.0f; //mEtaVCM;
			mMisVCWeightFactor = 1.f / mEtaVCM;


			Log(EInfo, "pathCount: %f", pathCount);
			Log(EInfo, "etaVCM: %f\nVM Weight: %f\nVC Weight: %f\n", mEtaVCM, mMisVMWeightFactor, mMisVCWeightFactor);
			//////////////////////////////////////////////////////////////////////////
			// Generate light paths
			//////////////////////////////////////////////////////////////////////////

			m_lightVertices.reserve(pathCount);
			m_pathEnds.reserve(pathCount);
			m_lightVertices.clear();
			m_pathEnds.clear();

			for (int i = 0; i < (int)pathCount; ++i) {
				mSampler->generate(Point2i(i));
				for (int advanceTimes = 0; advanceTimes < 20*iterationNum; advanceTimes++) {
					mSampler->advance();
				}

				SubPathState pathState;

				Ray ray = generateLightSample(pathState);

				for (;; ++pathState.mPathLength) {
					Intersection its;

					if (!mScene->rayIntersect(ray, its)) {
						break;
					}

					if (!its.isValid()) {
						break;
					}

					if (pathState.mPathLength > 1 || pathState.mIsFiniteLight) {
						pathState.dVCM *= its.t * its.t;
					}

					const BSDF *bsdf = its.getBSDF();
					Float cosTheta = std::abs(Frame::cosTheta(its.wi));

					pathState.dVCM /= cosTheta;
					pathState.dVC /= cosTheta;
					pathState.dVM /= cosTheta;

					if (!(bsdf->getType() & BSDF::EDelta)) {
						VCMVertex v(
								its.p,
								its,
								pathState.mThroughput,
								pathState.mPathLength,
								pathState.dVCM,
								pathState.dVC,
								pathState.dVM,
								its.getBSDF());

						m_lightVertices.push_back(v);
					}
					if (!(bsdf->getType() & BSDF::EDelta) &&
							(m_config.useVC || m_config.lightTraceOnly)) {
						connectToEye(pathState, its, bsdf, iterationNum);
					}

					if (pathState.mPathLength + 2 > m_config.maxDepth) {
						break;
					}

					if(!SampleScattering(ray, pathState, bsdf, its)) {
						break;
					}
				}

				m_pathEnds[i] = (int) m_lightVertices.size();
			}

			/////////////////////////////////////////////////////////////////////////
			// BUILD SEARCH STRUCT
			/////////////////////////////////////////////////////////////////////////
/*
			m_tree.reserve(pathCount);
			m_tree.clear();
			for (int i = 0; i < m_lightVertices.size(); i++) {
				m_tree.push_back(m_lightVertices[i]);
			}

			Log(EInfo, "Built the tree, it has %d vertices", m_tree.size());

*/
			/////////////////////////////////////////////////////////////////////////
			// GENERATE CAMERA PATHS
			/////////////////////////////////////////////////////////////////////////

			SubPathState pathState;
			for (int pathIdx = 0;
				(!m_config.lightTraceOnly && pathIdx < (int)pathCount);
				++pathIdx)
			{
				int currentX = pathIdx % res.x;
				int currentY = pathIdx / res.x;
				Point2i startPosition = Point2i(currentX, currentY);
				mSampler->generate(startPosition);

				for (int advanceTimes = 0; advanceTimes < 20*iterationNum; advanceTimes++) {
					mSampler->advance();
				}
				Ray ray = generateCameraSample(pathState, pathIdx, pathCount);

				Spectrum color(0.0f);
				for (;; pathState.mPathLength++) {
					Intersection its;

					if (!mScene->rayIntersect(ray, its)) {
						// TODO: add background emitter radiance
						break;
					}

					if (!its.isValid()) {
						break;
					}

					const BSDF *bsdf = its.getBSDF();

					Float cosTheta = std::abs(Frame::cosTheta(its.wi));
					pathState.dVCM *= its.t * its.t;
					pathState.dVCM /= cosTheta;
					pathState.dVC /= cosTheta;
					pathState.dVM /= cosTheta;

					// Case #1 - we hit the light
					if (its.isEmitter()) {
						Spectrum radiance = its.Le(its.wi);

						if (radiance.isZero()) {
							// if you don't try this direction then some lights
							// will be improperly displayed. Don't ask me why.
							// (try Veach mi.xml scene)
							radiance = its.Le(-its.wi);
							if (radiance.isZero())
								break;
						}

						if (pathState.mPathLength == 1) {
							color += pathState.mThroughput * radiance;
						} else {
							PositionSamplingRecord pRec(its);
							DirectionSamplingRecord dRec(-ray.d, ESolidAngle);
							const Emitter *emitter = its.shape->getEmitter();
							Float directPdf = emitter->pdfPosition(pRec);
							Float emissionPdf = directPdf * emitter->pdfDirection(dRec, pRec);

							// if emission is impossible, stop
							if (emissionPdf == 0)
								break;

							Float lightPickProb = LightPickProbability();
							directPdf *= lightPickProb;
							emissionPdf *= lightPickProb;
							const Float wSensor = directPdf * pathState.dVCM +
								emissionPdf * pathState.dVC;

							const Float misWeight = 1.f / (1.f + wSensor);

							color += pathState.mThroughput * misWeight * radiance;
						}
						break;
					}

					if (m_config.maxDepth != -1 && pathState.mPathLength + 1 > m_config.maxDepth) {
						break;
					}

					// Case #2 - Vertex Connection - connect to light
					// source
					if (!(bsdf->getType() & BSDF::EDelta)) {
						color += ConnectToLight(pathState, its, bsdf);
					}

					// Case #3 - Vertex Connection - Connect to light vertices
					if (!(bsdf->getType() & BSDF::EDelta)) {
						int startingLightVertex = pathIdx == 0 ? 0 : m_pathEnds[pathIdx-1];
						int lastLightVertex = m_pathEnds[pathIdx];

						for (int i = startingLightVertex; i < lastLightVertex; i++) {
							const VCMVertex &lightVertex = m_lightVertices[i];

							if (m_config.maxDepth != -1 &&
								lightVertex.mPathLength + pathState.mPathLength + 1 > m_config.maxDepth)
								break;

							color += pathState.mThroughput * lightVertex.mThroughput *
								ConnectVertices(pathState, lightVertex, its, bsdf);
						}
					}

					if(!SampleScattering(ray, pathState, bsdf, its)) {
						break;
					}
				}
				accumulateColor(currentX, currentY, color, iterationNum);
			}
			Log(EInfo, "Done iteration %d", iterationNum + 1);
			//accumulateMainBitmap(iterationNum, film, queue, job);
			film->addBitmap(mBitmaps[iterationNum], (Float) 1.0f / m_config.iterationCount);
			queue->signalRefresh(job);
		}

		Log(EInfo, "DONE!");
		return true;
	}

	Ray generateLightSample(SubPathState& pathState) {
		ref_vector<Emitter> emitters = mScene->getEmitters();
		int lightId = int(mSampler->next1D() * emitters.size());
		Float lightPickProb = (Float) 1.0f / emitters.size();
		ref<Emitter> emitter = emitters[lightId];


		PositionSamplingRecord pRec;
		DirectionSamplingRecord dRec;
		emitter->samplePosition(pRec, mSampler->next2D());
		emitter->sampleDirection(dRec, pRec, mSampler->next2D());

		Ray ray(pRec.p, dRec.d, 0); // assume no time sampling

		Float directPdfW = pRec.pdf;
		Float emissionPdfW = directPdfW * dRec.pdf;
		Float cosLight = dot(dRec.d, pRec.n);
		cosLight = std::max(cosLight, Epsilon);

		// this is a hack to get just the light radiance, may not work with any
		// type of light
		pathState.mThroughput = cosLight * emitter->evalPosition(pRec) / M_PI;

		directPdfW *= lightPickProb;
		emissionPdfW *= lightPickProb;

		pathState.mThroughput /= emissionPdfW;
		pathState.mPathLength = 1;
		pathState.mIsFiniteLight = emitter->needsDirectionSample();

		pathState.dVCM = directPdfW / emissionPdfW;

		if (!emitter->isDegenerate()) {
			Float usedCosLight = emitter->needsPositionSample() ? cosLight : 1.f;
			pathState.dVC = usedCosLight / emissionPdfW;
		} else {
			pathState.dVC = 0.0f;
		}

		pathState.dVM = pathState.dVC * mMisVCWeightFactor;

		return ray;
	}

	Ray generateCameraSample(SubPathState& pathState, int pixelIndex, Float pathCount) {
		Ray ray;
		const PerspectiveCamera* sensor =
			static_cast<const PerspectiveCamera*>(mScene->getSensor());
		const Vector2i res = sensor->getFilm()->getSize();
		int currentX = pixelIndex % res.x;
		int currentY = pixelIndex / res.x;
		Point2 startPosition(currentX, currentY);
		startPosition += mSampler->next2D() / 2;
		sensor->sampleRay(
				ray,
				startPosition,
				startPosition, // assume we don't need aperture sampling
				0			   // assume we don't need time sampling
				);

		PositionSamplingRecord pRec;
		sensor->samplePosition(pRec, mSampler->next2D());
		DirectionSamplingRecord dRec(ray.d, ESolidAngle);

		Float cameraPdfW = computeCameraPdfW(sensor, dRec, pRec);

		pathState.mPathLength = 1;
		pathState.mSpecularPath = true;
		pathState.mThroughput = Spectrum(1.0f);
		pathState.dVCM = pathCount / cameraPdfW;
		pathState.dVC = pathState.dVM = 0;

		//Log(EInfo, "CameraPdfW = %f, pathCount = %f, dVCM = %f", cameraPdfW, pathCount, pathState.dVCM);

		return ray;
	}

	Float computeCameraPdfW(
			const PerspectiveCamera* sensor,
			const DirectionSamplingRecord& dRec,
			const PositionSamplingRecord& pRec) {
		Float imagePlaneDist = sensor->getFilm()->getSize().x / (2.0f * std::tan(sensor->getXFov() * M_PI / 360));
		Float cameraPdfW = imagePlaneDist * imagePlaneDist * sensor->pdfDirection(dRec, pRec);
		return cameraPdfW;
	}

	void accumulateColor(int x, int y, const Spectrum& color, int iteration) {
		Spectrum *target = (Spectrum *) mBitmaps[iteration]->getUInt8Data();
		target[y * mBitmap->getWidth() + x] += color;
	}


	void accumulateMainBitmap(int iteration, Film *film, RenderQueue* queue, const RenderJob *job) {
		mBitmap->clear();
		Spectrum *target = (Spectrum*) mBitmap->getUInt8Data();
		Float iterationFactor = 1.0f / (1.0f+iteration);


		Log(EInfo, "Iteration factor: %f", iterationFactor);
		for (int i = 0; i <= iteration; i++) {
			Spectrum *currentItBitmap = (Spectrum*)mBitmaps[i]->getUInt8Data();
			for (int idx = 0; idx < mBitmap->getPixelCount(); idx++) {
				target[i] += currentItBitmap[i] * iterationFactor;
				if (idx % 31 == 0) {
					Log(EInfo, "accumulated color %s to target %s",
							currentItBitmap[i].toString().c_str(),
							target[i].toString().c_str());
				}
			}
		}
		film->setBitmap(mBitmap);
		queue->signalRefresh(job);
	}

	Spectrum ConnectVertices(SubPathState &pathState,
			const VCMVertex &lightVertex, const Intersection &its, const BSDF *bsdf) {

		Vector direction = lightVertex.mPosition - its.p;
		Float dist = direction.length();
		direction /= dist;

		Float cosCamera, cameraBsdfDirPdfW, cameraBsdfRevPdfW;
		BSDFSamplingRecord bsdfRec(its, its.toLocal(direction));
		Spectrum cameraBsdfFactor = bsdf->eval(bsdfRec);
		cosCamera = std::abs(Frame::cosTheta(bsdfRec.wo));

		if (cameraBsdfFactor.isZero()) {
			//Log(EInfo, "camera BSDF Factor is zero");
			return Spectrum(0.0f);
		}

		cameraBsdfDirPdfW = bsdf->pdf(bsdfRec);
		bsdfRec.reverse();
		cameraBsdfRevPdfW = bsdf->pdf(bsdfRec);

		Float continuationProbability = 1.0f;
		//cameraBsdfRevPdfW *= continuationProbability;
		//cameraBsdfDirPdfW *= continuationProbability;

		Float cosLight, lightBsdfDirPdfW, lightBsdfRevPdfW;
		BSDFSamplingRecord lightBsdfRec(lightVertex.mIntersect, lightVertex.mIntersect.toLocal(-direction));
		Spectrum lightBsdfFactor = lightVertex.mBsdf->eval(lightBsdfRec);
		cosLight = std::abs(Frame::cosTheta(lightBsdfRec.wo));

		lightBsdfDirPdfW = lightVertex.mBsdf->pdf(lightBsdfRec);
		lightBsdfRec.reverse();
		lightBsdfRevPdfW = lightVertex.mBsdf->pdf(lightBsdfRec);

		// it's not the proper value - we have to get it from light vertex :(
		Float lightCont = 1.0f;
		//lightBsdfDirPdfW *= lightCont;
		//lightBsdfRevPdfW *= lightCont;

		Float geometryTerm = 1.0f / (dist * dist);

		// Not sure how that's possible.
		if (cosCamera * cosLight < 0) {
			//Log(EInfo, "cosCamera: %f, cosLight: %f", cosCamera, cosLight);
			return Spectrum(0.0f);
		}

		Float cameraBsdfDirPdfA = pdfWtoA(cameraBsdfDirPdfW, dist, cosLight);
		Float lightBsdfDirPdfA = pdfWtoA(lightBsdfDirPdfW, dist, cosCamera);

		Float wLight = cameraBsdfDirPdfA * (mMisVMWeightFactor + lightVertex.dVCM
				+ lightVertex.dVC * lightBsdfRevPdfW);

		Float wCamera = lightBsdfDirPdfA * (mMisVMWeightFactor + pathState.dVCM
				+ pathState.dVC * cameraBsdfRevPdfW);

		Float misWeight = 1.0f / (wLight + 1.0f + wCamera);

		Spectrum contrib = (geometryTerm * misWeight) * cameraBsdfFactor * lightBsdfFactor;

		if (contrib.isZero() || occluded(its.p, direction, dist)) {
			//Log(EInfo, "contrib: %s, if non-zero, occluded!", contrib.toString().c_str());
			return Spectrum(0.0f);
		}
		//Log(EInfo, "Non zero contrib from vertex connection: %s, misWeight: %f", contrib.toString().c_str(), misWeight);

		return contrib;
	}

	bool occluded(const Point& point, const Vector& direction, Float distance) {
		Ray r(point, direction, Epsilon, distance * (1-ShadowEpsilon), 0);
		Intersection testIsect;

		return mScene->rayIntersect(r);
	}

	Float pdfWtoA(Float pdfW, Float dist, Float cosThere) {
		return pdfW * std::abs(cosThere) / (dist * dist);
	}

	Float LightPickProbability() {
		return 1.0f / mScene->getEmitters().size();
	}

	Spectrum ConnectToLight(SubPathState &pathState, const Intersection &its, const BSDF *bsdf) {
		Float lightPickProb = LightPickProbability();

		Point2 randomPoint = mSampler->next2D();

		Spectrum color(0.0f), radiance(0.0f);
		DirectSamplingRecord dRec(its);

		radiance = mScene->sampleEmitterDirect(dRec, randomPoint, true);
		radiance *= dRec.pdf;

		if (radiance.isZero()) {
			return color;
		}

		Float cosAtLight = dot(dRec.n, -dRec.d);
		const Emitter *sampledEmitter = static_cast<const Emitter *>(dRec.object);
		PositionSamplingRecord emitterPRec;
		sampledEmitter->samplePosition(emitterPRec, randomPoint);

		Float directPdfW = emitterPRec.pdf * dRec.dist * dRec.dist / cosAtLight;
		Float emissionPdfW = emitterPRec.pdf * cosAtLight * INV_PI;

		Float bsdfRevPdfW, bsdfDirPdfW, cosToLight;
		BSDFSamplingRecord bsdfRec(its, its.toLocal(dRec.d), ERadiance);
		cosToLight = std::abs(Frame::cosTheta(bsdfRec.wo));

		Spectrum bsdfFactor = bsdf->eval(bsdfRec);
		if (bsdfFactor.isZero()) {
			return color;
		}

		bsdfDirPdfW = bsdf->pdf(bsdfRec);
		bsdfRec.reverse();
		bsdfRevPdfW = bsdf->pdf(bsdfRec);

		Float continuationProbability = 1.0f ;

		if (sampledEmitter->isDegenerate()) {
			bsdfDirPdfW = 0.0f;
		}

		//bsdfDirPdfW *= continuationProbability;
		//bsdfRevPdfW *= continuationProbability;
		Float wLight = bsdfDirPdfW / (lightPickProb * directPdfW);
		Float wCamera = (emissionPdfW * cosToLight / (directPdfW * cosAtLight)
				* (mMisVMWeightFactor + pathState.dVCM + pathState.dVC * bsdfRevPdfW));
		//Log(EInfo, "\n  emissionPdfW: %f\n  cosToLight: %f\n  directPdfW: %f\n  cosAtLight: %f\n  mMisVMWeightFactor: %f\n  dVCM: %f\n  dVC: %f\n  bsdfRevPdfW: %f\n", emissionPdfW, cosToLight, directPdfW, cosAtLight, mMisVMWeightFactor, pathState.dVCM, pathState.dVC, bsdfRevPdfW);
		Float misWeight = 1.f / (wLight + 1.f + wCamera);
		//Log(EInfo, "wLight: %f, wCamera: %f, emissionPdfW: %f, directPdfW: %f\ndVCM: %f, dVC: %f, bsdfRevPdfW: %f",
		//	wLight, wCamera, emissionPdfW, directPdfW, pathState.dVCM, pathState.dVC, bsdfRevPdfW);
		//Log(EInfo, "wLight: %f, wCamera: %f, misWeight: %f", wLight, wCamera, misWeight);
		Spectrum contrib = (misWeight / (lightPickProb * directPdfW)) * (radiance * bsdfFactor);
		color = pathState.mThroughput * contrib;

		return color;
	}

	void connectToEye(SubPathState &pathState, const Intersection& its, const BSDF* bsdf, int iterationNum) {
		// Sample point on a sensor
		DirectSamplingRecord dRec(its.p, its.time);

		PerspectiveCamera *sensor = dynamic_cast<PerspectiveCamera *>(mScene->getSensor());
		if (sensor == 0) {
			Log(EError, "Only perspective camera supported");
		}

		sensor->sampleDirect(dRec, mSampler->next2D());

		// dRec.d is direction to camera
		const Float cosAtCamera = dot(dRec.n, -dRec.d);

		if (cosAtCamera <= 0.0f) {
			return;
		}

		// Compute pixel position of current point
		DirectionSamplingRecord dirRec(-dRec.d);
		Point2 imagePos;

		if (!(sensor->getSamplePosition(dRec, dirRec, imagePos))) {
			return;
		}

		if (dot(its.toLocal(dRec.d), its.wi) < 0) {
			//return;
		}
		// Compute BSDF
		BSDFSamplingRecord bsdfRec(its, its.toLocal(dRec.d));

		const Spectrum bsdfFactor = bsdf->eval(bsdfRec);

		if (bsdfFactor.isZero()) {
			return;
		}

		Float continuationProbability = 1; //std::min(pathState.mThroughput.max(), (Float) 0.95f);

		bsdfRec.reverse();
		Float bsdfRevPdfW = bsdf->pdf(bsdfRec) * continuationProbability;

		Float cosToCamera = std::abs(Frame::cosTheta(bsdfRec.wo));


		const Float imagePlaneDist = sensor->getFilm()->getSize().x / (2.0f * std::tan(sensor->getXFov() * M_PI / 360));
		const Float imagePointToCameraDist = imagePlaneDist / cosAtCamera;
		const Float imageToSolidAngleFactor = imagePointToCameraDist * imagePointToCameraDist/ cosAtCamera;
		const Float imageToSurfaceFactor = imageToSolidAngleFactor * std::abs(cosToCamera) / (dRec.dist * dRec.dist);

		const Float cameraPdfA = imageToSurfaceFactor;
		const Float wLight = (cameraPdfA / mLightSubPathCount) * (mMisVMWeightFactor + pathState.dVCM + pathState.dVC * bsdfRevPdfW);
		const Float misWeight = 1.f / (1.f + wLight);
		const Float surfaceToImageFactor = 1.f / imageToSurfaceFactor;
		const Spectrum contrib = misWeight * pathState.mThroughput * bsdfFactor / (mLightSubPathCount * surfaceToImageFactor);
		//Log(EInfo, "throughput: %s", pathState.mThroughput.toString().c_str());
		//Log(EInfo, "misWeight: %f, surfaceToImageFactor: %f, contrib: %s", misWeight, surfaceToImageFactor, contrib.toString().c_str());

		if (!contrib.isZero()) {
			if(occluded(its.p, dRec.d, dRec.dist)) {
				return;
			}

			accumulateColor(imagePos.x, imagePos.y, contrib, iterationNum);
		}
	}

	MTS_DECLARE_CLASS()

private:

	bool SampleScattering(
			Ray &ray,
			SubPathState &pathState,
			const BSDF *bsdf,
			const Intersection &its) {

		BSDFSamplingRecord bRec(its, mSampler, EImportance);
		Float bsdfDirPdfW, bsdfRevPdfW, cosThetaOut;
		bRec.wi = its.wi;

		bsdf->sample(bRec, bsdfDirPdfW, mSampler->next2D());
		Spectrum bsdfFactor = bsdf->eval(bRec, ESolidAngle);


		if (bsdfFactor.isZero())
			return false;

		Vector direction = its.toWorld(bRec.wo);
		ray = Ray(its.p, direction, its.time);

		bsdfDirPdfW = bsdf->pdf(bRec);
		bsdfRevPdfW = bsdfDirPdfW;
		if (bsdf->getType() & BSDF::ENonSymmetric) {
			bRec.reverse();
			bsdfRevPdfW = bsdf->pdf(bRec);
		}

		// TODO: Russian Roulette
		Float contProb = std::min(pathState.mThroughput.max(), (Float) 0.95f);

		if (mSampler->next1D() > contProb) {
			return false;
		}

		//bsdfDirPdfW /= contProb;
		//bsdfRevPdfW /= contProb;
		cosThetaOut = std::abs(Frame::cosTheta(bRec.wo));

		if (bsdf->getType() & BSDF::EDelta) {
			Assert(bsdfDirPdfW == bsdfRevPdfW);

			pathState.dVCM = 0.f;
			pathState.dVC *= cosThetaOut;
			pathState.dVM *= cosThetaOut;
			pathState.mSpecularPath = true;
		} else {
			pathState.dVC = (cosThetaOut / bsdfDirPdfW) * (
					pathState.dVC * bsdfRevPdfW +
					pathState.dVCM + mMisVMWeightFactor);

			pathState.dVM = (cosThetaOut / bsdfDirPdfW) * (
					pathState.dVM * bsdfRevPdfW +
					pathState.dVCM * mMisVCWeightFactor + 1.f);

			pathState.dVCM = 1.f / bsdfDirPdfW;
			pathState.mSpecularPath = false;
		}

		pathState.mThroughput *= bsdfFactor / bsdfDirPdfW;
		//Log(EInfo, "After scattering dVCM: %f, dVC: %f\n", pathState.dVCM, pathState.dVC);
		return true;
	}

	ref <ParallelProcess> m_process;
	std::vector <VCMVertex> m_lightVertices;
	std::vector<int> m_pathEnds;
	PointKDTree <VCMTreeEntry> m_tree;
	VCMConfiguration m_config;
	MemoryPool m_pool;
	ref<Bitmap> mBitmap;
	std::vector<Bitmap*> mBitmaps;
	Scene *mScene;
	Sampler *mSampler;

	Float mEtaVCM;
	Float mMisVMWeightFactor;
	Float mMisVCWeightFactor;
	Float mLightSubPathCount;
};

MTS_IMPLEMENT_CLASS_S(VCMIntegrator,false, Integrator)
MTS_EXPORT_PLUGIN(VCMIntegrator,"Vertex Connection And Merging Integrator");
MTS_NAMESPACE_END
