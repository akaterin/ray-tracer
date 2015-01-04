#include <mitsuba/bidir/vertex.h>
#include <mitsuba/core/kdtree.h>
#include <mitsuba/bidir/edge.h>
#include <mitsuba/bidir/path.h>
#include <vector>
#include <utility>
#include <cmath>


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
		mSampler->advance();
		int pathCount = res.x * res.y;
		mLightSubPathCount = float(res.x * res.y);

		mBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getSize());
		mBitmap->clear();


		const Float radiusAlpha = 0.75f;

		for (int iterationNum = 0; iterationNum < 1; iterationNum++) {

			Float radius = scene->getBSphere().radius * 0.0003f;
			radius /= std::pow(Float(iterationNum + 1), 0.5f * (1 - radiusAlpha));
			const Float radiusSqr = radius * radius;

			Log(EInfo, "Scene BSphere radius: %f", scene->getBSphere().radius);
			Log(EInfo, "VCM Radius %f", radius);

			mEtaVCM = (M_PI * radiusSqr) * pathCount;
			mMisVMWeightFactor = 0.0f; //mEtaVCM;
			mMisVCWeightFactor = 1.f / mEtaVCM;


			Log(EInfo, "pathCount: %d", pathCount);
			Log(EInfo, "etaVCM: %f\nVM Weight: %f\nVC Weight: %f\n", mEtaVCM, mMisVMWeightFactor, mMisVCWeightFactor);
			//////////////////////////////////////////////////////////////////////////
			// Generate light paths
			//////////////////////////////////////////////////////////////////////////

			/*m_lightVertices.reserve(pathCount);
			m_pathEnds.reserve(pathCount);
			m_lightVertices.clear();
			m_pathEnds.clear();
			for (int i = 0; i < pathCount; ++i) {
				Float time = i * 1000;
				Path *emitterPath = new Path();
				emitterPath->initialize(scene, time, EImportance, m_pool);
				emitterPath->randomWalk(scene, scene->getSampler(), m_config.maxDepth, m_config.rrDepth, EImportance, m_pool);
				SubPathState pathState;

				// skip Emitter Supernode
				for (int vertexIdx = 1; vertexIdx < emitterPath->vertexCount(); vertexIdx++) {
					PathVertexPtr vertex = emitterPath->vertex(vertexIdx);
					if (!(vertex->type & PathVertex::ENormal)) {
						continue;
					}

				// if it's on emitter, sample MIS values
				if (vertexIdx == 1 && emitterPath->vertexCount() > 2) {
					PathVertexPtr nextVertex = emitterPath->vertex(2);
					if (!(nextVertex->getType() & PathVertex::ESurfaceInteraction)) {
						break;
					}

						PositionSamplingRecord pRec = vertex->getPositionSamplingRecord();

						Assert(vertex->type & PathVertex::EEmitterSample);
						//PathEdge* path = emitterPath->edge(1);

						//Float emissionPdf = path->pdf[ERadiance];
						//Float directPdf = pRec.pdf; //scene->pdfEmitterDirect(dRec);
						const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
						Float lightPickProb = LightPickProbability();

						pathState.mIsFiniteLight = !(emitter->getType() & AbstractEmitter::EDeltaDirection);

						Vector wo = emitterPath->vertex(2)->getPosition() - pRec.p;
						Float dist = wo.length();
						wo /= dist;
						Float cosTheta = std::abs(dot(pRec.n, wo));

						DirectionSamplingRecord dRec(wo, ESolidAngle);

						Float directPdf = emitter->pdfPosition(pRec) * lightPickProb;
						Float emissionPdf = directPdf * emitter->pdfDirection(dRec, pRec) * lightPickProb;

						pathState.mThroughput = vertex->weight[ERadiance] / emissionPdf;

						pathState.dVCM = directPdf / emissionPdf;

						// TODO: handle delta and infinite lights
						pathState.dVC = cosTheta / emissionPdf;
						pathState.dVM = pathState.dVC * mMisVCWeightFactor;

						//Log(EInfo, "\n directPdf: %f\n emissionPdf: %f\n cosTheta: %f\n",
						//	directPdf, emissionPdf, cosTheta);
					}
					else {
						if (!(vertex->getType() & PathVertex::ESurfaceInteraction)) {
							//    Log(EInfo, "Other vertex type: %d", vertex->getType());
							break;
						}

						const Intersection &its = vertex->getIntersection();
						const BSDF *bsdf = its.getBSDF();

						Float cosTheta = std::abs(Frame::cosTheta(-its.wi));
						if (vertexIdx > 2 || pathState.mIsFiniteLight) {
							pathState.dVCM *= its.t * its.t;
						}

						pathState.dVCM /= cosTheta;
						pathState.dVC /= cosTheta;
						pathState.dVM /= cosTheta;

						if (!(bsdf->getType() & BSDF::EDelta)) {
							VCMVertex v(vertex->getPosition(),
									its,
									pathState.mThroughput,
									vertexIdx - 1,
									pathState.dVCM, pathState.dVC, pathState.dVM, its.getBSDF());
							m_lightVertices.push_back(v);
						}
						if (!(bsdf->getType() & BSDF::EDelta)) {
							//connectToEye(scene, time, emitterPath, pathState, vertexIdx, iterationNum);
						}

						SampleScattering(pathState, bsdf, its, emitterPath, vertexIdx);
					}

				}

				m_pathEnds[i] = (int) m_lightVertices.size();
				scene->getSampler()->advance();
			}

			Log(EInfo, "Starting to build...");
			/////////////////////////////////////////////////////////////////////////
			// BUILD SEARCH STRUCT
			/////////////////////////////////////////////////////////////////////////

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
			for (int pathIdx = 0; pathIdx < pathCount; ++pathIdx) {
				// Generate eye path
				Path *sensorPath = new Path();
				//sensorPath->initialize(scene, 0, ERadiance, m_pool);
				int currentX = pathIdx % res.x;
				int currentY = pathIdx / res.x;
				Point2i startPosition = Point2i(currentX, currentY);
				Ray ray = generateCameraSample(pathState, pathIdx, pathCount);

//				sensorPath->randomWalkFromPixel(scene, scene->getSampler(), m_config.maxDepth, startPosition, m_config.rrDepth, m_pool);

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

					//Log(EInfo, "Before scattering d: %f, dVCM: %f, dVC: %f\n", its.t, pathState.dVCM, pathState.dVC);
					Spectrum color(0.0f);
					// Case #1 - we hit the light
					if (its.isEmitter()) {
						Spectrum radiance = its.Le(its.toWorld(-its.wi));
						// we see light directly from camera
						// supernode->sensor->emitter
						if (pathState.mPathLength == 1) {
							color += Spectrum(1.0f);
							accumulateColor(currentX, currentY, color, iterationNum);
							//Log(EInfo, "Radiance emitter in opposite direction %s: %s", (its.wi).toString().c_str(),
							//	radiance.toString().c_str());
						} else {
							//color += pathState.mThroughput * radiance;
							//Log(EInfo, "Radiance emitter in direction %s: %s", (-its.wi).toString().c_str(),
							//		radiance.toString().c_str());
							// TO BE IMPLEMENTED
						}
						break;
					}

					if (m_config.maxDepth != -1 && pathState.mPathLength >= m_config.maxDepth) {
						break;
					}
					if (false && (its.isValid() && its.isEmitter())) {
						/*
						   PositionSamplingRecord pRec = vertex->getPositionSamplingRecord();
						//Log(EInfo, "Built pRec");
						//const Emitter *emitter = static_cast<const Emitter *>(dRec.object);
						//Log(EInfo, "Got emitter handle");
						//Log(EInfo, "%d", dRec.object);
						//Log(EInfo, dRec.object->getClass()->getName().c_str());
						Float lightPickProb = LightPickProbability();

						Vector wi = pRec.p - sensorPath->vertex(vertexIdx - 1)->getPosition();
						Float dist = wi.length();
						wi /= dist;

						DirectionSamplingRecord dRec(-wi, ESolidAngle);
						//Log(EInfo, "Built dRec");

						//Float directPdf = emitter->pdfPosition(pRec);
						//Float emissionPdf = emitter->pdfDirection(dRec, pRec);
						//Log(EInfo, "Computed pdfs");
						// WARNING: There be dragons
						//Spectrum radiance = GetLightRadiance(emitter, its, -wi);
						// we see light directly from camera
						// supernode->sensor->emitter
						if (vertexIdx == 2) {
						color = radiance;
						} else {
						directPdf *= lightPickProb;
						emissionPdf *= lightPickProb;

						const Float wSensor = directPdf * pathState.dVCM +
						emissionPdf * pathState.dVC;

						const Float misWeight = 1.f / (1.f + wSensor);

						color = pathState.mThroughput * misWeight * radiance;
						}
						target[currentY * mBitmap->getWidth() + currentX] += color;
						film->setBitmap(mBitmap);
						queue->signalRefresh(job);
						Log(EInfo, "Woot, we hit the light!!!");
						Log(EInfo, "Adding color %s on position (%d, %d)", color.toString().c_str(), currentX, currentY);*/
					}

					// TODO: terminate if sub-path too long for
					// connection/merging

					// Case #2 - Vertex Connection - connect to light
					// source
					if (!(bsdf->getType() & BSDF::EDelta)) {
						color += ConnectToLight(pathState, its, bsdf);
					} else {
						Log(EInfo, "Delta BSDF, cant do VC");
					}

					// Case #3 - Vertex Connection - Connect to light vertices
					if (false && !(bsdf->getType() & BSDF::EDelta)) {
						int startingLightVertex = pathIdx == 0 ? 0 : m_pathEnds[pathIdx-1];
						int lastLightVertex = m_pathEnds[pathIdx];

						for (int i = startingLightVertex; i < lastLightVertex; i++) {
							const VCMVertex &lightVertex = m_lightVertices[i];

							// TODO: check min path length
							//							if (lightVertex.mPathLength + vertexIdx
							//									< mMinPathLength)
							//								continue;

							// TODO: check max path length
							//if (lightVertex.mPathLength + vertexIdx > mMaxPathLength)
							//	break;

							//color += pathState.mThroughput * lightVertex.mThroughput *
							//	ConnectVertices(pathState, vertex, lightVertex, its, bsdf);
						}
					}

					accumulateColor(currentX, currentY, color, iterationNum);
					film->setBitmap(mBitmap);
					queue->signalRefresh(job);
					mSampler->advance();
					if(!SampleScattering2(ray, pathState, bsdf, its)) {
						break;
					}
				}
			}
			Log(EInfo, "Done iteration %d", iterationNum + 1);
			Log(EInfo, "0 radiance times: %d", zeroRadianceCount);
		}

		Log(EInfo, "DONE!");
		return true;
	}

	Ray generateCameraSample(SubPathState& pathState, int pixelIndex, int pathCount) {
		Ray ray;
		Sensor *sensor = mScene->getSensor();
		const Vector2i res = sensor->getFilm()->getSize();
		int currentX = pixelIndex % res.x;
		int currentY = pixelIndex / res.x;
		Point2 startPosition = Point2(currentX, currentY);
		sensor->sampleRay(
				ray,
				startPosition,
				startPosition, // assume we don't need aperture sampling
				0			   // assume we don't need time sampling
				);

		PositionSamplingRecord pRec;
		sensor->samplePosition(pRec, mSampler->next2D());
		mSampler->advance();
		DirectionSamplingRecord dRec(ray.d, ESolidAngle);
		Float cameraPdfW = sensor->pdfDirection(dRec, pRec);

		pathState.mPathLength = 1;
		pathState.mSpecularPath = true;
		pathState.mThroughput = Spectrum(1.0f);
		pathState.dVCM = pathCount / cameraPdfW;
		pathState.dVC = pathState.dVM = 0;

		return ray;
	}

	void accumulateColor(int x, int y, const Spectrum& color, int iteration) {
		Spectrum *target = (Spectrum *) mBitmap->getUInt8Data();
		Float iterationFactor = (iteration == 0 ? 0 : (Float)iteration / (iteration + 1));
		Float colorIterationFactor = (iteration == 0 ? 1 : (Float)1.0f / (iteration+1));
		target[y * mBitmap->getWidth() + x] =
			  target[y * mBitmap->getWidth() + x] * iterationFactor
			+ color * colorIterationFactor;
	}

	Spectrum ConnectVertices(SubPathState &pathState, PathVertexPtr vertex,
			const VCMVertex &lightVertex, const Intersection &its, const BSDF *bsdf) {

		Vector direction = lightVertex.mPosition - its.p;
		Float dist = direction.length();
		direction /= dist;

		Float cosCamera, cameraBsdfDirPdfW, cameraBsdfRevPdfW;
		BSDFSamplingRecord bsdfRec(its, its.toLocal(direction));
		Spectrum cameraBsdfFactor = bsdf->eval(bsdfRec);
		cosCamera = std::abs(Frame::cosTheta(bsdfRec.wo));

		if (cameraBsdfFactor.isZero()) {
			return Spectrum(0.0f);
		}

		cameraBsdfDirPdfW = bsdf->pdf(bsdfRec);
		bsdfRec.reverse();
		cameraBsdfRevPdfW = bsdf->pdf(bsdfRec);

		Float continuationProbability = 1.0f / vertex->rrWeight;
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
		Float lightCont = 1.0f / vertex->rrWeight;
		//lightBsdfDirPdfW *= lightCont;
		//lightBsdfRevPdfW *= lightCont;

		Float geometryTerm = 1.0f / (dist * dist);

		// Not sure how that's possible
		if (geometryTerm < 0) {
			return Spectrum(0.0f);
		}

		Float cameraBsdfDirPdfA = pdfWtoA(cameraBsdfDirPdfW, dist, cosLight);
		Float lightBsdfDirPdfA = pdfWtoA(lightBsdfDirPdfW, dist, cosCamera);

		Float wLight = cameraBsdfDirPdfA * (mMisVMWeightFactor + lightVertex.dVCM
				+ lightVertex.dVC * lightBsdfRevPdfW);

		Float wCamera = lightBsdfDirPdfA * (mMisVMWeightFactor + pathState.dVCM
				+ pathState.dVC * cameraBsdfRevPdfW);

		Float misWeight = 1.0f / (wLight + 1.0f + wCamera);

		Spectrum contrib = (misWeight * geometryTerm) * cameraBsdfFactor * lightBsdfFactor;

		if (contrib.isZero() || occluded(its.p, direction, dist)) {
			return Spectrum(0.0f);
		}

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

	int zeroRadianceCount = 0;
	Spectrum ConnectToLight(SubPathState &pathState, const Intersection &its, const BSDF *bsdf) {
		Float lightPickProb = LightPickProbability();

		Point2 randomPoint = mSampler->next2D();
		mSampler->advance();

		Spectrum color(0.0f), radiance(0.0f);
		DirectSamplingRecord dRec(its);
		Assert(dRec.ref == its.p);

		radiance = mScene->sampleEmitterDirect(dRec, randomPoint, true);

		if (radiance.isZero()) {
			//Log(EInfo, "dRec = %s", dRec.toString().c_str());
			//Log(EInfo, "Zero radiance, data: \n\n%s\n\n%s\n\n", its.toString().c_str(), dRec.toString().c_str());
		}
		radiance *= dRec.pdf;

		if (radiance.isZero()) {
			++zeroRadianceCount;
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
		Spectrum contrib = (/*misWeight*/ 1.f / (lightPickProb * directPdfW)) * (radiance * bsdfFactor);
		color = pathState.mThroughput * contrib;

		return color;
	}

	Spectrum GetLightRadiance(const Emitter *emitter, const Intersection &its, const Vector &d, Float *cosAtLight = NULL) {
		Spectrum radiance(0.0f);
		try {
			Ray r(its.p, -d, 0);
			Intersection lightIts;
			mScene->rayIntersect(r, lightIts);
			if (lightIts.isValid() && !lightIts.isEmitter()) {
				//Log(EInfo, "Point %s is not an emitter", lightIts.p.toString().c_str());
			}
			if (lightIts.isValid() && lightIts.isEmitter()) {
				radiance = emitter->eval(lightIts, d);
				if (cosAtLight) {
					*cosAtLight = dot(lightIts.shFrame.n, d);
				}
			}
		} catch (...) {
		}
		return radiance;
	}

	bool connectToEye(Scene *scene, Float time, Path *emitterPath, SubPathState &pathState, int vertexIdx, int iterationNum) {
		int prevSize = emitterPath->vertexCount();
		Path *sensorPath = new Path();
		Path *emitterPathCopy = new Path();
		emitterPath->clone(*emitterPathCopy, m_pool);
		sensorPath->initialize(scene, time, ERadiance, m_pool);
		sensorPath->randomWalk(scene, scene->getSampler(), 1, 0, ERadiance, m_pool);

		PathVertex *vertexOnEye = sensorPath->vertex(1);
		PathVertex *succOnSensor = sensorPath->vertex(0);
		PathEdge *succEdge = sensorPath->edge(0);
		PathEdge *lastEdge = emitterPathCopy->edge(vertexIdx - 1);
		PathVertex *predLastVertex = emitterPathCopy->vertexOrNull(vertexIdx - 1);
		PathVertex *lastVertex = emitterPathCopy->vertex(vertexIdx);
		PathEdge *newEgde = new PathEdge();

		bool succeded = PathVertex::connect(scene, predLastVertex, lastEdge, lastVertex, newEgde, vertexOnEye, succEdge, succOnSensor);
		if (!succeded) {
			emitterPathCopy->release(m_pool);
			delete sensorPath;
			delete emitterPathCopy;
			scene->getSampler()->advance();
			return false;
		}

		emitterPathCopy->append(newEgde, vertexOnEye);
		emitterPathCopy->append(succEdge, succOnSensor);

		const Intersection &its = lastVertex->getIntersection();
		const BSDF *bsdf = its.getBSDF();
		const DirectSamplingRecord dRec(its);
		Vector wo = vertexOnEye->getPosition() - dRec.p;
		Float dist = wo.length();
		wo /= dist;
		BSDFSamplingRecord bsdfRec(its, its.toLocal(wo));
		const Spectrum bsdfFactor = bsdf->eval(bsdfRec);
		if (bsdfFactor.isZero()) {
			return false;
		}
		//not sure if it shouldn't be bsdfRec.wo.z
		//Float cosToCamera = Frame::cosTheta(bsdfRec.wo);
		Float cosToCamera = 1;

		Float continuationProbability = std::min(pathState.mThroughput.max(), (Float) 0.95f);
		bsdfRec.reverse();
		Float bsdfRevPdfW = bsdf->pdf(bsdfRec) * continuationProbability;
		const Float cosAtCamera = std::abs(dot(wo, its.geoFrame.n));
		PerspectiveCamera *perspectiveCamera = dynamic_cast<PerspectiveCamera *>(scene->getSensor());
		if (perspectiveCamera == 0) {
			Log(EError, "Only perspective camera supported");
		}

		Transform viewTransform = perspectiveCamera->getViewTransform(0);
		Transform cameraToSampleTransform = getCaneraToSample(scene, perspectiveCamera);
		Point local = viewTransform(lastVertex->getPosition());
		Point imagePos = cameraToSampleTransform(local);
		const Vector2i res = scene->getSensor()->getFilm()->getSize();
		if(imagePos.x < 0 || imagePos.y <0 || local.z < 0) {
			return false;
		}

		imagePos.x = imagePos.x*res.x;
		imagePos.y = imagePos.y*res.y;

		const Float cameraImagePlaneDist = perspectiveCamera->getXFov();
		const Float imagePointToCameraDist = cameraImagePlaneDist / cosAtCamera;
		const Float imageToSolidAngleFactor = sqrt(imagePointToCameraDist) / cosAtCamera;
		const Float imageToSurfaceFactor = imageToSolidAngleFactor * std::abs(cosToCamera) / sqrt(dist);

		const Float cameraPdfA = imageToSurfaceFactor;
		const Float wLight = (cameraPdfA / mLightSubPathCount) * (mMisVMWeightFactor + pathState.dVCM + pathState.dVC * bsdfRevPdfW);
		const Float misWeight = 1.f / (1.f + wLight);
		const Float surfaceToImageFactor = 1.f / imageToSurfaceFactor;
		const Spectrum contrib = misWeight * pathState.mThroughput * bsdfFactor / (mLightSubPathCount * surfaceToImageFactor);
		//Log(EInfo, "thtoughput: %s", pathState.mThroughput.toString().c_str());
		//Log(EInfo, "misWeight: %f, surfaceToImageFactor: %f, bdsf: %s", misWeight, surfaceToImageFactor, contrib.toString().c_str());

		if (!contrib.isZero()) {
			//TODO is occuled

			//Log(EInfo, "x: %f, y: %f", imagePos.x, imagePos.y);
			accumulateColor(imagePos.x, imagePos.y, contrib, iterationNum);
		}

		emitterPathCopy->release(m_pool);
		delete sensorPath;
		delete emitterPathCopy;
		scene->getSampler()->advance();
		return succeded;
	}

	Transform getCaneraToSample(Scene* scene, const PerspectiveCamera* perspectiveCamera){
		const Vector2i &filmSize   = scene->getSensor()->getFilm()->getSize();
		const Vector2i &cropSize   = scene->getSensor()->getFilm()->getCropSize();
		const Point2i  &cropOffset = scene->getSensor()->getFilm()->getCropOffset();
		const Float aspect = scene->getSensor()->getAspect();
		ProjectiveCamera *projectiveCamera = dynamic_cast<ProjectiveCamera *>(scene->getSensor());
		Float nearClip = projectiveCamera->getNearClip();
		Float farClip = projectiveCamera->getFarClip();
		Float xFov = perspectiveCamera->getXFov();

		Vector2 relSize((Float) cropSize.x / (Float) filmSize.x,
				(Float) cropSize.y / (Float) filmSize.y);
		Point2 relOffset((Float) cropOffset.x / (Float) filmSize.x,
				(Float) cropOffset.y / (Float) filmSize.y);

		Transform cameraToSample =
				Transform::scale(Vector(1.0f / relSize.x, 1.0f / relSize.y, 1.0f))
						* Transform::translate(Vector(-relOffset.x, -relOffset.y, 0.0f))
						* Transform::scale(Vector(-0.5f, -0.5f*aspect, 1.0f))
						* Transform::translate(Vector(-1.0f, -1.0f/aspect, 0.0f))
						* Transform::perspective(xFov, nearClip, farClip);

		return cameraToSample;
 	}

	MTS_DECLARE_CLASS()

private:

	bool SampleScattering2(
			Ray &ray,
			SubPathState &pathState,
			const BSDF *bsdf,
			const Intersection &its) {

		BSDFSamplingRecord bRec(its, mSampler, EImportance);
		Float bsdfDirPdfW, bsdfRevPdfW, cosThetaOut;
		bRec.wi = its.wi;

		bsdf->sample(bRec, bsdfDirPdfW, mSampler->next2D());
		Spectrum bsdfFactor = bsdf->eval(bRec, ESolidAngle);

		mSampler->advance();

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
			Log(EInfo, "Monte Carlo break, contProb: %f", contProb);
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

	void SampleScattering(
			SubPathState &pathState,
			const BSDF *bsdf,
			const Intersection &its,
			const Path *path,
			int vertexIdx) {

		PathVertexPtr vertex = path->vertex(vertexIdx);
		if (path->vertexCount() > vertexIdx + 1) {
			PathVertexPtr nextVertex = path->vertex(vertexIdx + 1);
			Float bsdfRevPdfW, bsdfDirPdfW, cosThetaOut;
			Vector wo = nextVertex->getPosition() - vertex->getPosition();
			Float dist = wo.length();
			wo /= dist;
			BSDFSamplingRecord bsdfRec(its, its.toLocal(wo));
			cosThetaOut = std::abs(Frame::cosTheta(bsdfRec.wo));
			bsdfDirPdfW = bsdf->pdf(bsdfRec);

			// NOTE: It is already multiplied by cosThetaOut!
			Spectrum bsdfFactor = bsdf->eval(bsdfRec);

			if (bsdfFactor.isZero()) {
				Log(EInfo, "WHOOPSIE, BSDF IS 0");
				return;
			}

			// same for specular
			if (!(bsdf->getType() & BSDF::EDiffuse)) {
				bsdfRevPdfW = bsdfDirPdfW;
			} else { // differs for non-specular
				bsdfRec.reverse();
				bsdfRevPdfW = bsdf->pdf(bsdfRec);
			}

			// TODO: Russian roulette factor

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

			//Log(EInfo, "cosThetaOut: %f, bsdfDirPdfW: %f, bsdf: %s",
			//		cosThetaOut, bsdfDirPdfW, bsdf->eval(bsdfRec).toString().c_str());
			pathState.mThroughput *= bsdfFactor / bsdfDirPdfW;
			//Log(EInfo, "Throughput: %s", pathState.mThroughput.toString().c_str());
		}
	}

	ref <ParallelProcess> m_process;
	std::vector <VCMVertex> m_lightVertices;
	std::vector<int> m_pathEnds;
	PointKDTree <VCMTreeEntry> m_tree;
	VCMConfiguration m_config;
	MemoryPool m_pool;
	ref <Bitmap> mBitmap;
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
