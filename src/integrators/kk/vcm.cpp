#include <mitsuba/bidir/vertex.h>
#include <mitsuba/core/kdtree.h>
#include <mitsuba/bidir/edge.h>
#include <mitsuba/bidir/path.h>
#include <vector>
#include <utility>

#include "vcm.h"


MTS_NAMESPACE_BEGIN

typedef PathVertex
*
PathVertexPtr;
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

	VCMVertex(const Point &p, const Spectrum &throughput,
			uint32_t pathLength, float dvcm, float dvc, float dvm,
			const BSDF *bsdf) :
			mPosition(p), mThroughput(throughput),
			mPathLength(pathLength), dVCM(dvcm),
			dVC(dvc), dVM(dvm), mBsdf(bsdf) {
	}

	Point mPosition;
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
		m_config.maxDepth = props.getInteger("maxDepth", -1);
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

		const Float radius = 1e-4;
		const Float radiusSqr = radius * radius;

		scene->initializeBidirectional();
		ref <Sensor> sensor = scene->getSensor();
		Film *film = sensor->getFilm();
		const Vector2i res = film->getSize();
		scene->getSampler()->advance();
		int pathCount = res.x * res.y;

		mBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getSize());
		mBitmap->clear();

		mEtaVCM = (M_PI * radiusSqr) * pathCount;
		mMisVMWeightFactor = mEtaVCM;
		mMisVCWeightFactor = 1.f / mEtaVCM;


		//////////////////////////////////////////////////////////////////////////
		// Generate light paths
		//////////////////////////////////////////////////////////////////////////

		m_lightVertices.reserve(pathCount);
		m_lightVertices.clear();
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
					Float lightPickProb = LightPickProbability(scene);

					pathState.mIsFiniteLight = !(emitter->getType() & AbstractEmitter::EDeltaDirection);

					Vector wo = emitterPath->vertex(2)->getPosition() - pRec.p;
					Float dist = wo.length();
					wo /= dist;
					Float cosTheta = std::abs(dot(pRec.n, wo));

					DirectionSamplingRecord dRec(wo, ESolidAngle);

					// I think it's area pdf here and should be
					// solidAngle...
					Float directPdf = emitter->pdfPosition(pRec) * lightPickProb;
					Float emissionPdf = emitter->pdfDirection(dRec, pRec) * lightPickProb;
					//Log(EInfo, "Direct= %f, emission= %f", directPdf, emissionPdf);

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
					DirectSamplingRecord dRec(its);
					const BSDF *bsdf = its.getBSDF();
					// computing SampleScattering PDF values here since we
					// already have a computed path
					SampleScattering(pathState, bsdf, dRec, its, emitterPath, vertexIdx);

					// is this the right cos angle to compute?
					Float cosTheta = std::abs(dot(its.toWorld(-its.wi), its.geoFrame.n));
					if (vertexIdx > 2 || pathState.mIsFiniteLight) {
						pathState.dVCM *= its.t * its.t;
					}

					pathState.dVCM /= cosTheta;
					pathState.dVC /= cosTheta;
					pathState.dVM /= cosTheta;


					// add vertex iff the bsdf is not purely specular
					// not sure if thats how you check it ;)
					if (!(bsdf->getType() & BSDF::EDelta)) {
						VCMVertex v(vertex->getPosition(),
								pathState.mThroughput,
								vertexIdx - 1,
								pathState.dVCM, pathState.dVC, pathState.dVM, its.getBSDF());
						m_lightVertices.push_back(v);
					}
					if (!(bsdf->getType() & BSDF::EDelta)) {
						connectToEye(scene, time, emitterPath, pathState, vertexIdx);
					}
				}

			}

			scene->getSampler()->advance();
		}

		Log(EInfo, "Starting to build...");

		/////////////////////////////////////////////////////////////////////////
		// BUILD SEARCH STRUCT
		/////////////////////////////////////////////////////////////////////////

		m_tree.reserve(pathCount);
		for (int i = 0; i < m_lightVertices.size(); i++) {
			m_tree.push_back(m_lightVertices[i]);
		}

		Log(EInfo, "Built the tree");


		/////////////////////////////////////////////////////////////////////////
		// GENERATE CAMERA PATHS
		/////////////////////////////////////////////////////////////////////////

		SubPathState pathState;
		Spectrum *target = (Spectrum *) mBitmap->getUInt8Data();
		for (int pathIdx = 0; pathIdx < pathCount; ++pathIdx) {
			// Generate eye path
			Path *sensorPath = new Path();
			sensorPath->initialize(scene, 0, ERadiance, m_pool);
			int currentX = pathIdx % res.x;
			int currentY = pathIdx / res.x;
			Point2i startPosition = Point2i(currentX, currentY);
			sensorPath->randomWalkFromPixel(scene, scene->getSampler(), m_config.maxDepth, startPosition, m_config.rrDepth, m_pool);

			// watch out for the "-1" part. I don't count the last vertex
			// since it doesn't have w_0 vector
			for (int vertexIdx = 1; vertexIdx < sensorPath->vertexCount() - 1; vertexIdx++) {
				PathVertexPtr vertex = sensorPath->vertex(vertexIdx);

				if (vertexIdx == 1 && sensorPath->vertexCount() > 2) {
					PathVertexPtr nextVertex = sensorPath->vertex(2);
					if (!(nextVertex->getType() & PathVertex::ESurfaceInteraction)) {
						Log(EInfo, "Other vertex type: %x", vertex->getType());
						break;
					}
					PositionSamplingRecord pRec = vertex->getPositionSamplingRecord();

					Assert(vertex->type & PathVertex::ESensorSample);

					Vector wo = nextVertex->getPosition() - pRec.p;
					Float dist = wo.length();
					wo /= dist;

					DirectionSamplingRecord dRec(wo, ESolidAngle);

					Float cameraPdfW = dRec.pdf;

					pathState.mSpecularPath = true;
					pathState.mThroughput = Spectrum(1.0f);
					pathState.dVCM = pathCount / cameraPdfW;
					pathState.dVC = pathState.dVM = 0;
				} else {
					if (!(vertex->getType() & PathVertex::ESurfaceInteraction ||
							vertex->isEmitterSample())) {
						Log(EInfo, "Other vertex type and not emitter: %d", vertex->getType());
						break;
					}

					const Intersection &its = vertex->getIntersection();
					DirectSamplingRecord dRec(its);
					const BSDF *bsdf = its.getBSDF();
					SampleScattering(pathState, bsdf, dRec, its, sensorPath, vertexIdx);

					Float cosTheta = std::abs(dot(its.toWorld(-its.wi), its.geoFrame.n));

					pathState.dVCM *= its.t * its.t;
					pathState.dVCM /= cosTheta;
					pathState.dVC /= cosTheta;
					pathState.dVM /= cosTheta;

					Spectrum color(0.0f);
					// Case #1 - we hit the light
					if (vertex->isEmitterSample()) {
						PositionSamplingRecord pRec = vertex->getPositionSamplingRecord();
						const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
						Float lightPickProb = LightPickProbability(scene);

						Vector wi = pRec.p - sensorPath->vertex(vertexIdx - 1)->getPosition();
						Float dist = wi.length();
						wi /= dist;

						DirectionSamplingRecord dRec(-wi, ESolidAngle);

						Float directPdf = emitter->pdfPosition(pRec);
						Float emissionPdf = emitter->pdfDirection(dRec, pRec);
						// WARNING: There be dragons
						Spectrum radiance = GetLightRadiance(emitter, its, -wi);
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
						queue->signalRefresh(job);
						Log(EInfo, "Woot, we hit the light!!!");
						Log(EInfo, "Adding color %s on position (%d, %d)", color.toString().c_str(), currentX, currentY);
					}
					// Case #2 - Vertex Connection - connect to light
					// source
					if (!(bsdf->getType() & BSDF::EDelta)) {
						Sampler *sampler = scene->getSampler();
						ref_vector <Emitter> &emitters = scene->getEmitters();
						ref <Emitter> sampledEmitter = emitters[emitters.size() * sampler->next1D()];
						Float lightPickProb = LightPickProbability(scene);

						Point2 randomPoint = sampler->next2D();
						PositionSamplingRecord emitterPRec;
						sampledEmitter->samplePosition(emitterPRec, randomPoint);
						Point emitterPoint = emitterPRec.p;
						Vector dirToLight = emitterPoint - its.p;
						Float dist = dirToLight.length();
						dirToLight /= dist;

						Spectrum radiance = GetLightRadiance(sampledEmitter, its, dirToLight);
						//Log(EInfo, "Radiance %s", radiance.toString().c_str());

						if (radiance != Spectrum(0.0f)) {
							Float cosNormalDir = dot(emitterPRec.n, -dirToLight);

							Float directPdfW = emitterPRec.pdf * dist * dist / cosNormalDir;
							Float emissionPdfW = emitterPRec.pdf * cosNormalDir * INV_PI;

							PathVertexPtr nextVertex = sensorPath->vertex(vertexIdx + 1);
							Float bsdfRevPdfW, bsdfDirPdfW, cosThetaOut;
							BSDFSamplingRecord bsdfRec(its, its.toLocal(dirToLight));
							cosThetaOut = Frame::cosTheta(bsdfRec.wo);

							Spectrum bsdfFactor = bsdf->eval(bsdfRec);
							if (bsdfFactor == Spectrum(0.0f)) {
								//Log(EInfo, "bsdf factor 0, breaking");
								continue;
							}
							bsdfDirPdfW = bsdf->pdf(bsdfRec);

							bsdfRec.reverse();
							bsdfRevPdfW = bsdf->pdf(bsdfRec);
							Float continuationProbability = 1.0f / vertex->rrWeight;

							if (sampledEmitter->isDegenerate()) {
								bsdfDirPdfW = 0.0f;
							}

							bsdfDirPdfW *= continuationProbability;
							bsdfRevPdfW *= continuationProbability;
							Float wLight = bsdfDirPdfW / (lightPickProb * directPdfW);
							Float wCamera = (emissionPdfW * cosThetaOut / (directPdfW * cosNormalDir)
									* (mMisVMWeightFactor + pathState.dVCM + pathState.dVC * bsdfRevPdfW));
							//Log(EInfo, "\n  emissionPdfW: %f\n  cosThetaOut: %f\n  directPdfW: %f\n  cosNormalDir: %f\n  mMisVMWeightFactor: %f\n  dVCM: %f\n  dVC: %f\n  bsdfRevPdfW: %f\n", emissionPdfW, cosThetaOut, directPdfW, cosNormalDir, mMisVMWeightFactor, pathState.dVCM, pathState.dVC, bsdfRevPdfW);
							Float misWeight = 1.f / (wLight + 1.f + wCamera);
							Spectrum contrib = (misWeight * cosThetaOut / (lightPickProb * directPdfW)) * (radiance * bsdfFactor);
							color = pathState.mThroughput * contrib;
							target[currentY * mBitmap->getWidth() + currentX] += color;
							//Log(EInfo, "\n  wLight: %f\n  wCamera: %f\n  misWeight: %f\n  cosThetaOut: %f\n  lightPickProb: %f\n  directPdfW: %f\n  radiance: %s\n  bsdfFactor: %s\n", wLight, wCamera, misWeight, cosThetaOut, lightPickProb, directPdfW, radiance.toString().c_str(), bsdfFactor.toString().c_str());
							//Log(EInfo, "bsdfDirPdfW: %f\n bsdfRevPdfW: %f\n contProb: %f\n lightPickProb: %f\n, directPdfW: %f\n emissionPdfW: %f\n cosThetaOut: %f\n cosNormalDir: %f\n", bsdfDirPdfW, bsdfRevPdfW, continuationProbability, lightPickProb, directPdfW, emissionPdfW, cosThetaOut, cosNormalDir);
							//Log(EInfo, "Adding color %s, contrib: %s, througput: %s on position (%d, %d)", color.toString().c_str(), contrib.toString().c_str(), pathState.mThroughput.toString().c_str(), currentX, currentY);
							queue->signalRefresh(job);
						}

					}
				}

			}
			scene->getSampler()->advance();
		}


		Log(EInfo, "DONE!");
		film->setBitmap(mBitmap);
		queue->signalRefresh(job);

		return true;
	}

	Float LightPickProbability(const Scene *scene) {
		return 1.0f / scene->getEmitters().size();
	}

	Spectrum GetLightRadiance(const Emitter *emitter, const Intersection &its, const Vector &d) {
		Spectrum radiance(0.0f);
		try {
			radiance = emitter->eval(its, d);
		} catch (...) {
		}
		return radiance;
	}

	bool connectToEye(Scene *scene, Float time, Path *emitterPath, SubPathState &pathState, int vertexIdx) {
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
		DirectSamplingRecord dRec(its);
		Vector wo = vertexOnEye->getPosition() - dRec.p;
		Float dist = wo.length();
		wo /= dist;
		BSDFSamplingRecord bsdfRec(its, wo);
		const Spectrum bsdfFactor = bsdf->eval(bsdfRec);
		if (bsdfFactor.isZero())
			return false;

		Float continuationProbability = std::min(pathState.mThroughput.max(), (Float) 0.95f);
		bsdfRec.reverse();
		Float bsdfRevPdfW = bsdf->pdf(bsdfRec) * continuationProbability;

		emitterPathCopy->release(m_pool);
		delete sensorPath;
		delete emitterPathCopy;
		scene->getSampler()->advance();
		return succeded;
	}

	MTS_DECLARE_CLASS()

private:

	void SampleScattering(
			SubPathState &pathState,
			const BSDF *bsdf,
			const DirectSamplingRecord &dRec,
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
			cosThetaOut = Frame::cosTheta(bsdfRec.wo);
			bsdfDirPdfW = bsdf->pdf(bsdfRec);

			// same for specular
			if (!(bsdf->getType() & BSDF::EDiffuse)) {
				bsdfRevPdfW = bsdfDirPdfW;
			} else { // differs for non-specular
				bsdfRec.reverse();
				bsdfRevPdfW = bsdf->pdf(bsdfRec);
			}

			// TODO: Russian roulette factor

			if (!(bsdf->getType() & BSDF::EDiffuse)) {
				pathState.dVCM = 0.f;
				assert(bsdfDirPdfW == bsdfRevPdfW);
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

			Log(EInfo, "cosThetaOut: %f, bsdfDirPdfW: %f, bsdf: %s",
					cosThetaOut, bsdfDirPdfW, bsdf->eval(bsdfRec).toString().c_str());
			pathState.mThroughput *= bsdf->eval(bsdfRec) * (cosThetaOut / bsdfDirPdfW);
			//Log(EInfo, "Throughput: %s", pathState.mThroughput.toString().c_str());
		}
	}

	ref <ParallelProcess> m_process;
	std::vector <VCMVertex> m_lightVertices;
	PointKDTree <VCMTreeEntry> m_tree;
	VCMConfiguration m_config;
	MemoryPool m_pool;
	ref <Bitmap> mBitmap;

	Float mEtaVCM;
	Float mMisVMWeightFactor;
	Float mMisVCWeightFactor;
};

MTS_IMPLEMENT_CLASS_S(VCMIntegrator,
false, Integrator)
MTS_EXPORT_PLUGIN(VCMIntegrator,
"Vertex Connection And Merging Integrator");
MTS_NAMESPACE_END
