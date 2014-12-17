
#include <mitsuba/bidir/vertex.h>
#include <mitsuba/core/kdtree.h>
#include <mitsuba/bidir/edge.h>
#include <mitsuba/bidir/path.h>
#include <vector>
#include <utility>
#include "vcm.h"

MTS_NAMESPACE_BEGIN

typedef PathVertex*		  PathVertexPtr;
//    typedef VCMKDTree::IndexType    IndexType;
//    typedef VCMKDTree::SearchResult SearchResult;

enum EVCMVertexData {
    EdVCMData = 0,
    EdVCData = 1,
    EdVMData = 2
};

struct MTS_EXPORT_RENDER VCMTreeEntry :
	public SimpleKDNode<Point, PathVertexPtr> {
public:
	/// Dummy constructor
	inline VCMTreeEntry() {}
	inline VCMTreeEntry(PathVertexPtr pathVertexPtr) {
	   position = pathVertexPtr->getPosition();
	}

	/// Return a string representation (for debugging)
	std::string toString() const {}
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

		ref<Sensor> sensor = scene->getSensor();
		const Film *film = sensor->getFilm();
		const Vector2i res = film->getSize();
		scene->getSampler()->advance();
		int pathCount = res.x*res.y;

		const Float etaVCM = (M_PI * radiusSqr) * pathCount;
		const Float misVMWeightFactor = etaVCM;
		const Float misVCWeightFactor = 1.f / etaVCM;


		//////////////////////////////////////////////////////////////////////////
		// Generate light paths
		//////////////////////////////////////////////////////////////////////////

		std::vector<Path* > paths;
		paths.reserve(pathCount);
		m_lightVertices.reserve(pathCount);
		m_lightVertices.clear();
		for(int i = 0; i<pathCount; ++i) {
			Float time = i*1000;
			//Log(EInfo, "%i", i);
			Path* emitterPath = new Path();
			//Log(EInfo, "Attempting initialization");
			emitterPath->initialize(scene, time, EImportance, m_pool);
			//Log(EInfo, "path initialized");
			//m_config.dump()
			emitterPath->randomWalk(scene, scene->getSampler(), m_config.maxDepth, m_config.rrDepth, EImportance, m_pool );

			Float dVCM = 0;
			Float dVC = 0;
			Float dVM = 0;
			Spectrum throughput;

			// skip Emitter Supernode
			for(int vertexIdx = 1; vertexIdx < emitterPath->vertexCount(); vertexIdx++) {

				PathVertexPtr vertex = emitterPath->vertex(vertexIdx);
				Log(EInfo, "Vertex type %d", vertex->type);
				if (!(vertex->type & PathVertex::ENormal)) {
				    continue;
				}


				const Intersection &its = vertex->getIntersection();
				DirectSamplingRecord dRec(its);
				BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));
				// if it's on sensor, sample MIS values
				if (vertexIdx == 1) {
				    // not sure if it's exactly this
				    // if not, look around
				    // src/libbidir/vertex.cpp line 806-824
				    throughput = vertex->weight[ERadiance];
				    DirectionSamplingRecord dirRec(its);
				    Float emissionPdf = dirRec.pdf;

				    dVCM = 1 / emissionPdf;

				    // TODO: handle delta and infinite lights
				    dVC = vertex->getGeometricNormal().z / emissionPdf;
				    dVM = dVC * misVCWeightFactor;
				}
				else {
				    // TODO: handle infinite light
				    dVCM *= its.t * its.t;

				    dVCM /= std::abs(bRec.wi.z);
				    dVC  /= std::abs(bRec.wi.z);
				    dVM  /= std::abs(bRec.wi.z);

				    // TODO: don't store those values if BSDF
				    // is purely specular
				    vertex->data[EdVCMData] = dVCM;
				    vertex->data[EdVCData]  = dVC;
				    vertex->data[EdVMData]  = dVM;
				    m_lightVertices.push_back(vertex);
				}
			}


			scene->getSampler()->advance();
			//Log(EInfo, "walk done");
			paths.push_back(emitterPath);
		}


		/////////////////////////////////////////////////////////////////////////
		// BUILD SEARCH STRUCT
		/////////////////////////////////////////////////////////////////////////

		m_tree.reserve(pathCount);
		for(int i=0; i<m_lightVertices.size(); i++) {
		    m_tree.push_back(m_lightVertices[i]);
		}

		for(int i = 0; i<pathCount; ++i){
			paths[i]->release(m_pool);
			delete paths[i];
		}



		return true;
	}

	MTS_DECLARE_CLASS()
private:
	ref<ParallelProcess> m_process;
	std::vector<PathVertexPtr> m_lightVertices;
	PointKDTree<VCMTreeEntry> m_tree;
	VCMConfiguration m_config;
	MemoryPool m_pool;
};

MTS_IMPLEMENT_CLASS_S(VCMIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(VCMIntegrator, "Vertex Connection And Merging Integrator");
MTS_NAMESPACE_END
