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

struct VCMVertex {
    VCMVertex() {}
    VCMVertex(const Point &p, const Spectrum& throughput,
	    uint32_t pathLength, float dvcm, float dvc, float dvm,
	    const BSDF* bsdf) :
	mPosition(p), mThroughput(throughput),
	mPathLength(pathLength), dVCM(dvcm),
	dVC(dvc), dVM(dvm), mBsdf(bsdf)
    {}
    Point mPosition;
    Spectrum mThroughput; // Path throughput (including emission)
    uint32_t  mPathLength; // Number of segments between source and vertex

    // Stores all required local information, including incoming direction.
    // TODO: will it boom as it's a pointer?
    const BSDF* mBsdf;
    Float dVCM; // MIS quantity used for vertex connection and merging
    Float dVC;  // MIS quantity used for vertex connection
    Float dVM;  // MIS quantity used for vertex merging
};

// to make easier to carry those values around
struct SubPathState
{
    Spectrum mThroughput;         // Path throughput
    bool  mIsFiniteLight;
    bool  mSpecularPath;

    Float dVCM; // MIS quantity used for vertex connection and merging
    Float dVC;  // MIS quantity used for vertex connection
    Float dVM;  // MIS quantity used for vertex merging
};
struct VCMTreeEntry :
    public SimpleKDNode<Point, VCMVertex> {
	public:
	    /// Dummy constructor
	    inline VCMTreeEntry() {}
	    inline VCMTreeEntry(const VCMVertex& v) {
		position = v.mPosition;
		data = v;
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
		Path* emitterPath = new Path();
		emitterPath->initialize(scene, time, EImportance, m_pool);
		emitterPath->randomWalk(scene, scene->getSampler(), m_config.maxDepth, m_config.rrDepth, EImportance, m_pool );
		SubPathState pathState;

		// skip Emitter Supernode
		for(int vertexIdx = 1; vertexIdx < emitterPath->vertexCount(); vertexIdx++) {
		    PathVertexPtr vertex = emitterPath->vertex(vertexIdx);
		    if (!(vertex->type & PathVertex::ENormal)) {
			continue;
		    }

		    // if it's on emitter, sample MIS values
		    if (vertexIdx == 1 && emitterPath->vertexCount() > 2) {
			PathVertexPtr nextVertex = emitterPath->vertex(2);
			if (! (nextVertex -> getType() & PathVertex::ESurfaceInteraction) ) {
			    break;
			}
			PositionSamplingRecord pRec = vertex->getPositionSamplingRecord();

			Assert( vertex->type & PathVertex::EEmitterSample );
			//PathEdge* path = emitterPath->edge(1);

			//Float emissionPdf = path->pdf[ERadiance];
			//Float directPdf = pRec.pdf; //scene->pdfEmitterDirect(dRec);
			const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
			Float lightPickProb = scene->pdfEmitterDiscrete(emitter);

			pathState.mIsFiniteLight = !(emitter->getType() & AbstractEmitter::EDeltaDirection);

			Vector wo = emitterPath->vertex(2)->getPosition() - pRec.p;
			Float dist = wo.length(); wo /= dist;
			Float cosTheta = std::abs(dot(pRec.n, wo));

			DirectionSamplingRecord dRec(wo, ESolidAngle);

			Float directPdf = emitter->pdfPosition(pRec) * lightPickProb;
			Float emissionPdf = emitter->pdfDirection(dRec, pRec) * lightPickProb;
			//Log(EInfo, "Direct= %f, emission= %f", directPdf, emissionPdf);

			pathState.mThroughput = vertex->weight[ERadiance] / emissionPdf;

			pathState.dVCM = directPdf / emissionPdf;

			// TODO: handle delta and infinite lights
			pathState.dVC = cosTheta / emissionPdf;
			pathState.dVM = pathState.dVC * misVCWeightFactor;
		    }
		    else {
			// there's much more to computing MIS values here
			// (SampleScattering @ SmallVCM)
			if (! (vertex->getType() & PathVertex::ESurfaceInteraction) ) {
			//    Log(EInfo, "Other vertex type: %d", vertex->getType());
			    continue;
			}
			const Intersection &its = vertex->getIntersection();
			DirectSamplingRecord dRec(its);

			// is this the right cos angle to compute?
			// should it be dRec.n instead?
			Float cosTheta = std::abs(dot(dRec.d, dRec.refN));

			if (vertexIdx > 2 || pathState.mIsFiniteLight) {
			    pathState.dVCM *= its.t * its.t;
			}

			pathState.dVCM /= cosTheta;
			pathState.dVC  /= cosTheta;
			pathState.dVM  /= cosTheta;

			const BSDF* bsdf = its.getBSDF();

			// add vertex iff the bsdf is not purely specular
			// not sure if thats how you check it ;)
			if ( !(bsdf->getType() & BSDF::EDiffuseReflection
				    || bsdf->getType() & BSDF::EDiffuseTransmission) ) {
			    VCMVertex v( vertex->getPosition(),
				 pathState.mThroughput,
				 vertexIdx - 1,
				 pathState.dVCM, pathState.dVC, pathState.dVM, its.getBSDF() );
			    m_lightVertices.push_back(v);
			}
		    }
		}

		connectToEye(scene, time, emitterPath);

		scene->getSampler()->advance();
		//Pytanie czy ścieżki, których się nie udało połączyć z okiem też tu wrzucać?
		paths.push_back(emitterPath);
	    }

	    Log(EInfo, "Starting to build...");

	    /////////////////////////////////////////////////////////////////////////
	    // BUILD SEARCH STRUCT
	    /////////////////////////////////////////////////////////////////////////

	    m_tree.reserve(pathCount);
	    for(int i=0; i<m_lightVertices.size(); i++) {
		m_tree.push_back(m_lightVertices[i]);
	    }


	    for(int i = 0; i<pathCount; ++i) {
		// Generate eye path
		Path *sensorPath = new Path();
		sensorPath->initialize(scene, 0, ERadiance, m_pool);
		Point2i startPosition = Point2i(i % res.x, (i - i % res.x) / res.x);
		//Log(EInfo,"%i %i %i", i, startPosition.x, startPosition.y);
		//Sleep(100);
		sensorPath->randomWalkFromPixel(scene, scene->getSampler(), m_config.maxDepth, startPosition, m_config.rrDepth, m_pool);
		scene->getSampler()->advance();
	    }


	    for(int i = 0; i<pathCount; ++i){
		paths[i]->release(m_pool);
		delete paths[i];
	    }

	    return true;
	}

	bool connectToEye(Scene *scene, Float time, Path *emitterPath) {
	    int prevSize = emitterPath->vertexCount();
	    Path* sensorPath = new Path();
	    sensorPath->initialize(scene, time, ERadiance, m_pool);
	    sensorPath->randomWalk(scene, scene->getSampler(), 1, 0, ERadiance, m_pool);
	    PathVertex* vertexOnEye = sensorPath->vertex(1);
	    PathVertex* succOnSensor = sensorPath->vertex(0);
	    PathEdge* succEdge = sensorPath->edge(0);
	    PathEdge* lastEdge = emitterPath->edge(emitterPath->edgeCount()-1);
	    PathVertex* predLastVertex = emitterPath->vertexOrNull(emitterPath->vertexCount()-2);
	    PathVertex* lastVertex = emitterPath->vertex(emitterPath->vertexCount()-1);
	    PathEdge* newEgde = new PathEdge();
	    bool succeded = PathVertex::connect(scene, predLastVertex, lastEdge, lastVertex, newEgde, vertexOnEye, succEdge, succOnSensor);
	    if( !succeded )
		return false;
	    emitterPath->append(newEgde, vertexOnEye);
	    emitterPath->append(succEdge, succOnSensor);
	    return succeded;
	}

	MTS_DECLARE_CLASS()
    private:
	    ref<ParallelProcess> m_process;
	    std::vector<VCMVertex> m_lightVertices;
	    PointKDTree<VCMTreeEntry> m_tree;
	    VCMConfiguration m_config;
	    MemoryPool m_pool;
};

MTS_IMPLEMENT_CLASS_S(VCMIntegrator, false, Integrator)
    MTS_EXPORT_PLUGIN(VCMIntegrator, "Vertex Connection And Merging Integrator");
    MTS_NAMESPACE_END
