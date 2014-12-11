
#include <mitsuba/bidir/vertex.h>
#include <mitsuba/bidir/edge.h>
#include <mitsuba/bidir/path.h>
#include <vector>
#include <utility>
#include <windows.h>
#include "vcm.h"

MTS_NAMESPACE_BEGIN

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

		ref<Sensor> sensor = scene->getSensor();
		const Film *film = sensor->getFilm();
		const Vector2i res = film->getSize();
		scene->getSampler()->advance();
		int pathCount = res.x*res.y;


		//////////////////////////////////////////////////////////////////////////
		// Generate light paths
		//////////////////////////////////////////////////////////////////////////

		std::vector<Path* > paths;
		paths.reserve(pathCount);
		for(int i = 0; i<pathCount; ++i) {
			Float time = i*1000;
			Path* emitterPath = new Path();
			emitterPath->initialize(scene, time, EImportance, m_pool);
			emitterPath->randomWalk(scene, scene->getSampler(), m_config.maxDepth, m_config.rrDepth, EImportance, m_pool );

			connectToEye(scene, time, emitterPath);

			scene->getSampler()->advance();
			//Pytanie czy ścieżki, których się nie udało połączyć z okiem też tu wrzucać?
			paths.push_back(emitterPath);
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
	VCMConfiguration m_config;
	MemoryPool m_pool;
};

MTS_IMPLEMENT_CLASS_S(VCMIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(VCMIntegrator, "Vertex Connection And Merging Integrator");
MTS_NAMESPACE_END
