
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
			//Log(EInfo, "%i", i);
			//Sleep(50);
			Path* emitterPath = new Path();
			//Log(EInfo, "Attempting initialization");
			emitterPath->initialize(scene, time, EImportance, m_pool);
			//Log(EInfo, "path initialized");
			//m_config.dump()
			emitterPath->randomWalk(scene, scene->getSampler(), m_config.maxDepth, m_config.rrDepth, EImportance, m_pool );
			scene->getSampler()->advance();
			//Log(EInfo, "walk done");
			paths.push_back(emitterPath);
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
	VCMConfiguration m_config;
	MemoryPool m_pool;
};

MTS_IMPLEMENT_CLASS_S(VCMIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(VCMIntegrator, "Vertex Connection And Merging Integrator");
MTS_NAMESPACE_END
