
#if !defined(__VCM_H)
#define __VCM_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/renderproc.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/core/bitmap.h>


MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                         Configuration storage                        */
/* ==================================================================== */

/**
 * \brief Stores all configuration parameters of the
 * vertex connection and merging tracer
 */
struct VCMConfiguration {
	int maxDepth;
	int rrDepth;
	int iterationCount;
	bool useVC;
	bool useVM;
	bool lightTraceOnly;


	inline VCMConfiguration() { }

	inline VCMConfiguration(Stream *stream){
		maxDepth = stream->readInt();
		rrDepth = stream->readInt();
		iterationCount = stream->readInt();
		useVC = stream->readBool();
		useVM = stream->readBool();
		lightTraceOnly = stream->readBool();
	}

	inline void serialize(Stream *stream) const {
		stream->writeInt(maxDepth);
		stream->writeInt(rrDepth);
		stream->writeInt(iterationCount);
		stream->writeBool(useVC);
		stream->writeBool(useVM);
		stream->writeBool(lightTraceOnly);
	}

	void dump() const {
		SLog(EInfo, "Vertex connection and merging configuration:");
		SLog(EInfo, "   Maximum path depth          : %i", maxDepth);
		SLog(EInfo, "   Russian roulette depth      : %i", rrDepth);
		SLog(EInfo, "   Iteration count             : %i", iterationCount);
	}
};

MTS_NAMESPACE_END

#endif /* __VCM_H */
