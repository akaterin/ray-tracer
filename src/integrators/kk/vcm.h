
#if !defined(__VCM_H)
#define __VCM_H

#include <mitsuba/mitsuba.h>


MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                         Configuration storage                        */
/* ==================================================================== */

/**
 * \brief Stores all configuration parameters of the
 * vertex connection and merging tracer
 */
struct VCMConfiguration {
	inline VCMConfiguration() { }

	inline VCMConfiguration(Stream *stream){

	}

	inline void serialize(Stream *stream) const {
	}

	void dump() const {
	}
};

MTS_NAMESPACE_END

#endif /* __VCM_H */
