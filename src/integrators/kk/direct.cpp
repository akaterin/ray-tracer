#include <mitsuba/render/scene.h>
#include "direct.h"
#include <vector>
#include <iostream>

MTS_NAMESPACE_BEGIN

struct SamplingValue{
	Spectrum value;
	Spectrum bsdfValue;
	Float pdf;
	Float otherPdf;
};

class DirectIntegrator : public SamplingIntegrator {
public:
	DirectIntegrator(const Properties &props) : SamplingIntegrator(props) {
		brdfSamplesNo = props.getSize("bsdfSamples", 1);
		lightSamplesNo = props.getSize("lightSamples", 1);
		hideEmitters = props.getBoolean("hideEmitters", false);
	}

	DirectIntegrator(Stream *stream, InstanceManager *manager)
	 : SamplingIntegrator(stream, manager) {
		brdfSamplesNo = stream->readSize();
		lightSamplesNo = stream-> readSize();
		hideEmitters = stream->readBool();
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		SamplingIntegrator::serialize(stream, manager);
		stream->writeSize(brdfSamplesNo);
		stream->writeSize(lightSamplesNo);
		stream->writeBool(hideEmitters);
	}

	void configure() {
		SamplingIntegrator::configure();
	}

	void configureSampler(const Scene *scene, Sampler *sampler){
		SamplingIntegrator::configureSampler(scene, sampler);
		if( brdfSamplesNo > 1 )
			sampler->request2DArray(brdfSamplesNo);
		if( lightSamplesNo > 1 )
			sampler->request2DArray(lightSamplesNo);
	}

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		RayDifferential ray(r);
		const Scene *scene = rRec.scene;
		Spectrum Li(0.0f);

		if (!rRec.rayIntersect(ray)) {
			if (rRec.type & RadianceQueryRecord::EEmittedRadiance && !hideEmitters)
				return scene->evalEnvironment(ray);
			else
				return Spectrum(0.0f);
		}


		if (rRec.its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
			Li += rRec.its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

		if (rRec.its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance) && !hideEmitters)
			Li += rRec.its.Le(-ray.d);

		std::vector<SamplingValue> brdfSamples = integrateBrdf(r,rRec);
		std::vector<SamplingValue> lightSamples = integrateLight(r, rRec);
		return computeMIS(brdfSamples, lightSamples, Li);
	}

	std::vector<SamplingValue> integrateBrdf(const RayDifferential &r, RadianceQueryRecord &rRec) const{
		Intersection &its = rRec.its, bsdfIts;
		const BSDF *bsdf = its.getBSDF(r);
		DirectSamplingRecord directRec(its);
		std::vector<SamplingValue> samplingValues;
			
		Point2* samplesArray;
		if( brdfSamplesNo > 1 ){
			samplesArray = rRec.sampler->next2DArray(brdfSamplesNo);
		}else{
			samplesArray = &(rRec.sampler->next2D());
		}

		for (size_t i=0; i<brdfSamplesNo; ++i) {
			Float bsdfPdf;

			BSDFSamplingRecord bsdfRec(its, rRec.sampler, ERadiance);
			Spectrum bsdfVal = bsdf->sample(bsdfRec, bsdfPdf, samplesArray[i]);

			if (bsdfVal.isZero())
				continue;

			const Vector wo = its.toWorld(bsdfRec.wo);
			Float woDotGeoN = dot(its.geoFrame.n, wo);
			if (woDotGeoN * Frame::cosTheta(bsdfRec.wo) <= 0)
				continue;

			Ray bsdfRay(its.p, wo, r.time);

			Spectrum value;
			if (rRec.scene->rayIntersect(bsdfRay, bsdfIts)) {
				if (!bsdfIts.isEmitter())
					continue;

				value = bsdfIts.Le(-bsdfRay.d);
				directRec.setQuery(bsdfRay, bsdfIts);
			} else {
				const Emitter *env = rRec.scene->getEnvironmentEmitter();

				if (!env || (hideEmitters && bsdfRec.sampledType == BSDF::ENull))
					continue;

				value = env->evalEnvironment(RayDifferential(bsdfRay));
				if (!env->fillDirectSamplingRecord(directRec, bsdfRay))
					continue;
			}
		
			Float lightPdf = (!(bsdfRec.sampledType & BSDF::EDelta)) ? rRec.scene->pdfEmitterDirect(directRec) : 0;
			samplingValues.push_back(createSamplingRecord(value, bsdfVal, bsdfPdf, lightPdf));
		}
		return samplingValues;
	}

	std::vector<SamplingValue> integrateLight(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		std::vector<SamplingValue> samplingValues;
		Intersection its = rRec.its;
		const BSDF *bsdf = its.getBSDF(r);

		Point2* samples;
		if( lightSamplesNo > 1 ){
			samples = rRec.sampler->next2DArray(lightSamplesNo);
		}else{
			samples = &(rRec.sampler->next2D());
		}

		DirectSamplingRecord dRec(its);
		if (!bsdf->getType() || !BSDF::ESmooth) {
			return samplingValues;
		}

		for (size_t i=0; i<lightSamplesNo; ++i) {
				Spectrum value = rRec.scene->sampleEmitterDirect(dRec, samples[i]);
				if (!value.isZero()) {
					const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

					BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));

					const Spectrum bsdfVal = bsdf->eval(bRec);
					Float bsdfPdf = emitter->isOnSurface() ? bsdf->pdf(bRec) : 0;

					if (!bsdfVal.isZero() && (dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {
						samplingValues.push_back(createSamplingRecord(value, bsdfVal, dRec.pdf, bsdfPdf));
					}
				}
		}
		return samplingValues;
	}

	Spectrum computeMIS(std::vector<SamplingValue> brdfSamples, std::vector<SamplingValue> lightSamples, Spectrum Li) const{
		Float cBrdf = brdfSamplesNo/(Float)(brdfSamplesNo + lightSamplesNo);
		Float cLight = lightSamplesNo/(Float)(brdfSamplesNo + lightSamplesNo);
		for( int i = 0; i < brdfSamples.size(); ++i ){
			Float weight = cBrdf*brdfSamples[i].pdf/( cBrdf*brdfSamples[i].pdf + cLight*brdfSamples[i].otherPdf );
			Li += weight*brdfSamples[i].value * brdfSamples[i].bsdfValue;
		}

		for( int i = 0; i < lightSamples.size(); ++i ){
			Float weight = cLight*lightSamples[i].pdf/( cBrdf*lightSamples[i].otherPdf + cLight*lightSamples[i].pdf );
			Li += weight*lightSamples[i].value * lightSamples[i].bsdfValue;
		}
		return Li;
	}
	
	MTS_DECLARE_CLASS()
private:

	SamplingValue createSamplingRecord(Spectrum value, Spectrum bsdfValue, Float pdf, Float otherPdf) const{
		SamplingValue record;
		record.value = value;
		record.bsdfValue = bsdfValue;
		record.pdf = pdf;
		record.otherPdf = otherPdf;
		return record;
	}

	size_t brdfSamplesNo;
	size_t lightSamplesNo;
	bool hideEmitters;

};

MTS_IMPLEMENT_CLASS_S(DirectIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(DirectIntegrator, "KK - Direct illumination integrator");

MTS_NAMESPACE_END