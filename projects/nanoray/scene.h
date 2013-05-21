#pragma once

#include "core/core.h"
#include "core/maths.h"
#include "core/skylight.h"

#include "primitives.h"

#include <vector>

struct Light
{
	enum Type
	{
		ePoint = 0,
		eSpot = 1,
		eArea = 2
	};

	Light() 
		: m_type(eArea)
		, m_numSamples(0)
		, m_primitive(NULL)
	{
	}
	
	Light(Type t, uint32_t numSamples=1, const Primitive* primitive=NULL) 
		: m_type(t)
		, m_numSamples(numSamples)
		, m_primitive(primitive) 
	{}

	Type m_type;

	// shape to use for emission
	const Primitive* m_primitive;
	// number of samples to take
	uint32_t m_numSamples;
};

class Volume
{
public:

	virtual ~Volume() {};

	virtual void Eval(const Point3& p, const Vector3& dir, float maxT, Colour& transmittance, Colour& emission, Colour& inScatter) const = 0;

};

inline float PowerHeuristic(uint32_t nf, float fPdf, uint32_t ng, float gPdf)
{
	float f = nf*fPdf;
	float g = ng*gPdf;

	float ff = f*f;
	float gg = g*g;

	float r = 1.0f / (ff + gg);

	if (r > 0.0f)
		return (ff) *r;
	else
		return 1.0f;
}

inline float BalanceHeuristic(uint32_t nf, float fPdf, uint32_t ng, float gPdf)
{
	float f = nf*fPdf;
	float g = ng*gPdf;

	return (f) / (f + g);
}

#if 0

#define validate3(c)\
{\
	if (!_finite(c.r) ||\
		!_finite(c.g) ||\
		!_finite(c.b) ||\
		!_finite(c.a) ||\
		_isnan(c.r) ||\
		_isnan(c.g) ||\
		_isnan(c.b) ||\
		_isnan(c.a))\
	{\
	std::cout << "Colour failed validation " << c.r << ", " << c.g << ", " << c.b << #c << "line: " << __LINE__ << std::endl;\
	}\
}\


#define validate(f)\
{\
	if (!_finite(f) ||\
		_isnan(f))\
	{\
		std::cout << "float failed to validate" << #f << " = " << f << __LINE__ << std::endl;\
	}\
}

#else

#define validate3(c)
#define validate(f)

#endif


class Scene
{
public:

	// contiguous buffer for the data + 1 sentinel uint8_t at the end
	typedef std::vector<const Primitive*> PrimitiveArray;
	typedef std::vector<const Volume*> VolumeArray;
	typedef std::vector<const Light*> LightArray;

	PrimitiveArray m_primitives;	
	VolumeArray m_volumes;
	LightArray m_lights;
	
	float m_skyTheta;
	float m_skyPhi;
	float m_skyTurbidity;

	Scene() 
		: m_skyTheta(kPi/2.1f)
		, m_skyPhi(kPi/1.5f)
		, m_skyTurbidity(2.0f)
	{		
	}
	
	void AddPrimitive(const Primitive* p)
	{
		assert(p);
		m_primitives.push_back(p);
	}

	void AddVolume(const Volume* v)
	{
		assert(v);
		m_volumes.push_back(v);
	}

	void AddLight(const Light* l)
	{
		assert(l);
		m_lights.push_back(l);
	}	

	void SetSkyParams(float theta, float phi, float turbidity)
	{
		m_skyTheta = theta;
		m_skyPhi = phi;
		m_skyTurbidity = turbidity;
	}
			
	// trace a ray against the scene returning the closest intersection
	bool Trace(const Point3& rayOrigin, const Vector3& rayDir, float& outT, Vector3& outNormal, const Primitive** outPrimitive) const
	{
		// disgard hits closer than this distance to avoid self intersection artifacts
		const float kEpsilon = 0.001f;

		float minT = FLT_MAX;
		const Primitive* closestPrimitive = NULL;
		Vector3 closestNormal(0.0f);

		for (PrimitiveArray::const_iterator iter=m_primitives.begin(), end=m_primitives.end(); iter != end; ++iter)
		{
			float t;
			Vector3 n;

			const Primitive* p = *iter;

			if (p->Intersect(rayOrigin, rayDir, t, &n))
			{
				if (t < minT && t > kEpsilon)
				{
					minT = t;
					closestPrimitive = p;
					closestNormal = n;
				}
			}
		}
		
		outT = minT;		
		outNormal = closestNormal;
		*outPrimitive = closestPrimitive;

		return closestPrimitive != NULL;
	}

	// illum output represents contribution from inscattering and self-illumination
	void SampleVolumes(const Point3& p, const Vector3& dir, float maxT, Colour& attenuation, Colour& illum) const
	{
		attenuation = Colour(1.0f, 1.0f, 1.0f, 1.0f);
		illum = Colour(0.0f, 0.0f, 0.0f, 0.0f);

		for (uint32_t i=0; i < m_volumes.size(); ++i)
		{
			Colour e;
			Colour s;
			Colour a;

			m_volumes[i]->Eval(p, dir, maxT, a, e, s);

			attenuation += a;
			illum += e + s;
		}
	}

	// returns the radiance arriving in the direction wi at the sample point
	Colour SampleLights(const Point3& p, const Vector3& n, const Vector3& wo, const BRDF* brdf) const 
	{	
		Colour sum(0.0f);

		for (uint32_t i=0; i < m_lights.size(); ++i)
		{
			const Light& l = *m_lights[i];
					
			assert(l.m_primitive);

			// assume all lights are area lights for now
			const Primitive* shape = l.m_primitive;
			
			Colour L(0.0f);

			for (uint32_t s=0; s < 1; ++s)
			{
				// sample light source
				Point3 lp = shape->Sample(NULL);
				Vector3 wi = Normalize(lp-p);

				// check visibility
				float t;
				Vector3 ln;
				const Primitive* h;
				if (Trace(p, wi, t, ln, &h))
				{
					// did we hit the light prim?
					if (h == shape)
					{
						const Colour f = brdf->F(wo, wi);
						validate3(f);

						// light pdf
						const float nl = Clamp(Dot3(ln, -wi), 0.0f, 1.0f);
						validate(nl);
						
						if (nl > 0.0 && f.a  != Colour::kBlack)
						{
							const float lightPdf = (t*t) / ( nl * h->Area() );
							validate(lightPdf);

							// calculate brdf pdf
							const float brdfPdf = brdf->Pdf(brdf->WorldToLocal(wo), brdf->WorldToLocal(wi));
							validate(brdfPdf);							

							// MIS weight
							const float weight = PowerHeuristic(1, lightPdf, 1, brdfPdf);							
							validate(weight);

							L += f * shape->m_emission * Clamp(Dot3(wi, n), 0.0f, 1.0f) * weight / lightPdf;
							validate3(L);
						}
					}
				}		
				
				// sample brdf
				float brdfPdf;
				brdf->Sample(wo, wi, brdfPdf);

				validate(brdfPdf);

				// test visibility from point to light using brdf sample
				if (Trace(p, wi, t, ln, &h))
				{
					if (h == shape)
					{
						const Colour f = brdf->F(wo, wi);

						validate3(f);

						// light pdf
						const float nl = Clamp(Dot3(ln, -wi), 0.0f, 1.0f);
						validate(nl);

						if (nl > 0.0f)
						{
							const float lightPdf = (t*t) / ( nl * h->Area() );
							validate(lightPdf);						

							// MIS weight
							const float weight = PowerHeuristic(1, brdfPdf, 1, lightPdf);
							validate(weight);
	
							L += f * shape->m_emission * Clamp(Dot3(wi, n), 0.0f, 1.0f) * weight / brdfPdf;			
							validate3(L);
						}
					}
				}					
			}
		
			//sum += L / l.m_numSamples;
			sum += L;
		}

		return sum;
	}

	Colour SampleSky(const Vector3& rayDir) const
	{
		if (m_skyTurbidity < 0.0f)
			return Colour(0.0f);
		else
			return SkyLight((1.0f-Clamp(rayDir.y, 0.0f, 1.0f))*kPi*0.5f, atan2(rayDir.z, rayDir.x), m_skyTheta, m_skyPhi, m_skyTurbidity);
	}
};
