#pragma once

#include "core/core.h"
#include "core/maths.h"
#include "core/memory.h"

#include "textures.h"

extern _declspec(thread) MemoryArena* g_memArena;

class BRDF
{
public:

	BRDF(const Point3& p, const Vector3& n)
	{
		m_localToWorld = TransformFromVector(n, p);
		m_worldToLocal = AffineInverse(m_localToWorld);
	}
	
	virtual ~BRDF(){};
	
	// all vectors in world space
	virtual Colour F(const Vector3& wo, const Vector3& wi) const=0;
	// all vectors in world space
	virtual void Sample(const Vector3& wo, Vector3& wi, float& pdf) const=0;
	
	virtual float Pdf(const Vector3& wo, const Vector3& wi) const=0;

	Vector3 WorldToLocal(const Vector3& v) const { return m_worldToLocal * v; }
	Vector3 LocalToWorld(const Vector3& v) const { return m_localToWorld * v; }

	Matrix44 m_localToWorld;	
	Matrix44 m_worldToLocal;
};

class Lambert : public BRDF
{
public:

	Lambert(const Point3& p, const Vector3& n, const Colour& r) : BRDF(p, n), m_rOverPi(r*kInvPi) {}

	virtual Colour F(const Vector3& wi, const Vector3& wo) const	
	{
		return m_rOverPi;
	}

	virtual void Sample(const Vector3& woWorld, Vector3& wiWorld, float& pdf) const
	{
		// generate a sample on the hemisphere weighted by cosTheta term
		Vector3 wiLocal = CosineSampleHemisphere();

		pdf = Pdf(Vector3(0.0), wiLocal);

		wiWorld = m_localToWorld * wiLocal;
	}

	virtual float Pdf(const Vector3& wo, const Vector3& wi) const
	{
		float pdf = wi.z * kInvPi;

		if (pdf == 0.0f)
			pdf = FLT_MAX;	

		return pdf;
	}

	Colour m_rOverPi;

};

// cosTheta should be the angle between the wi and wh
inline Colour Schlick(const Colour& c, float cosTheta)
{
	return c + (Colour(1.0f, 1.0f, 1.0f)-c)*powf(1.0f-cosTheta, 5.0f);
}

class Blinn : public BRDF
{
public:

	Blinn(const Point3& p, const Vector3& n, const Colour& r, float e) : BRDF(p, n), m_reflectance(r), m_exponent(e)
	{

	}
	
	virtual Colour F(const Vector3& wi, const Vector3& wo) const
	{
		Vector3 n = m_localToWorld.GetCol(2);

		// calculate half-angle
		Vector3 wh = Normalize(wi+wo);

		float NdotWh = Abs(Dot(wh, n));
		float NdotWo = Abs(Dot(wo, n));
		float NdotWi = Abs(Dot(wi, n));
		float WodotWh = Abs(Dot(wo, wh));

		//if (Dot(wo, wi) < 0.0f)
			//return Colour::kBlack;

		Colour f = Schlick(m_reflectance, WodotWh);

		// geometric term
		float g = Min(1.0f, Min((2.0f * NdotWh * NdotWo / WodotWh),
						     	(2.0f * NdotWh *NdotWi / WodotWh)));
	
		float d = (m_exponent + 2.0f) * kInv2Pi * powf(Abs(Dot(wh, n)), m_exponent);

		return f  * d * g / (4.0f * NdotWi * NdotWo + 1.0e-4f);
	}

	Vector3 SphericalDirection(float sinTheta, float cosTheta, float phi) const
	{
		return Vector3(sinTheta * cosf(phi),
	              	   sinTheta * sinf(phi),
				  	   cosTheta);
	}
	
	/*
	// uniform sampling reference
	virtual void Sample(const Vector3& woWorld, Vector3& wiWorld, float& pdf) const
	{
		Vector3 wiLocal = UniformSampleSphere();
		if (wiLocal.z < 0.0)
			wiLocal.z *= -1.0f;

		wiWorld = m_localToWorld * wiLocal;
		pdf = kInv2Pi;
	}

	virtual float Pdf(const Vector3& wo, const Vector3& wi) const
	{
		return kInv2Pi;
	}
	*/

	// importance sampling the BRDF
	virtual void Sample(const Vector3& woWorld, Vector3& wiWorld, float& pdf) const
	{
		Vector3 woLocal = m_worldToLocal*woWorld;
		Vector3 wiLocal;

		float u1 = Randf();
		float u2 = Randf();

		float costheta = powf(u1, 1.f / (m_exponent+1.f));
		float sintheta = sqrtf(Max(0.f, 1.f - costheta*costheta));
		float phi = u2 * 2.f * kPi;
		Vector3 H = SphericalDirection(sintheta, costheta, phi);
		
		if (Dot(woLocal, H) < 0.f)
			H.z *= -1.0f;

		// Compute incident direction by reflecting about $\wh$
		wiLocal = -woLocal + 2.f * Dot(woLocal, H) * H;

		float blinn_pdf = ((m_exponent + 1.f) *  powf(costheta, m_exponent)) / (2.f * kPi * 4.f * Max(1.0e-6f, Dot(woLocal, H)));
		pdf = blinn_pdf;

		wiWorld = m_localToWorld * wiLocal;
	}

	virtual float Pdf(const Vector3& wo, const Vector3& wi) const
	{
		Vector3 H = Normalize(wo + wi);
		float cosTheta = Abs(H.z);

		// Compute PDF for wi from Blinn distribution
		float blinn_pdf = ((m_exponent + 1.f) * powf(cosTheta, m_exponent)) / (2.f * kPi * 4.f *  Abs(Dot(wo, H)));
		return blinn_pdf;
	}
	
	Colour m_reflectance;
	float m_exponent;
};

class Material
{
public:

	virtual ~Material() = 0 {} ;
	virtual BRDF* GetBRDF(const Point3& p, const Vector3& n) const = 0;	
};

// pure Lambertian surface
class MatteMaterial : public Material
{
public:
	
	MatteMaterial(const Texture* r) : m_reflectance(r)
	{
	}

	virtual BRDF* GetBRDF(const Point3& p, const Vector3& n) const
	{
		return new (*g_memArena) Lambert(p, n, m_reflectance->Evaluate(p));
	}

	const Texture* m_reflectance;

};

class PlasticMaterial : public Material
{
public:

	PlasticMaterial(const Texture* r, float gloss=30.0f) : m_reflectance(r), m_gloss(gloss)
	{
	}

	virtual BRDF* GetBRDF(const Point3& p, const Vector3& n) const
	{
		return new (*g_memArena) Blinn(p, n, m_reflectance->Evaluate(p), m_gloss);
	}

	float m_gloss;
	const Texture* m_reflectance;
};
