#pragma once

#include "core/core.h"
#include "core/maths.h"
#include "core/perlin.h"

class Texture
{
public:

	virtual ~Texture() {}
	virtual Colour Evaluate(const Point3& p) const = 0;

};

class ConstantTexture : public Texture
{
public:

	ConstantTexture(const Colour& c) : m_colour(c) {}

	virtual Colour Evaluate(const Point3& p) const
	{
		return m_colour;
	}

	Colour m_colour;
};

class MixTexture : public Texture
{
public:

	MixTexture(const Texture* tex1, const Texture* tex2, const Texture* mixTex) : m_tex1(tex1), m_tex2(tex2), m_mixTex(mixTex) {}

	virtual Colour Evaluate(const Point3& p) const
	{
		return Lerp(m_tex1->Evaluate(p), m_tex2->Evaluate(p), Clamp(Abs(m_mixTex->Evaluate(p).r), 0.0f, 1.0f));
	}
	
	const Texture* m_tex1;
	const Texture* m_tex2;
	const Texture* m_mixTex;
};

class PerlinTexture : public Texture
{
public:

	PerlinTexture(float frequency, float persistance=1.0f, uint32_t octaves=1) : m_frequency(frequency),
																			   m_persistance(persistance),
																			   m_octaves(octaves)
	{
	}

	virtual Colour Evaluate(const Point3& p) const
	{
		float n = 0.5f + 0.5f*(Perlin3D(p.x*m_frequency, p.y*m_frequency, p.z*m_frequency, m_octaves, m_persistance));

		return Colour(n, n, n);
	}

	float m_frequency;
	float m_persistance;
	uint32_t m_octaves;
	
};

class CheckerboardTexture : public Texture
{
public:

	CheckerboardTexture(float size, const Colour& c1, const Colour& c2) : m_invSize(1.0f/size), m_colour1(c1), m_colour2(c2)
	{
	}

	Colour Evaluate(const Point3& p) const
	{
		int x = (int)floorf(p.x*m_invSize);
		int y = (int)floorf(p.y*m_invSize);
		int z = (int)floorf(p.z*m_invSize);
		
		if ((x+y+z)&0x1)
		{
			return m_colour1;
		}
		else
		{
			return m_colour2;
		}		
	}

	Colour m_colour1;
	Colour m_colour2;
	float m_invSize;	
};
