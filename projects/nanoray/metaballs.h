#pragma once

#include "core/maths.h"
#include "core/core.h"

class Metaballs
{
public:

	struct Ball
	{
		Point3 m_pos;
		float m_radius;
	};

	Metaballs(const Ball* balls, uint32_t n)
	{
		m_numBalls = n;
		m_balls = new Ball[n];
		memcpy(m_balls, balls, sizeof(Ball)*n);		

		// calculate the sum of the radii used to approximate the distance field
		for (uint32_t i=0; i < n; ++i)
		{
			m_radiiSum += balls[i].m_radius;
		}
	}

	~Metaballs()
	{
		delete[] m_balls;
	}


	// evaluate the distance of a point to the iso-surface of the surface (lower-bound)
	float Distance(const Point3& p, float threshold)
	{
		float density = 0.0f;
		for (uint32_t i=0; i < m_numBalls; ++i)
		{
			density += Density(Length(m_balls[i].m_pos-p), m_balls[i].m_radius); 
		}

		return 2.0f/3.0f * m_radiiSum * (threshold - density);
	}

private:

	float Density(float r, float rmax)
	{
		if (r > rmax)
			return 0.0f;
		else
		{
			float t = (r/rmax);
			return 2.0f*(t*t*t) - 3.0f*(t*t) + 1.0f;
		}
	}

	float m_radiiSum;
	uint32_t m_numBalls;

	Ball* m_balls;
};