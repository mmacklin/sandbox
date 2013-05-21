/*
 *  Grid3D.h
 *  Fluid
 *
 *  Created by Miles Macklin on 10/27/10.
 *  Copyright 2010 None. All rights reserved.
 *
 */

#pragma once

#include "core/maths.h"

#define TIMERS 1
#if TIMERS

struct ScopedTimer
{
	ScopedTimer(const char* name) : m_startTime(GetSeconds()), m_name(name)
	{
	}
	
	~ScopedTimer()
	{
		cout << m_name << " took: " << (GetSeconds()-m_startTime)*1000.0f << "ms" << endl;
	}
	
	double m_startTime;
	const char* m_name;
};

#else

struct ScopedTimer
{
	ScopedTimer(const char*) {}
};

#endif


// interpolate between b and c using derivatives defined as (c-a) and (d-b)
inline float CubicInterpolate(float a, float b, float c, float d, float t)
{
	const float tt = t*t;
	const float ttt = tt*t;
	
	float dk1 = 0.5f*(c-a);
	float dk2 = 0.5f*(d-b);

	float r = (2.0f*ttt -3.0f*tt + 1.0f)*b + (ttt - 2.0f*tt + t)*dk1 + (-2.0f*ttt + 3.0f*tt)*c + (ttt-tt)*dk2;
	
	//assert(fabsf(r-r2) < 0.001f);
	
	r = Clamp(r, b, c);

	return r;
	
}

/*
 float texture3DCubic(sampler3D sampler, vec3 uvw)
 {
 float i = floor(uvw.x);
 float j = floor(uvw.y);
 float k = floor(uvw.z);
 
 // tricubic interpolation
 float tx = uvw.x-i;
 float ty = uvw.y-j;
 float tz = uvw.z-k;
 
 float a0 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j-1, k-1)).x, texture3D(sampler, vec3(i, j-1, k-1)).x, texture3D(sampler, vec3(i+1, j-1, k-1)).x, texture3D(sampler, vec3(i+2, j-1, k-1)).x, tx);				
 float a1 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j+0, k-1)).x, texture3D(sampler, vec3(i, j+0, k-1)).x, texture3D(sampler, vec3(i+1, j+0, k-1)).x, texture3D(sampler, vec3(i+2, j+0, k-1)).x, tx);				
 float a2 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j+1, k-1)).x, texture3D(sampler, vec3(i, j+1, k-1)).x, texture3D(sampler, vec3(i+1, j+1, k-1)).x, texture3D(sampler, vec3(i+2, j+1, k-1)).x, tx);				
 float a3 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j+2, k-1)).x, texture3D(sampler, vec3(i, j+2, k-1)).x, texture3D(sampler, vec3(i+1, j+2, k-1)).x, texture3D(sampler, vec3(i+2, j+2, k-1)).x, tx);				
 
 float a4 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j-1, k+0)).x, texture3D(sampler, vec3(i, j-1, k+0)).x, texture3D(sampler, vec3(i+1, j-1, k+0)).x, texture3D(sampler, vec3(i+2, j-1, k+0)).x, tx);				
 float a5 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j+0, k+0)).x, texture3D(sampler, vec3(i, j+0, k+0)).x, texture3D(sampler, vec3(i+1, j+0, k+0)).x, texture3D(sampler, vec3(i+2, j+0, k+0)).x, tx);				
 float a6 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j+1, k+0)).x, texture3D(sampler, vec3(i, j+1, k+0)).x, texture3D(sampler, vec3(i+1, j+1, k+0)).x, texture3D(sampler, vec3(i+2, j+1, k+0)).x, tx);				
 float a7 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j+2, k+0)).x, texture3D(sampler, vec3(i, j+2, k+0)).x, texture3D(sampler, vec3(i+1, j+2, k+0)).x, texture3D(sampler, vec3(i+2, j+2, k+0)).x, tx);				
 
 float a8 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j-1, k+1)).x, texture3D(sampler, vec3(i, j-1, k+1)).x, texture3D(sampler, vec3(i+1, j-1, k+1)).x, texture3D(sampler, vec3(i+2, j-1, k+1)).x, tx);				
 float a9 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j+0, k+1)).x, texture3D(sampler, vec3(i, j+0, k+1)).x, texture3D(sampler, vec3(i+1, j+0, k+1)).x, texture3D(sampler, vec3(i+2, j+0, k+1)).x, tx);				
 float a10 = CubicInterpolate(texture3D(sampler, vec3(i-1, j+1, k+1)).x, texture3D(sampler, vec3(i, j+1, k+1)).x, texture3D(sampler, vec3(i+1, j+1, k+1)).x, texture3D(sampler, vec3(i+2, j+1, k+1)).x, tx);				
 float a11 = CubicInterpolate(texture3D(sampler, vec3(i-1, j+2, k+1)).x, texture3D(sampler, vec3(i, j+2, k+1)).x, texture3D(sampler, vec3(i+1, j+2, k+1)).x, texture3D(sampler, vec3(i+2, j+2, k+1)).x, tx);				
 
 float a12 = CubicInterpolate(texture3D(sampler, vec3(i-1, j-1, k+2)).x, texture3D(sampler, vec3(i, j-1, k+2)).x, texture3D(sampler, vec3(i+1, j-1, k+2)).x, texture3D(sampler, vec3(i+2, j-1, k+2)).x, tx);				
 float a13 = CubicInterpolate(texture3D(sampler, vec3(i-1, j+0, k+2)).x, texture3D(sampler, vec3(i, j+0, k+2)).x, texture3D(sampler, vec3(i+1, j+0, k+2)).x, texture3D(sampler, vec3(i+2, j+0, k+2)).x, tx);				
 float a14 = CubicInterpolate(texture3D(sampler, vec3(i-1, j+1, k+2)).x, texture3D(sampler, vec3(i, j+1, k+2)).x, texture3D(sampler, vec3(i+1, j+1, k+2)).x, texture3D(sampler, vec3(i+2, j+1, k+2)).x, tx);				
 float a15 = CubicInterpolate(texture3D(sampler, vec3(i-1, j+2, k+2)).x, texture3D(sampler, vec3(i, j+2, k+2)).x, texture3D(sampler, vec3(i+1, j+2, k+2)).x, texture3D(sampler, vec3(i+2, j+2, k+2)).x, tx);				
 
 
 float b0 = CubicInterpolate(a0, a1, a2, a3, ty);
 float b1 = CubicInterpolate(a4, a5, a6, a7, ty);
 float b2 = CubicInterpolate(a8, a9, a10, a11, ty);
 float b3 = CubicInterpolate(a12, a13, a14, a15, ty);
 
 float c0 = CubicInterpolate(b0, b1, b2, b3, tz);
 
 return c0;
 }
*/ 

/*
inline uint32_t Part1By2(uint32_t n)
{
	n = (n ^ (n << 16)) & 0xff0000ff; 
	n = (n ^ (n << 8)) & 0x0300f00f;
	n = (n ^ (n << 4)) & 0x030c30c3;
	n = (n ^ (n << 2)) & 0x09249249;
	
	return n;
}

inline int Index(uint32_t x, uint32_t y, uint32_t z)
{
	return int((Part1By2(z) << 2) + (Part1By2(y) << 1) + Part1By2(x));
}
*/

// linear index
template <int width, int height, int depth>
inline int Index(int x, int y, int z)
{
	return z*(width*height) + y*width + x;
}

template <int width, int height, int depth>
class XGrid3D
{
public:
	
	enum 
	{
		kWidth = width,
		kHeight = height,
		kDepth = depth
	};
	
	XGrid3D() 
	{
		m_data = new float[width*height*depth];
		Reset();
	}
	
	~XGrid3D()
	{
		delete[] m_data;
	}
	
	float CubicInterp(float x, float y, float z) const
	{
		int i = floorf(x);
		int j = floorf(y);
		int k = floorf(z);
		
		// tricubic interpolation
		float tx = x-i;
		float ty = y-j;
		float tz = z-k;

		// collapse to a 4x4 patch
		float patch[4][4];
		
		int row = 0;
		
		for (int pz=k-1; pz < (k+3); ++pz)
		{
			int col = 0;
			
			for (int py=j-1; py < (j+3); ++py)
			{
				//assert(row >= 0 && row < 4);
				//assert(col >= 0 && col < 4);
				
				patch[row][col] = CubicInterpolate(GetClamped(x-1, py, pz), GetClamped(x, py, pz), GetClamped(x+1, py, pz), GetClamped(x+2, py, pz), tx);				
				
				col++;
			}
			row++;
		}
	
		// collapse patch to a 4x1 segment
		float segment[4];
		
		for (int r=0; r < 4; ++r)
		{
			segment[r] = CubicInterpolate(patch[r][0], patch[r][1], patch[r][2], patch[r][3], ty);
		}
	
		
		// collapse to a point
		float g = CubicInterpolate(segment[0], segment[1], segment[2], segment[3], tz);
		return g;
		
	}
	
	float LinearInterp(float x, float y, float z) const
    {
		int i = floorf(x);
		int j = floorf(y);
		int k = floorf(z);
		
		// trilinear interpolation
		float tx = x-i;
		float ty = y-j;
		float tz = z-k;
		
		float a = Lerp(GetClamped(i, j, k), GetClamped(i, j, k+1), tz);
		float b = Lerp(GetClamped(i+1, j, k), GetClamped(i+1, j, k+1), tz);
		float c = Lerp(GetClamped(i, j+1, k), GetClamped(i, j+1, k+1), tz);		
		float d = Lerp(GetClamped(i+1, j+1, k), GetClamped(i+1, j+1, k+1), tz);
		
		float e = Lerp(a, b, tx);
		float f = Lerp(c, d, tx);
		
		float g = Lerp(e, f, ty);
		
		return g;
	}
	
	float GetUnsafe(int x, int y, int z) const
	{
		return m_data[Index<width, height, depth>(x, y, z)];
	}
	
	float Get(int x, int y, int z) const
    {
        if (x < 0 || x >= width || y < 0 || y >= height || z < 0 || z >= depth)
        {
            return 0.0f;
        }
		else
		{
	        return m_data[Index<width, height, depth>(x, y, z)];
		}
    }
	
	float GetClamped(int x, int y, int z) const
	{
		//return Get(x, y, z);
		
		x = std::max(0, std::min(x, width-1));
		y = std::max(0, std::min(y, height-1));
		z = std::max(0, std::min(z, depth-1));
		
        return m_data[Index<width, height, depth>(x, y, z)];
	}
	
    void Set(int x, int y, int z, float v)
    {
		//assert(!isnan(v));
        //assert(x >= 0 && x < m_width);
        //assert(y >= 0 && y < m_height);
		//assert(z >= 0 && z < m_depth);
		
        m_data[Index<width, height, depth>(x, y, z)] = v;
    }
	
    void Add(int x, int y, int z, float v)
    {
		//assert(!isnan(v));
		
     //   assert(x >= 0 && x < m_width);
       // assert(y >= 0 && y < m_height);
		//assert(z >= 0 && z < m_depth);
		if (x < 0 || x > width-1 ||
			y < 0 || y > height-1 ||
			z < 0 || z > depth-1)
		{
			return;
		}
		
        m_data[Index<width, height, depth>(x, y, z)] += v;
    }
	
	// useful for debugging
    float Sum() const
    {
        float s = 0.0f;
		
        for (int i=0; i < width*height*depth; ++i)
        {
            s += fabsf(m_data[i]);
        }
		
        return s;
    }
	
	float Average() const
	{
		return Sum() / (width*height*depth);
	}
	
    void Reset()
    {
        memset(m_data, 0, sizeof(float)*width*height*depth);
    }
	
	
	float* m_data;
};

template <bool cubic, typename Grid3D>
void AdvectQ(Grid3D& dest, const Grid3D& q, const Grid3D& u, const Grid3D& v, const Grid3D& w, float dt)
{
	for (int z=0; z < Grid3D::kDepth; ++z)
	{
		for (int y=0; y < Grid3D::kHeight; ++y)
		{
			for (int x=0; x < Grid3D::kWidth; ++x)
			{
				
				float px = x - u.Get(x, y, z)*dt;
				float py = y - v.Get(x, y, z)*dt;
				float pz = z - w.Get(x, y, z)*dt;
				
				// interpolated q	
				float nq = (cubic)?q.CubicInterp(px, py, pz):q.LinearInterp(px, py, pz);
				
				dest.Set(x, y, z, nq);
			}
		}
	}
}

template <typename Grid3D>
void CalculateDivergence(Grid3D& d, const Grid3D& u, const Grid3D& v, const Grid3D& w, float dt)
{
	ScopedTimer timer(__FUNCTION__);
	
	for (int z=0; z < Grid3D::kDepth; ++z)
	{
		for (int y=0; y < Grid3D::kHeight; ++y)
		{
			for (int x=0; x < Grid3D::kWidth; ++x)
			{
				
				float dx = u.Get(x+1, y, z)-u.Get(x-1, y, z);
				float dy = v.Get(x, y+1, z)-v.Get(x, y-1, z);
				float dz = w.Get(x, y, z+1)-w.Get(x, y, z-1);
				
				float divergence = 0.5f*(dx + dy + dz);
				
				d.Set(x, y, z, divergence);
			}
		}
	}
}

template <typename Grid3D>
void PressureSolve(Grid3D& out, const Grid3D& divergence, float dt)
{
	ScopedTimer timer(__FUNCTION__);

	Grid3D scratch;
	
    Grid3D* read = &out;
    Grid3D* write = read;//&scratch;
    
    const float scale = dt;
	const float rcpScale = 1.0f / scale;
	
	for (int i=0; i < 40; ++i)
    {
		for (int z=0; z < Grid3D::kDepth; ++z)
		{
			for (int y=0; y < Grid3D::kHeight; ++y)
			{
				for (int x=0; x < Grid3D::kWidth; ++x)
				{
					
					// Jacobi solver
					float p = -divergence.Get(x, y, z);
					
					float a = read->Get(x-1, y, z);
					float b = read->Get(x+1, y, z);
					float c = read->Get(x, y+1, z);
					float d = read->Get(x, y-1, z);
					float e = read->Get(x, y, z-1);
					float f = read->Get(x, y, z+1);
					
					float v = (p+(a+b+c+d+e+f))/6.0f;
					
					write->Set(x, y, z, v);
				}
			}
		}
		
		std::swap(read, write);		
	}
}


template <typename Grid3D>
void PressureApply(Grid3D& u, Grid3D& v, Grid3D& w, const Grid3D& pressure, float dt)
{
	ScopedTimer timer(__FUNCTION__);

	for (int z=0; z < Grid3D::kDepth; ++z)
	{
		for (int y=0; y < Grid3D::kHeight; ++y)
		{
			for (int x=0; x < Grid3D::kWidth; ++x)
			{
				
				// calculate pressure gradient
				float dpdx = 0.5f*(pressure.Get(x+1, y, z)-pressure.Get(x-1, y, z));
				float dpdy = 0.5f*(pressure.Get(x, y+1, z)-pressure.Get(x, y-1, z));
				float dpdz = 0.5f*(pressure.Get(x, y, z+1)-pressure.Get(x, y, z-1));
				
				u.Add(x, y, z, -dpdx);
				v.Add(x, y, z, -dpdy);
				w.Add(x, y, z, -dpdz);
			}
		}
	}
}












