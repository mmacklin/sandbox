/*
 *  Fluid2D.cpp
 *  Fluid
 *
 *  Created by Miles Macklin on 10/19/10.
 *  Copyright 2010 None. All rights reserved.
 *
 */

#include <stdio.h>
#include <memory.h>
#include <assert.h>
#include <float.h>
#include <math.h>

#include <iostream>
#include <algorithm>
#include <vector>

#include <core/maths.h>
#include <core/shader.h>

using namespace std;

#include "blackbody.h"

const int kWidth = 128;
const int kHeight = 128;

float g_zoomx = 1.0f;
float g_zoomy = 1.0f;

//#define TIMERS 1
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
	float e = c-b;
	
	// modified bicubic interpolation to avoid overshoot, check that central difference signs 
	// agree with forward differences otherwise set derivative to zero, guarantees no overshoot from Visual Simulation of Smoke (Fedkiw)
	/*
	if (Sign(dk1) != Sign(e))
		dk1 = 0.0f;
	
	if (Sign(dk2) != Sign(e))
		dk2 = 0.0f;
	*/
	float r = (2.0f*ttt -3.0f*tt + 1.0f)*b + (ttt - 2.0f*tt + t)*dk1 + (-2.0f*ttt + 3.0f*tt)*c + (ttt-tt)*dk2;

	return Clamp(r, b, c);
}


class Grid2D
{
public:
	
	Grid2D(int width, int height) : m_width(width)
								  , m_height(height)
	{
		m_data = new float[width*height];
		
		Reset();
	}

	~Grid2D()
	{
		delete[] m_data;
	}

	float LinearInterp(float x, float y) const
	{
		int i = floorf(x);
		int j = floorf(y);
	
		// bilinear interpolation
		float tx = x-i;
		float ty = y-j;
		
		float a = Lerp(Get(i, j), Get(i+1, j), tx);
		float b = Lerp(Get(i, j+1), Get(i+1, j+1), tx);

		float r = Lerp(a, b, ty);
		return r;
	}

	float CubicInterp(float x, float y) const
	{
		int i = floorf(x);
		int j = floorf(y);

		// bicubic interpolation
		float tx = x-i;
		float ty = y-j;
		
		// interpolate 4 rows horizontally first
		float row = j-1;
		
		float a = CubicInterpolate(GetClamped(i-1, row), GetClamped(i, row), GetClamped(i+1, row), GetClamped(i+2, row), tx); ++row;
		float b = CubicInterpolate(GetClamped(i-1, row), GetClamped(i, row), GetClamped(i+1, row), GetClamped(i+2, row), tx); ++row;
		float c = CubicInterpolate(GetClamped(i-1, row), GetClamped(i, row), GetClamped(i+1, row), GetClamped(i+2, row), tx); ++row;
		float d = CubicInterpolate(GetClamped(i-1, row), GetClamped(i, row), GetClamped(i+1, row), GetClamped(i+2, row), tx);
		
		// interpoalte vertically
		float r = CubicInterpolate(a, b, c, d, ty);
		
		/*
		float minimum = FLT_MAX;
		float maximum = 0.0f;
		
		for (int x=i; x < i+2; ++x)
		{
			for (int y=j; y < j+2; ++y)
			{
				minimum = std::min(minimum, GetClamped(x, y));
				maximum = std::max(maximum, GetClamped(x, y));
			}
		}
		minimum = std::max(0.0f, minimum);
		
		r = Clamp(r, minimum, maximum);
		*/

		return r;		
	}
	
	float Get(int x, int y) const
	{
		if (x < 0 || x >= m_width || y < 0 || y >= m_height)
		{
			return 0.0f;
		}
		else
		{
			return m_data[y*m_width+x];
		}
	}
	
	float GetClamped(int x, int y) const
	{
		//return Get(x, y);
		
		x = std::max(0, std::min(x, m_width-1));
		y = std::max(0, std::min(y, m_height-1));
		
		return m_data[y*m_width+x];
	}
	
	void Set(int x, int y, float v)
	{
		assert(x >= 0 && x < m_width);
		assert(y >= 0 && y < m_height);

		m_data[y*m_width+x] = v;
	}

	void Add(int x, int y, float v)
	{
		assert(x >= 0 && x < m_width);
		assert(y >= 0 && y < m_height);
		
		m_data[y*m_width+x] += v;
	}
	
	// useful for debugging
	float Sum() const
	{
		float s = 0.0f;

		for (int i=0; i < m_width*m_height; ++i)
		{
			s += fabsf(m_data[i]);
		}

		return s;
	}

	void Reset()
	{
		memset(m_data, 0, sizeof(float)*m_width*m_height);
	}

	float* m_data;

	int m_width;
	int m_height;
};


template <bool cubic>
void AdvectQ(Grid2D& newq, const Grid2D& q, const Grid2D& u, const Grid2D& v, float dt)
{
	ScopedTimer timer("AdvectQ");
	
	// for each cell trace backwards and interpolate the quantity
	for (int y=0; y < q.m_height; ++y)
	{
		for (int x=0; x < q.m_width; ++x)
		{
			float px = x - dt*u.LinearInterp(x+0.5f, y);
			float py = y - dt*v.LinearInterp(x, y+0.5f);

			// interpolate the quantity at p 
			if (cubic)
			{
				newq.Set(x, y, q.CubicInterp(px, py));
			}
			else
			{
				newq.Set(x, y, q.LinearInterp(px, py));
			}
		}
	}
}


void AdvectU(Grid2D& newu, const Grid2D& u, const Grid2D& v, float dt)
{
	// for each cell trace backwards and interpolate the quantity
	for (int y=0; y < kHeight; ++y)
	{
		for (int x=0; x < kWidth+1; ++x)
		{
			float iu = u.Get(x, y);
			float iv = v.LinearInterp(x-0.5f, y+0.5f);
		
			float px = x - dt*iu;
			float py = y - dt*iv;
			
			// interpolate the quantity at p 
			newu.Set(x, y, u.LinearInterp(px, py));
		}
	}
}

void AdvectV(Grid2D& newv, const Grid2D& u, const Grid2D& v, float dt)
{	
	// for each cell trace backwards and interpolate the quantity
	for (int y=0; y < kHeight+1; ++y)
	{
		for (int x=0; x < kWidth; ++x)
		{
			float iu = u.LinearInterp(x+0.5f, y-0.5f);
			float iv = v.Get(x, y);
			
			float px = x - dt*iu;
			float py = y - dt*iv;
			
			// interpolate the quantity at p 
			newv.Set(x, y, v.LinearInterp(px, py));
		}
	}
}

Grid2D g_u1(kWidth+1, kHeight);
Grid2D g_v1(kWidth, kHeight+1);

Grid2D g_u2(kWidth+1, kHeight);
Grid2D g_v2(kWidth, kHeight+1);

Grid2D g_t1(kWidth, kHeight);
Grid2D g_t2(kWidth, kHeight);

Grid2D g_s1(kWidth, kHeight);
Grid2D g_s2(kWidth, kHeight);

Grid2D* g_currentU = &g_u1;
Grid2D* g_currentV = &g_v1;

Grid2D* g_lastU = &g_u2;
Grid2D* g_lastV = &g_v2;

Grid2D* g_currentT = &g_t1;
Grid2D* g_lastT = &g_t2;

Grid2D* g_currentS = &g_s1;
Grid2D* g_lastS = &g_s2;

Grid2D g_divergence(kWidth, kHeight);
Grid2D g_pressure(kWidth, kHeight);

struct Particle
{
	Vec2 x;
	float w;
};

std::vector<Particle> g_particles;

float g_screen[kHeight][kWidth][4];

float g_blackBodyMinT = 1000.0f;
float g_blackBodyMaxT = 4000.0f;

const uint32_t kBlackBodySamples = 1024;
Colour g_blackBodyTable[kBlackBodySamples];


void Reset()
{
	g_u1.Reset();
	g_u2.Reset();
	g_v1.Reset();
	g_v2.Reset();
	g_t1.Reset();
	g_t2.Reset();
	g_s1.Reset();
	g_s2.Reset();
	g_divergence.Reset();
	g_pressure.Reset();

	g_particles.resize(0);

	
}

void Init()
{
	Reset();
	
	for (int i=0 ; i < kBlackBodySamples; ++i)
	{
		Colour c;
		BlackBodyXYZ(g_blackBodyMinT + (g_blackBodyMaxT-g_blackBodyMinT)*(float(i)/kBlackBodySamples), (float*)&c);

		c = XYZToLinear(c.r, c.g, c.b);
		c.r = Max(0.0f, c.r);
		c.g = Max(0.0f, c.g);
		c.b = Max(0.0f, c.b);
		
		//float m = Max(c.r, Max(c.g, c.b));		
		//c /= m;
		
		g_blackBodyTable[i] = c;
	}
}

void Divergence(Grid2D& out, const Grid2D& u, const Grid2D& v)
{
	for (int y=0; y < kHeight; ++y)
	{
		for (int x=0; x < kWidth; ++x)
		{
			//out.Set(x, y, (u.dX(x, y) + v.dY(x,y)));
			float du = u.Get(x+1, y)-u.Get(x, y);
			float dv = v.Get(x, y+1)-v.Get(x, y);
			
			out.Set(x, y, du+dv);
		}
	}
}

const float dt = 1.0f / 60.0f;

void PressureSolve(Grid2D& out, const Grid2D& divergence)
{
	Grid2D scratch(kWidth, kHeight);

	Grid2D* read = &out;
	Grid2D* write = &scratch;
	
	const float scale = dt;
	const float rcpScale = 1.0f / scale;
	
	for (int i=0; i < 20; ++i)
	{
		for (int y=0; y < kHeight; ++y)
		{
			for (int x=0; x < kWidth; ++x)
			{
				// Jacobi solver
				float p = -divergence.Get(x, y);

				float a = read->Get(x-1, y);
				float b = read->Get(x+1, y);
				float c = read->Get(x, y+1);
				float d = read->Get(x, y-1);
				
				float v = 0.25f*(p + scale*(a+b+c+d))*rcpScale;
				//float v = (p+a+b+c+d)/4.0f;

				//if (!isfinite(v) || isnan(v))
				//   printf("blah");

				write->Set(x, y, v);
			}
		}

		std::swap(read, write);
	}
}

void PressureApply(Grid2D& u, Grid2D& v, const Grid2D& pressure)
{
	ScopedTimer timer("PressureApply");
	
	/*
	for (int y=0; y < kHeight; ++y)
	{
		for (int x=0; x < kWidth; ++x)
		{
			float p = dt*pressure.Get(x, y);
			
			u.Add(x, y, -p);
			u.Add(x+1, y, p);
			
			v.Add(x, y, -p);
			v.Add(x, y+1, p);
		}
	}
	*/
	
	for (int y=0; y < kHeight; ++y)
	{
		for (int x=0; x < kWidth+1; ++x)
		{
			// calculate pressure differential for this cell
			float pu = pressure.Get(x, y)-pressure.Get(x-1, y);
													   
			float newu = u.Get(x, y);
			newu -= dt*pu;            
			u.Set(x, y, newu);
		}
	}
	
	for (int y=0; y < kHeight+1; ++y)
	{
		for (int x=0; x < kWidth; ++x)
		{
			// calculate pressure differential for this cell
			float pv = pressure.Get(x, y)-pressure.Get(x, y-1);
						
			float newv = v.Get(x, y);
			newv -= dt*pv;
			v.Set(x, y, newv);
		}
	}
}

void UpdateParticles(Grid2D& newu, Grid2D& newv, const Grid2D& u, const Grid2D& v, float dt)
{
	std::vector<Particle> newParticles;
	newParticles.reserve(g_particles.size());

	// advect particles
	for (size_t i=0; i < g_particles.size(); ++i)
	{
		Particle& p = g_particles[i];

		float ux = u.LinearInterp(p.x.x, p.x.y);
		float uy = v.LinearInterp(p.x.x, p.x.y);

		p.x.x += dt*ux;
		p.x.y += dt*uy;

		if (p.x.x < 0.0f || p.x.x > kWidth || p.x.y < 0.0f || p.x.y > kHeight)
		{
		}
		else
			newParticles.push_back(p);
	}

	g_particles.swap(newParticles);

	// 
	const float kEpsilon = 1.0f*dt;

	// copy vorticity to the grid u values
	for (size_t i=0; i < g_particles.size(); ++i)
	{
		Particle& p = g_particles[i];

		// apply to u
		{
			int lx = max(0.0f, floorf(p.x.x)-1);
			int ux = min(kWidth, lx+3);

			int ly = max(0.0f, floorf(p.x.y)-1);
			int uy = min(kHeight-1, ly+2);

			for (int y=ly; y <= uy; ++y)
			{
				for (int x=lx; x <= ux; ++x)
				{
					float cx = x;
					float cy = y + 0.5f;

					float dx = p.x.x - cx;
					float dy = p.x.y - cy;

					float s = kEpsilon*max(0.0f, 1.0f - 0.5f*sqrtf(dx*dx + dy*dy));

					Vec2 f = PerpCCW(Normalize(Vec2(dx, dy)));

					newu.Add(x, y, f.x*s*p.w);					
				}
			}			
		}

		// apply to v
		{
			int lx = max(0.0f, floorf(p.x.x)-1);
			int ux = min(kWidth-1, lx+2);

			int ly = max(0.0f, floorf(p.x.y)-1);
			int uy = min(kHeight, ly+3);

			for (int y=ly; y <= uy; ++y)
			{
				for (int x=lx; x <= ux; ++x)
				{
					float cx = x + 0.5f;
					float cy = y;

					float dx = p.x.x - cx;
					float dy = p.x.y - cy;

					float s = kEpsilon*max(0.0f, 1.0f - 0.5f*sqrtf(dx*dx + dy*dy));

					Vec2 f = PerpCCW(Normalize(Vec2(dx, dy)));

					newv.Add(x, y, f.y*s*p.w);	
				}
			}			
		}

	}
}

void InjectVorticity(Grid2D& newu, Grid2D& newv, const Grid2D& u, const Grid2D& v, float dt)
{
	Grid2D curl(kWidth, kHeight);

	for (int y=0; y < kHeight; ++y)
	{
		for (int x=0; x < kWidth; ++x)
		{
			float dvdx = v.Get(x+1, y)-v.Get(x, y);
			float dudy = u.Get(x, y+1)-u.Get(x, y);
			
			curl.Set(x, y, dvdx-dudy);
		}
	}
	
	for (int y=0; y < kHeight; ++y)
	{
		for (int x=0; x < kWidth; ++x)
		{
			float c = curl.Get(x, y);
			float e = 2.0f;//2.0f;//1.5f;
			
			float dcdx = 0.5f*(fabsf(curl.Get(x+1, y))-fabsf(curl.Get(x-1, y)));
			float dcdy = 0.5f*(fabsf(curl.Get(x, y+1))-fabsf(curl.Get(x, y-1)));
		
			float l = sqrt(dcdx*dcdx + dcdy*dcdy) + 0.001f;
			
			dcdx = dt*e*(dcdx / l);
			dcdy = dt*e*(dcdy / l);
			
			newu.Add(x, y, c*dcdy);
			newv.Add(x, y, -dcdx*c);
		}
	}
}

void ForcesApply(Grid2D& u, Grid2D& v)
{
   for (int y=0; y < kHeight; ++y)
	{
		for (int x=0; x < kWidth; ++x)
		{
			float s = g_lastS->Get(x, y);
			float t = g_lastT->Get(x, y)*50.0f;
			
			float b = dt*(-1.0*s + t); 
			
			v.Add(x, y, b*0.5f);
			v.Add(x, y+1, b*0.5f);
		}
	}
}

void ApplyBoundaries(Grid2D& u, Grid2D& v)
{
	for (int y=0; y < kHeight; ++y)
	{
		for (int x=0; x < kWidth; ++x)
		{
			if (x == 0 || x == kWidth-1)
				u.Set(x, y, 0.0f);
			if (y == 0 || y == kHeight-1)
				v.Set(x, y, 0.0f);
		}
	}
}

bool g_explode = false;

void InitSolid(float dt)
{
	/*
	// initialize the velocity field
	for (int y=0; y < 10; ++y)
	{
		for (int x=0; x < 11; ++x)
		{
			if (g_explode)
			{
				g_lastV->Set(kWidth/2-5 + x, y, 20.0f);
			}

			if (y == 0)
			{
				g_lastT->Set(kWidth/2-5 + x, y, 1.0f);
			}
		}
	}    
	*/
	
	static float t=0.0f;
	if (g_explode)
		t = 0.0f;
	
	for (int y=0; y < 6; ++y)
	{
		for (int x=0; x < 6; ++x)
		{
			int ix = kWidth/2 + x - 3;

			g_lastS->Set(ix, y, 1.0f);		
			//g_lastT->Set(ix, y, 0.5f);
			
			if (t == 0.0f)
			{
				g_lastT->Set(ix, y, Randf(0.0, 0.8f));
				g_divergence.Add(ix, y, -500.0f*Randf());

			}

			if (t < 0.2f)
			{
				for (int i=0; i < 50; ++i)
				{
					// add a vortex particle
					//if (Randf() > 0.95f)
					{
						Particle p;
						p.w = Randf(-1.0f, 1.0f);
						p.x = Vec2(ix, y) + Vec2(Randf(-0.5f, 0.5f), Randf(-0.5f, 0.5f));

						g_particles.push_back(p);
					}
				}
			}
			//g_lastV->Set(ix, y, 200.0f);

		}
	}

	if (t < 0.5)
	{
		//g_divergence.Add(kWidth/2, 0, -200.0f);
		//g_divergence.Add(kWidth/2+1, 0, -500.0f);
	}
	
	t += dt;
}
bool g_drawVelocity = false;
bool g_drawParticles = false;


void GLUTUpdate()
{
	ForcesApply(*g_lastU, *g_lastV);
	
	// calculate divergence
	Divergence(g_divergence, *g_lastU, *g_lastV);

	InitSolid(dt);	
	
	// solve for the pressure whose gradient when applied to the current velocity will yield a field with zero divergence
	PressureSolve(g_pressure, g_divergence);
	PressureApply(*g_lastU, *g_lastV, g_pressure);

	ApplyBoundaries(*g_lastU, *g_lastV);
	
	/*
	
	Grid2D scratch(kWidth, kHeight);
	Divergence(scratch, *g_currentU, *g_currentV);

	// subtract pressure gradient
	printf("%f : %f\n", g_divergence.Sum(), scratch.Sum());
	*/

	// advect temperature
	AdvectQ<false>(*g_currentT, *g_lastT, *g_lastU, *g_lastV, dt);
	AdvectQ<true>(*g_currentS, *g_lastS, *g_lastU, *g_lastV, dt); 
				  
	// advect velocity
	AdvectU(*g_currentU, *g_lastU, *g_lastV, dt);
	AdvectV(*g_currentV, *g_lastU, *g_lastV, dt);

	UpdateParticles(*g_currentU, *g_currentV, *g_lastU, *g_lastV, dt);

	//InjectVorticity(*g_currentU, *g_currentV, *g_lastU, *g_lastV, dt);
	
	// add forces

	// update screen
	//g_currentT->Set(64, 64, 15.0f);
//	g_currentT->Set(64, 62, 1.0f);
//	g_currentT->Set(62, 64, 1.0f);
//	g_currentT->Set(62, 62, 1.0f);
//	
	static float maxI = 1.0f;
	float newMaxI = 0.0f;
	
	if (!g_drawVelocity)
	{
		for (int y=0; y < kHeight; ++y)
		{
			for (int x=0; x < kWidth; ++x)
			{
				//float t = powf(fabsf(g_currentT->Get(x, y)), 1.0f/1.0f);//fabsf(g_currentT->CubicInterp(60 + x*0.05f, 60.0f + y*0.05f));
				uint32_t s = g_currentT->Get(x, y)*kBlackBodySamples;
				
				Colour c = g_blackBodyTable[Clamp(s, 0U, kBlackBodySamples-1)];
				newMaxI = Max(g_blackBodyTable[300].r, Max(newMaxI, Max(c.r, Max(c.g, c.b))));
				c /=  maxI*0.01f;
				
				//float b = (x/16+y/16)&1?0.1f:0.3f;
				float b = 0.0f;
				
				Colour f = Lerp(Colour(b, b, b), c, g_currentS->Get(x, y));
				f = LinearToSrgb(f);
				
				/*
				g_screen[y][x][0] = f.r;//fabsf(g_currentU->Get(x, y)*0.05f);
				g_screen[y][x][1] = f.g;//fabsf(g_currentV->Get(x, y)*0.05f);
				g_screen[y][x][2] = f.b;//abs(g_divergence.Get(x, y));//0.0f;//g_currentT->Get(x, y);//0.0f;
				g_screen[y][x][3] = 1.0f;
				*/

				g_screen[y][x][0] = g_currentS->Get(x, y);
				g_screen[y][x][1] = g_currentS->Get(x, y);
				g_screen[y][x][2] = g_currentS->Get(x, y);
				g_screen[y][x][3] = 1.0f;
			}
		}
		
		maxI = newMaxI;// Lerp(maxI, newMaxI, 0.1f);
	}
	else 
	{
		for (int y=0; y < kHeight; ++y)
		{
			for (int x=0; x < kWidth; ++x)
			{
				g_screen[y][x][0] = fabsf(g_currentU->Get(x, y)*0.05f);
				g_screen[y][x][1] = fabsf(g_currentV->Get(x, y)*0.05f);
				g_screen[y][x][2] = g_currentS->Get(x, y);//abs(g_divergence.Get(x, y));//0.0f;//g_currentT->Get(x, y);//0.0f;
				g_screen[y][x][3] = 1.0f;
			}
		}
	}

		
	glDisable(GL_BLEND);
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	glPixelZoom(g_zoomx, g_zoomy);
	glDrawPixels(kWidth,kHeight,GL_RGBA,GL_FLOAT, g_screen);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0f, kWidth, 0.0f, kHeight);

	
	if (g_drawParticles)
	{
	
		glBegin(GL_POINTS);
		glColor3f(0.0f, 0.0f, 1.0f);

		for (size_t i=0; i < g_particles.size(); ++i)
		{
			glVertex2fv(g_particles[i].x);
		}
	
		glEnd();
	}

	// flip
	glutSwapBuffers();

	// swap physics buffers
	std::swap(g_currentU, g_lastU);
	std::swap(g_currentV, g_lastV);
	std::swap(g_currentT, g_lastT);
	std::swap(g_currentS, g_lastS);

	
}

void GLUTReshape(int width, int height)
{
	g_zoomx = float(width) / kWidth;
	g_zoomy = float(height) / kHeight;
}

void GLUTArrowKeys(int key, int x, int y)
{
}

void GLUTArrowKeysUp(int key, int x, int y)
{
}

void GLUTKeyboardDown(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'e':
		{
			g_explode = true;
			break;
		}
		case 'r':
		{
			Reset();
			break;
		}
		case 't':
		{
			g_drawVelocity = !g_drawVelocity;
			break;
		}
		case 'p':
		{
			g_drawParticles = !g_drawParticles;
			break;
		}
		
	case 27:
		exit(0);
		break;
	};
}

void GLUTKeyboardUp(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'e':
		{
			g_explode = false;
			break;
		}
	}
}

static int lastx;
static int lasty;

void GLUTMouseFunc(int b, int state, int x, int y)
{
	y = (kHeight*g_zoomy)-y;

	switch (state)
	{
	case GLUT_UP:
		{
			lastx = x;
			lasty = y;			
		}
	case GLUT_DOWN:
		{
			lastx = x;
			lasty = y;
		}
	}
}

void GLUTMotionFunc(int x, int y)
{
	y = (kHeight*g_zoomy)-y;

	int dx = x-lastx;
	int dy = y-lasty;

	lastx = x;
	lasty = y;
	
	int i = x/g_zoomx;
	int j = y/g_zoomy;
	
	const int width = 2;
	
	for (int x=i-width; x < i+width; ++x)
	{
		for (int y=j-width; y < j+width; ++y)
		{
			if (x > 0 && x < kWidth && y > 0 && y < kHeight)
			{
				//g_lastU->Set(x, y, dx*100);
				//g_lastV->Set(x, y, dy*10);
				g_lastT->Set(x, y, Randf(0.5f, 0.6f));
				g_lastS->Set(x, y, 1.0f);//Randf(0.0f, 1.0f));
			}
		}
	}
}


int main(int argc, char* argv[])
{	
	RandInit();
	
	// init gl
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH);

	glutInitWindowSize(kWidth*4, kHeight*4);
	glutCreateWindow("Fluid2D");
	glutPositionWindow(200, 200);

	Init();

	glutMouseFunc(GLUTMouseFunc);
	glutReshapeFunc(GLUTReshape);
	glutDisplayFunc(GLUTUpdate);
	glutKeyboardFunc(GLUTKeyboardDown);
	glutKeyboardUpFunc(GLUTKeyboardUp);
	glutIdleFunc(GLUTUpdate);	
	glutSpecialFunc(GLUTArrowKeys);
	glutSpecialUpFunc(GLUTArrowKeysUp);
	glutMotionFunc(GLUTMotionFunc);

	glutMainLoop();
}


