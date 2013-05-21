/*
 *  Fluid3D.cpp
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

using namespace std;

#include "core/perlin.h"
#include "core/platform.h"
#include "core/maths.h"
#include "core/tga.h"
#include "core/shader.h"

#include "grid3d.h"
#include "blackbody.h"

const int kWidth = 128;
const int kHeight = 256;
const int kDepth = 128;

typedef XGrid3D<kWidth, kHeight, kDepth> Grid3D;

GLuint CreateVolumeTexture(int width, int height, int depth);
void UpdateVolumeTexture(GLuint texid, const Grid3D& g, float scale=1.0f);
void UpdateSmokeTexture(GLuint texid, const Grid3D& a, Grid3D& b);

GLuint g_noiseTexture;
GLuint g_densityTexture;
GLuint g_temperatureTexture;
GLuint g_blackBodyTexture;
GLuint g_volumeShader;
GLuint g_shadowShader;

int g_screenWidth = 512;
int g_screenHeight = 1024;

float g_xoff = 10.0f;
float g_yoff = 5.0f;
float g_zoff = 100.0f;

float g_angle = 0.0f;
float g_angularVelocity = 0.0f;

float g_lightTheta = DegToRad(60.0f);
float g_lightPhi = DegToRad(-20.0f);
float g_lightRadius = 3.0f;
float g_lightIntensity = 20.0f;

float g_absorption = 30.5f;
float g_scatter = 15.5f;

float g_minT = 500.0f;
float g_maxT = 4000.0f;

bool g_pause = false;
bool g_drawParticles = false;

struct Particle
{
	Vec3 w;
	Vec3 p;
};

std::vector<Particle> g_particles;

Grid3D g_u1;
Grid3D g_v1;
Grid3D g_w1;

Grid3D g_u2;
Grid3D g_v2;
Grid3D g_w2;

Grid3D g_t1;
Grid3D g_t2;

Grid3D g_s1;
Grid3D g_s2;

Grid3D g_f1;
Grid3D g_f2;

Grid3D* g_currentU = &g_u1;
Grid3D* g_currentV = &g_v1;
Grid3D* g_currentW = &g_w1;

Grid3D* g_lastU = &g_u2;
Grid3D* g_lastV = &g_v2;
Grid3D* g_lastW = &g_w2;

Grid3D* g_currentS = &g_s1;
Grid3D* g_lastS = &g_s2;

Grid3D* g_currentT = &g_t1;
Grid3D* g_lastT = &g_t2;

Grid3D* g_currentF = &g_f1;
Grid3D* g_lastF = &g_f2;

Grid3D g_divergence;
Grid3D g_pressure;


void FluidApplyForces(float dt)
{
	ScopedTimer timer(__FUNCTION__);
	
	for (int z=0; z < kDepth; ++z)
	{	
		for (int y=0; y < kHeight; ++y)
		{
			for(int x=0; x < kWidth; x++)
			{
				float s = g_lastF->Get(x, y, z);
				float t = g_lastT->Get(x, y, z)*0.135f;
				
				float b = dt*(-40.0f*g_lastS->Get(x, y, z) + t); 
				
				g_lastV->Add(x, y, z, b*1.5f);				
			}
		}
	}
}

void FluidApplyBoundaries()
{
	ScopedTimer timer(__FUNCTION__);

	for (int z=0; z < kDepth; ++z)
	{	
		for (int y=0; y < kHeight; ++y)
		{
			for(int x=0; x < kWidth; x++)
			{				
				// boundary conditions
				if (z == 0 || z == kDepth-1)
				{
					g_lastW->Set(x, y, z, 0.0f);
				}
				if (y == 0 || y == kHeight-1)
				{
					g_lastV->Set(x, y, z, 0.0f);
				}
				if (x == 0 || x == kWidth-1)
				{
					g_lastU->Set(x, y, z, 0.0f);
				}		
			}
		}
	}
}

Vector3 GetCurl(const Grid3D& u, const Grid3D& v, const Grid3D& w, int x, int y, int z)
{
	float dwdy = 0.5f*(w.Get(x, y+1, z)-w.Get(x, y-1, z));
	float dwdx = 0.5f*(w.Get(x+1, y, z)-w.Get(x-1, y, z));
	
	float dudz = 0.5f*(u.Get(x, y, z+1)-u.Get(x, y, z-1));
	float dudy = 0.5f*(u.Get(x, y+1, z)-u.Get(x, y-1, z));
	
	float dvdz = 0.5f*(v.Get(x, y, z+1)-v.Get(x, y, z-1));
	float dvdx = 0.5f*(v.Get(x+1, y, z)-v.Get(x-1, y, z));
	
	return Vector3(dwdy-dvdz, dudz-dwdx, dvdx-dudy);
}

void FluidInjectVorticity(Grid3D& newu, Grid3D& newv, Grid3D& neww, const Grid3D& u, const Grid3D& v, const Grid3D& w, float dt)
{	
	ScopedTimer timer(__FUNCTION__);
	
	for (int z=0; z < kDepth; ++z)
	{
		for (int y=0; y < kHeight; ++y)
		{
			for (int x=0; x < kWidth; ++x)
			{
				// calculate gradient of curl magnitude
				float dcdx = 0.5f*(Length(GetCurl(u, v, w, x+1, y, z)) - Length(GetCurl(u, v, w, x-1, y, z)));
				float dcdy = 0.5f*(Length(GetCurl(u, v, w, x, y+1, z)) - Length(GetCurl(u, v, w, x, y-1, z)));
				float dcdz = 0.5f*(Length(GetCurl(u, v, w, x, y, z+1)) - Length(GetCurl(u, v, w, x, y, z-1)));
				
				
				Vector3 n = SafeNormalize(Vector3(dcdx, dcdy, dcdz));
				
				float e = 55.25f;
				
				Vector3 f = e*dt*Cross(n, GetCurl(u, v, w, x, y, z));
				
				newu.Add(x, y, z, f.x);
				newv.Add(x, y, z, f.y);
				neww.Add(x, y, z, f.z);
			}
		}
	}
}

void FluidCombust(Grid3D& fuel, Grid3D& divergence, Grid3D& smoke, Grid3D& temperature, float dt)
{
	dt = max(1.0f/60.0f, dt);

	const float kIgnition = 300;
	const float kBurnSpeed = 0.25f;
	const float kSmoke = 1.2f;
	const float kHeat = 3500.0f;
	const float kExpansion = 10.5f / dt;
	const float kDecay = expf(-0.4f * dt);
	const float kSmokeDecay = expf(-0.5f * dt);
	
	for (int z=0; z < kDepth; ++z)
	{
		for (int y=0; y < kHeight; ++y)
		{
			for (int x=0; x < kWidth; ++x)
			{
				float f = fuel.Get(x, y, z);
				float t = temperature.Get(x, y, z);
				float s = smoke.Get(x, y, z);
				
				if (f > 0.0f && t >= kIgnition)
				{
					float nf = max(0.0f, f-kBurnSpeed*dt);
					float df = f - nf;
					
					// update fuel
					fuel.Set(x, y, z, nf);
					
					// emitted smoke density starts low and increases as combustion occurs
					float smokeRamp = 1.0f-nf;
					
					// other side effects
					s += df*kSmoke*smokeRamp;
					
					// generate heat
					t += df*kHeat;
					
					//divergence.Add(x, y, z, -df*kExpansion);					
				}
				
				temperature.Set(x, y, z, t*kDecay);
				smoke.Set(x, y, z, s*kSmokeDecay);
			}
		}
	}
}

Vec3 grad(const Grid3D& u, Vec3 p)
{
	float h = 1.0f;

	Vec3 r(u.LinearInterp(p.x + h, p.y, p.z) - u.LinearInterp(p.x - h, p.y, p.z),
	 	   u.LinearInterp(p.x, p.y + h, p.z) - u.LinearInterp(p.x, p.y - h, p.z),
		   u.LinearInterp(p.x, p.y, p.z + h) - u.LinearInterp(p.x, p.y, p.z - h));

	return r*0.5f;
}

float gauss(Vec3 delta, float radius)
{
	return expf(-Dot(delta, delta) / (2.0f*radius*radius)) / (radius*radius*radius*15.74960f);
}

void FluidUpdateParticles(Grid3D& newu, Grid3D& newv, Grid3D& neww, const Grid3D& u, const Grid3D& v, const Grid3D& w, float dt)
{
	ScopedTimer timer("UpdateParticles");

	std::vector<Particle> newParticles;
	newParticles.reserve(g_particles.size());
	
	for (size_t i=0; i < g_particles.size(); ++i)
	{
		Particle& p = g_particles[i];

		p.p.x += u.LinearInterp(p.p.x, p.p.y, p.p.z)*dt;
		p.p.y += v.LinearInterp(p.p.x, p.p.y, p.p.z)*dt;
		p.p.z += w.LinearInterp(p.p.x, p.p.y, p.p.z)*dt;

		if (p.p.x <= 0.0f || p.p.x >= kWidth  ||
			p.p.y <= 0.0f || p.p.y >= kHeight ||
			p.p.z <= 0.0f || p.p.z >= kDepth)
			continue;

		// calculate gradient of velocity field
		Vec3 gradu = grad(u, p.p);
		Vec3 gradv = grad(v, p.p);
		Vec3 gradw = grad(w, p.p);

		// rotate particle vorticity vector
		Vec3 dw = p.w.x * gradu + p.w.y * gradv + p.w.z * gradw;
		//Vec3 w = p.w.x * (gradu.x + gradv.x + gradw.x) + p.w.y * (gradu.y + gradv.y + gradw.y) + p.w.z * (gradu.z + gradv.z + gradw.z);
		
		p.w = Normalize(p.w + dt*dw);
		
		newParticles.push_back(p);

		const int kSupport = 4;

		int lz = max(0, int(floorf(p.p.z))-kSupport);
		int uz = min(kDepth, lz+kSupport);

		int ly = max(0, int(floorf(p.p.y))-kSupport);
		int uy = min(kHeight, ly+kSupport);

		int lx = max(0, int(floorf(p.p.x))-kSupport);
		int ux = min(kWidth, lx+kSupport);

		const float kEpsilon = 300000.0f*dt*3.0f;

		for (int z=lz; z <= uz; ++z)
		{
			for (int y=ly; y <= uy; ++y)
			{
				for (int x=lx; x <= ux; ++x)
				{
					Vec3 delta = p.p - Vec3(x, y, z);
					
					float l = Length(delta);

					if (l <= 4.0f && l > 0.0f)
					{
						float s = gauss(delta, 4.0f);

						Vec3 f = Cross(p.w, delta/l)*kEpsilon*s;

						newu.Add(x, y, z, f.x);
						newv.Add(x, y, z, f.y);
						neww.Add(x, y, z, f.z);
					}
				}
			}
		}
	}

	g_particles.swap(newParticles);
}

void FluidInit(float dt)
{
	static float t = 0.0f;
	
	const int kSeedSize = 24;

	for (int y=0; y < kSeedSize; ++y)
	{
		for (int z=0; z < kSeedSize; ++z)
		{
			for (int x=0; x < kSeedSize; ++x)
			{	
				int ix = kWidth/2 + x - kSeedSize/2;
				int iy = y;
				int iz = kDepth/2 + z - kSeedSize/2;
				
				if (t < 0.05f)
				{					 
					g_divergence.Add(ix, iy, iz, -150.0f*(Randf()+0.1f)*(0.5f-t)*2.0f);
				}

				if (t < 0.5f)
				{
					// add particles
					for (int i=0; i < 1; ++i)
					{
						Particle p;
						p.p = Vec3(ix, iy, iz) + Vec3(Randf(-2.0f, 2.0f), Randf(0.0f, 1.0f), Randf(-2.0f, 2.0f));
						
						float dx = kWidth/2.0f - p.p.x;
						float dz = kDepth/2.0f - p.p.z;

						p.w = Normalize(Vec3(dz, 0.0f, -dx));

						if (Randf() > 0.95f)
							g_particles.push_back(p);
					}

					g_lastS->Set(ix, iy, iz, 1.0f);	
					g_lastT->Add(ix, iy, iz, (t*2.0f + 0.5f)*200.0f);

				}
			}
		}
	}
	
	t += dt;
}

void FluidInit2(float dt)
{
	static float t = 0.0f;
	
	if (t >= 3.0f)
	{
		return;
	}

	t += dt;

	const int kSeedSize = 3;
		
	for (int x=0; x < 6; ++x)
	{
		for (int y=0; y < kSeedSize; ++y)
		{
			for (int z=0; z < kSeedSize; ++z)
			{
				int ix = x;
				int iy = kHeight/2 - kSeedSize/2 + y;
				int iz = kDepth/2 + z - kSeedSize/2;
				
				g_lastF->Add(ix, iy, iz, 1.0f*Randf(0.5f, 1.5f));
			}
		}
	}

	for (int x=0; x < 6; ++x)
	{
		for (int y=-1; y < kSeedSize+1; ++y)
		{
			for (int z=-1; z < kSeedSize+1; ++z)
			{
				int ix = x;
				int iy = kHeight/2 - kSeedSize/2 + y;
				int iz = kDepth/2 + z - kSeedSize/2;

				
				g_lastT->Add(ix, iy, iz, 4000.0f*Randf(0.5f, 1.5f));

				if (y == -1 || z == -1 || y == kSeedSize || z == kSeedSize)
				{
					g_lastU->Set(ix, iy, iz, 500.0f*Randf(0.5f, 1.5f));
					//g_lastV->Set(ix, iy, iz, -100.0f*Randf(0.5f, 1.0f));
				}
				else
				{
					g_lastU->Set(ix, iy, iz, 800.0f*Randf(0.5f, 1.5f));
					//g_lastV->Set(ix, iy, iz, -100.0f*Randf(0.5f, 1.0f));
				}
			}
		}
	}
}


void FluidStep(float dt)
{
	
	if (!g_pause)
	{
		FluidApplyForces(dt);

		CalculateDivergence(g_divergence, *g_lastU, *g_lastV, *g_lastW, dt);

		FluidInit(dt);		
		//FluidCombust(*g_lastF, g_divergence, *g_lastS, *g_lastT, dt);
					 
		//std::cout << "Divergance before: " << g_divergence.Sum() << std::endl;
		
		PressureSolve(g_pressure, g_divergence, dt);
		PressureApply(*g_lastU, *g_lastV, *g_lastW, g_pressure, dt);

		FluidApplyBoundaries();

		/*
		Grid3D scratch(kWidth, kHeight, kDepth);
		CalculateDivergence(scratch, *g_lastU, *g_lastV, *g_lastW, dt);
		
		std::cout << "Divergance after: " << scratch.Sum() << std::endl;
		*/
		
		{ ScopedTimer timer("Advection");
		
			// advect velocity
			AdvectQ<false>(*g_currentU, *g_lastU, *g_lastU, *g_lastV, *g_lastW, dt);
			AdvectQ<false>(*g_currentV, *g_lastV, *g_lastU, *g_lastV, *g_lastW, dt);
			AdvectQ<false>(*g_currentW, *g_lastW, *g_lastU, *g_lastV, *g_lastW, dt);
			
			// advect temperature
			AdvectQ<false>(*g_currentT, *g_lastT, *g_lastU, *g_lastV, *g_lastW, dt);
			
			// advect smoke (using cubic interpolation)
			AdvectQ<true>(*g_currentS, *g_lastS, *g_lastU, *g_lastV, *g_lastW, dt);			
			
			// advect fuel
			//AdvectQ<true>(*g_currentF, *g_lastF, *g_lastU, *g_lastV, *g_lastW, dt);
		}
			
		//FluidInjectVorticity(*g_currentU, *g_currentV, *g_currentW, *g_lastU, *g_lastV, *g_lastW, dt);
		FluidUpdateParticles(*g_currentU, *g_currentV, *g_currentW, *g_lastU, *g_lastV, *g_lastW, dt);
	}
		
	cout << "Average T: " << g_currentT->Average() << endl;
	
	// swap buffers
	std::swap(g_currentU, g_lastU);
	std::swap(g_currentV, g_lastV);
	std::swap(g_currentW, g_lastW);
	std::swap(g_currentS, g_lastS);
	std::swap(g_currentT, g_lastT);
	std::swap(g_currentF, g_lastF);
}


GLuint CreateBlackBodyTexture(float minT, float maxT)
{
	const uint32_t kWidth = 1024;
	
	Colour* data = new Colour[kWidth];
	
	for (uint32_t i=0; i < kWidth; ++i)
	{
		float T = minT + (maxT-minT)*(float(i)/kWidth);
		
		Colour c;
		BlackBodyXYZ(T, (float*)&c);
		
		c = XYZToLinear(c.r, c.g, c.b);
		
		c.r = max(0.0f, c.r);
		c.g = max(0.0f, c.g);
		c.b = max(0.0f, c.b);
		
		data[i] = c;
	}
		
	GLuint texid;
	glGenTextures(1, &texid);
	
	glVerify(glBindTexture(GL_TEXTURE_2D, texid));
	
	glVerify(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR));
	glVerify(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR));
	
	glVerify(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE));
	glVerify(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE));
	glVerify(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE));
	
	glVerify(glTexImage2D(GL_TEXTURE_2D,
				 0,
				 GL_RGBA32F_ARB,
				 kWidth,
				 1,
				 0,
				 GL_RGBA,
				 GL_FLOAT,
				data));
	
	delete[] data;
	
	return texid;
}


GLuint CreateVolumeTexture(int width, int height, int depth)
{
	GLuint texid;
	glGenTextures(1, &texid);
	GLenum target = GL_TEXTURE_3D;

	glVerify(glBindTexture(target, texid));

	glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	//glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	//glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	glTexParameteri(target, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(target, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
	glTexParameteri(target, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);
	//glTexParameteri(target, GL_GENERATE_MIPMAP, GL_TRUE);

	glVerify(glPixelStorei(GL_UNPACK_ALIGNMENT, 1));	
	
	return texid;
}

GLuint CreateNoiseTexture(int width, int height, int depth)
{
	GLuint texid;
	glGenTextures(1, &texid);
	GLenum target = GL_TEXTURE_3D;
	
	glVerify(glBindTexture(target, texid));
	
	glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	
	//glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	//glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	
	glTexParameteri(target, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(target, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(target, GL_TEXTURE_WRAP_R, GL_REPEAT);
	//glTexParameteri(target, GL_GENERATE_MIPMAP, GL_TRUE);
	
	glVerify(glPixelStorei(GL_UNPACK_ALIGNMENT, 1));	

	// update texture
	uint8_t *data = new uint8_t[width*height*depth];
	uint8_t *ptr = data;
	
	for (int z=0; z < depth; ++z)
	{
		for (int y=0; y < height; ++y)
		{
			for(int x=0; x < width; x++)
			{					
				//*(ptr++) = fabsf(Clamp(g.GetUnsafe(x, y, z)*scale, 0.0f, 1.0f))*255.0f;
				*(ptr++) = Randf()*255;
			}
		}
	}

	glTexImage3D(target, 0, GL_LUMINANCE, width, height, depth, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, data);
	delete[] data;		
	
	return texid;
}

void UpdateVolumeTexture(GLuint texid, const Grid3D& g, float scale)
{
	// update texture
	uint8_t *data = new uint8_t[kWidth*kHeight*kDepth];
	uint8_t *ptr = data;
	
	{ ScopedTimer timer("TextureUpdate");
		
		for (int z=0; z < kDepth; ++z)
		{
			for (int y=0; y < kHeight; ++y)
			{
				for(int x=0; x < kWidth; x++)
				{					
					*(ptr++) = fabsf(Clamp(g.GetUnsafe(x, y, z)*scale, 0.0f, 1.0f))*255.0f;
					//*(ptr++) = (Length(Vec3(x, y, z)-Vec3(64,64,64)) < 32.0f)?128:0;
				}
			}
		}
	}
	
	{ ScopedTimer timer("TextureUpload");
		
		GLenum target = GL_TEXTURE_3D;
		glBindTexture(target, texid);
		
		glTexImage3D(target, 0, GL_LUMINANCE, kWidth, kHeight, kDepth, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, data);
		
		//cout << "Fluid upload took: " << (GetSeconds()-start)*1000.0f << "ms" << endl;
		
		delete[] data;		
	}	
}

void UpdateSmokeTexture(GLuint texid, const Grid3D& a, Grid3D& b)
{
	// update texture
	uint8_t *data = new uint8_t[kWidth*kHeight*kDepth];
	uint8_t *ptr = data;
	
	{ ScopedTimer timer("TextureUpdate");
		
		for (int z=0; z < kDepth; ++z)
		{
			for (int y=0; y < kHeight; ++y)
			{
				for(int x=0; x < kWidth; x++)
				{					
					*(ptr++) = fabsf(Clamp(a.GetUnsafe(x, y, z)+b.GetUnsafe(x, y, z), 0.0f, 1.0f))*255.0f;
				}
			}
		}
	}
	
	{ ScopedTimer timer("TextureUpload");
		
		GLenum target = GL_TEXTURE_3D;
		glBindTexture(target, texid);
		
		glTexImage3D(target, 0, GL_LUMINANCE, kWidth, kHeight, kDepth, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, data);
		
		//cout << "Fluid upload took: " << (GetSeconds()-start)*1000.0f << "ms" << endl;
		
		delete[] data;		
	}	
}

void DrawVolumeTexture(GLuint texid)
{
	Vector3 lightPos = g_lightRadius * SphericalToXYZ(g_lightTheta, g_lightPhi);

	Vector3 kAbsorption = Vector3(1.0f, 1.0f, 1.0f)*g_absorption;
	Vector3 kScatter = Vector3(1.0f, 1.0f, 0.90f)*g_scatter;
	Vector3 kLightIntensity = Vector3(1.0f, 1.0f, 1.0f)*g_lightIntensity;
	
	glDisable(GL_BLEND);
	glDisable(GL_TEXTURE_3D);
	glDisable(GL_TEXTURE_2D);
	
	glPushMatrix();
	glTranslatef(lightPos.x, lightPos.y, lightPos.z);
	glColor3f(1.0f, 1.0f, 1.0f);
	glutSolidSphere(0.05f, 10, 10);
	glPopMatrix();

	glEnable(GL_TEXTURE_3D);
	glEnable(GL_TEXTURE_2D);
	
	glVerify(glActiveTexture(GL_TEXTURE0));
	glVerify(glBindTexture(GL_TEXTURE_3D, g_densityTexture));

	glVerify(glActiveTexture(GL_TEXTURE1));
	glVerify(glBindTexture(GL_TEXTURE_3D, g_temperatureTexture));
	
	glVerify(glActiveTexture(GL_TEXTURE2));
	glVerify(glBindTexture(GL_TEXTURE_2D, g_blackBodyTexture));
	
	glVerify(glActiveTexture(GL_TEXTURE3));
	glVerify(glBindTexture(GL_TEXTURE_3D, g_noiseTexture));

	glVerify(glUseProgram(g_shadowShader));	
	{
		GLuint paramDensityTex = glGetUniformLocation(g_shadowShader, "g_densityTexture");
		glVerify(glUniform1i(paramDensityTex, (GLint)0)); // texunit 0
		
		GLuint paramLightPos = glGetUniformLocation(g_shadowShader, "g_lightPos");
		glVerify(glUniform3fv(paramLightPos, 1, lightPos));
		
		GLuint paramLightIntensity = glGetUniformLocation(g_shadowShader, "g_lightIntensity");
		glVerify(glUniform3fv(paramLightIntensity, 1, kLightIntensity));
		
		GLuint paramAbsorption = glGetUniformLocation(g_shadowShader, "g_absorption");
		glVerify(glUniform3fv(paramAbsorption, 1, kAbsorption));
		
		GLuint paramScatter = glGetUniformLocation(g_shadowShader, "g_scatter");
		glVerify(glUniform3fv(paramScatter, 1, kScatter));
	}
	
	const float kQuadSize = 20.0f;
	glDisable(GL_CULL_FACE);
	glBegin(GL_QUADS);
	glVertex3f(-kQuadSize, -1.0f, -kQuadSize);
	glVertex3f(kQuadSize, -1.0f, -kQuadSize);
	glVertex3f(kQuadSize, -1.0f, kQuadSize);
	glVertex3f(-kQuadSize, -1.0f, kQuadSize);
	glEnd();
		
	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	
	glVerify(glUseProgram(g_volumeShader));
	{
		GLuint paramTex = glGetUniformLocation(g_volumeShader, "g_densityTexture");
		glVerify(glUniform1i(paramTex, (GLint)0)); // texunit 0
		
		GLuint paramTemperatureTex = glGetUniformLocation(g_volumeShader, "g_temperatureTexture");
		glVerify(glUniform1i(paramTemperatureTex, (GLint)1)); // texunit 1
		
		GLuint paramBlackBodyTex = glGetUniformLocation(g_volumeShader, "g_blackBodyTexture");
		glVerify(glUniform1i(paramBlackBodyTex, (GLint)2)); // texunit 2
		
		GLuint paramNoiseTex = glGetUniformLocation(g_volumeShader, "g_noiseTexture");
		glVerify(glUniform1i(paramNoiseTex, (GLint)3)); // texuint 3
		
		GLuint paramLightPos = glGetUniformLocation(g_volumeShader, "g_lightPos");
		glVerify(glUniform3fv(paramLightPos, 1, lightPos));
		
		GLuint paramLightIntensity = glGetUniformLocation(g_volumeShader, "g_lightIntensity");
		glVerify(glUniform3fv(paramLightIntensity, 1, kLightIntensity));
		
		GLuint paramAbsorption = glGetUniformLocation(g_volumeShader, "g_absorption");
		glVerify(glUniform3fv(paramAbsorption, 1, kAbsorption));
		
		GLuint paramScatter = glGetUniformLocation(g_volumeShader, "g_scatter");
		glVerify(glUniform3fv(paramScatter, 1, kScatter));
	}
	
	glPushMatrix();
	glTranslatef(0.0f, 1.0f, 0.0f);
	glScalef(1.0f, 2.0f, 1.0f);
	glutSolidCube(2.0f);
	glPopMatrix();

	glUseProgram(0);
	
}

void DrawParticles(const std::vector<Particle>& particles)
{
	glUseProgram(0);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
	glDisable(GL_CULL_FACE);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_3D);
	glDisable(GL_TEXTURE_2D);

	glVerify(glActiveTexture(GL_TEXTURE0));
	glVerify(glBindTexture(GL_TEXTURE_3D, 0));

	glVerify(glActiveTexture(GL_TEXTURE1));
	glVerify(glBindTexture(GL_TEXTURE_3D, 0));
	
	glVerify(glActiveTexture(GL_TEXTURE2));
	glVerify(glBindTexture(GL_TEXTURE_2D, 0));
	
	glVerify(glActiveTexture(GL_TEXTURE3));
	glVerify(glBindTexture(GL_TEXTURE_3D, 0));


	glBegin(GL_POINTS);
	glColor3f(0.0f, 1.0f, 0.0f);

	for (size_t i=0; i < particles.size(); ++i)
	{
		glVertex3fv(Vec3(2.0f/kWidth, 2.0f/kHeight, 2.0f/kDepth)*particles[i].p - Vec3(1.0f, 1.0f, 1.0f));
	}

	glEnd();

}

void InitShaders()
{
	g_volumeShader = CompileProgramFromFile("RayMarchVolumeVS.glsl", "RayMarchVolumePS.glsl");
	g_shadowShader = CompileProgramFromFile("RayMarchVolumeVS.glsl", "RayMarchShadowPS.glsl");	
}

void Init()
{
	RandInit();
	
	// create volume textures
	g_densityTexture = CreateVolumeTexture(kWidth, kHeight, kDepth);	
	g_noiseTexture = CreateNoiseTexture(32, 32, 32);
	g_temperatureTexture = CreateVolumeTexture(kWidth, kHeight, kDepth);
	g_blackBodyTexture = CreateBlackBodyTexture(g_minT, g_maxT);

	InitShaders();
}

bool g_doCapture = false;

void GLUTUpdate()
{	
	static uint32_t frameCounter = 0;
	static double totalSimTime = 0.0;
	
	double startTime = GetSeconds();

	if (frameCounter == 0)
	{
		FluidStep(0.0f);
	}
	else
	{
		const uint32_t kNumSubsteps=1;
		
		for (int i=0; i < kNumSubsteps; ++i)
		{
			// for the first step just initialize
			FluidStep(1.0f/60.0f/kNumSubsteps);
		}
	}
	
	double endTime = GetSeconds();
	totalSimTime += (endTime-startTime);
	
	//UpdateSmokeTexture(g_densityTexture, *g_lastF, *g_lastS);
	UpdateVolumeTexture(g_densityTexture, *g_lastS);
	UpdateVolumeTexture(g_temperatureTexture, *g_lastT, 1.0f/4000.0f);

	cout << "Frame: " << frameCounter << " took: " << (endTime-startTime)*1000.0f << " total time: " << totalSimTime*1000.0f << endl;
	cout << "Num Particles: " << g_particles.size() << endl;
	++frameCounter;
	
	g_angle += g_angularVelocity*1.0f/60.0f;
	g_angularVelocity -= g_angularVelocity*(2.5f/60.0f);
	
	glViewport(0, 0, g_screenWidth, g_screenHeight);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0f, g_screenWidth/float(g_screenHeight), 0.01f, 1000.0f);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	gluLookAt(sin(g_angle)*6.0f, 1.2f, cos(g_angle)*6.0f, 0.0f, 1.2f, 0.0f, 0.0f, 1.0f, 0.0f);

	DrawVolumeTexture(g_densityTexture);
	//DrawVolumeTexture(g_noiseTexture);
	
	if (g_drawParticles)
		DrawParticles(g_particles);

	// flip
	glutSwapBuffers();
		
	if (!g_pause)
	{
		ScopedTimer timer("DumpImage");
		
		static int i=0;
		char buffer[255];
		sprintf(buffer, "dump/frame%d.tga", ++i);
		
		TgaImage img;
		img.m_width = g_screenWidth;
		img.m_height = g_screenHeight;
		img.m_data = new uint32_t[g_screenWidth*g_screenHeight];
		
		glReadPixels(0, 0, g_screenWidth, g_screenHeight, GL_RGBA, GL_UNSIGNED_BYTE, img.m_data);
		
		TgaSave(buffer, img);
		
		delete[] img.m_data;
	}
	
	cout << "----------------------------------------------" << endl;
}

void GLUTReshape(int width, int height)
{
}

void GLUTArrowKeys(int key, int x, int y)
{
}

void GLUTArrowKeysUp(int key, int x, int y)
{
}

bool g_spaceDown = false;

void GLUTKeyboardDown(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'e':
		{
			break;
		}
		case 's':
		{
			InitShaders();
			break;
		}
		case 't':
		{
			FluidStep(1.0f/30.0f);
			break;
		}
		case 'y':
		{
			g_pause = !g_pause;
			break;
		}
		case 'b':
		{
			break;
		}
		case ' ':
		{
			g_spaceDown = true;
			break;
		}
		case 'r':
		{
			g_xoff = Randf()*10000.f;
			//UpdateNoiseTexture(g_noiseTexture, kWidth, kHeight, kDepth, g_radius);
			break;
		}
		case 'd':
		{
			g_doCapture = true;
			break;
		}
		case '[':
		{
			g_lightRadius -= 0.1f;
			break;
		}
		case ']':
		{
			g_lightRadius += 0.1f;
			break;
		}
		case 'o':
		{
			g_lightIntensity += 0.1f;
			break;
		}
		case 'p':
		{
			g_lightIntensity -= 0.1f;
			break;
		}
		
		case 'k':
		{
			g_scatter += 0.1f;
			break;
		}
		case 'l':
		{
			g_scatter -= 0.1f;
			break;
		}
			
		case 'n':
		{
			g_absorption += 0.5f;
			break;
		}
		case 'm':
		{
			g_absorption -= 0.5f;
			break;
		}
		case 'c':
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
			break;
		}
		case 'd':
		{
			g_doCapture = false;
			break;
		}
		case ' ':
		{
			g_spaceDown = false;
			break;
		}
	}
}

static int lastx;
static int lasty;

void GLUTMouseFunc(int b, int state, int x, int y)
{	
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
	int dx = x-lastx;
	int dy = y-lasty;
	
	lastx = x;
	lasty = y;
	
	g_angularVelocity = -dx*60.0f*kPi/800.0;	
}

void GLUTPassiveMotionFunc(int x, int y)
{
	int dx = x-lastx;
	int dy = y-lasty;
	
	lastx = x;
	lasty = y;
	
	if (g_spaceDown)
	{
		g_lightTheta += DegToRad(dy);
		g_lightPhi += DegToRad(dx);
	}
}

void GenerateBlackBodyTestImage();

int main(int argc, char* argv[])
{	

	// init gl
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH);
	
	glutInitWindowSize(g_screenWidth, g_screenHeight);
	glutCreateWindow("Fluid3D");
	glutPositionWindow(350, 100);
	
#if WIN32
	glewInit();
#endif

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
	glutPassiveMotionFunc(GLUTPassiveMotionFunc);
	
	glutMainLoop();
}

