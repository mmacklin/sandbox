#if 0

#include <core/maths.h>
#include <core/shader.h>

typedef Vec2 float2;
typedef Vec3 float3;

#include "solve.h"

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

const int kMaxContactsPerSphere = 9; 

struct Contact
{
	Vec2 point;
	Vec2 normal;
	int j;	
	float e;
	float lambda;
	float total;
};

struct GrainSystem
{
public:
	
	float2* mPositions;
	float2* mVelocities;
	float* mRadii;
	float* mMass;

	float2* mCandidatePositions;
	float2* mCandidateVelocities;
	float*  mNewMass;

	float2* mForces;

	Matrix22* mStress;

	unsigned int* mCellStarts;
	unsigned int* mCellEnds;
	unsigned int* mIndices;

	uint32_t* mSpringIndices;
	float* mSpringLengths;
	int mNumSprings;

	Contact* mContacts;
	int* mContactCounts;	

	float mMaxRadius;
	
	int mNumGrains;
	GrainParams mParams;
};

float invCellEdge = 1.0f/0.1f;

// transform a world space coordinate into cell coordinate
unsigned int GridCoord(float x, float invCellEdge)
{
	// offset to handle negative numbers
	float l = x+1000.0f;
	
	int c = (unsigned int)(floorf(l*invCellEdge));
	
	return c;
}


unsigned int GridHash(int x, int y)
{
	const unsigned int kDim = 128;
	
	unsigned int cx = x & (kDim-1);
	unsigned int cy = y & (kDim-1);
	
	return cy*kDim + cx;
}

struct IndexPair
{
	IndexPair(unsigned int c, unsigned int p) : mCellId(c), mParticleId(p) {}
	
	unsigned int mCellId;
	unsigned int mParticleId;
	
	bool operator<(const IndexPair& a) const { return mCellId < a.mCellId; }
};

void ConstructGrid(float invCellEdge, int numGrains, const float2* positions, 
				   unsigned int* indices, unsigned int* cellStarts, unsigned int* cellEnds)
{	
	std::vector<IndexPair> indexPairs;
	indexPairs.reserve(numGrains);
	
	for (int i=0; i < numGrains; ++i)
	{
		indexPairs.push_back(IndexPair(GridHash(GridCoord(positions[i].x, invCellEdge),
												GridCoord(positions[i].y, invCellEdge)), i));
	}
	
	// sort the indices based on the cell id
	std::sort(indexPairs.begin(), indexPairs.end());
	
	// scan the particle-cell array to find the start and end
	for (int i=0; i < numGrains; ++i)
	{
		IndexPair c = indexPairs[i];
		
		if (i == 0)
		{
			cellStarts[c.mCellId] = i;
		}
		else
		{
			IndexPair p = indexPairs[i-1];

			if (c.mCellId != p.mCellId)
			{
				cellStarts[c.mCellId] = i;
				cellEnds[p.mCellId] = i;			
			}
		}
		
		if (i == numGrains-1)
		{
			cellEnds[c.mCellId] = i+1;
		}
		
		indices[i] = c.mParticleId;
	}
}


inline float CollideSweptSpherePlane(const Vec2& a, const Vec2& b, const Vec3& p, float radius)
{
	float t = (radius - (a.x*p.x + a.y*p.y + p.z)) / Dot(b-a, Vec2(p.x, p.y)); 
	return t;
}

inline float norm(const Matrix22& m) { return sqrtf(Dot(m.cols[0], m.cols[0]) + Dot(m.cols[1], m.cols[1])); }
inline float sqr(float x) { return x*x; }
inline float kernel(float x) { return x; } //return x*sqr(max(1.0f-x*20.0f, 0.0f)); } 

inline int Collide(
		float2 xi, float ri,
		const unsigned int* cellStarts, 
		const unsigned int* cellEnds, 
		const unsigned int* indices,
		const float2* positions,
		const float* radii,
		const float3* planes,
		const int numPlanes,
		Contact* contacts,
		int maxContacts) 
{
	// collide particles
	int cx = GridCoord(xi.x, invCellEdge);
	int cy = GridCoord(xi.y, invCellEdge);
	
	int numContacts = 0;

	for (int i=cx-1; i <= cx+1; ++i)
	{
		for (int j=cy-1; j <= cy+1; ++j)
		{
			const unsigned int cellIndex = GridHash(i, j);
			const unsigned int cellStart = cellStarts[cellIndex];
			const unsigned int cellEnd = cellEnds[cellIndex];
					
			// iterate over cell
			for (unsigned int i=cellStart; i < cellEnd; ++i)
			{
				const unsigned int particleIndex = indices[i];
				
				const float2 xj = positions[particleIndex];
				const float rj = radii[particleIndex];

				// distance to sphere
				const float2 xij = xi - xj; 
				
				const float dSq = LengthSq(xij);
				const float rsum = ri + rj;
			
				if (dSq < sqr(rsum) && dSq > 0.001f)
				{	
					/*
					contacts[numContacts++] = particleIndex;	

					if (numContacts == maxContacts)
						return numContacts;
					*/

					if (numContacts == maxContacts)
					{
						//printf("Dropping contacts!\n");
						return numContacts;
					}

					contacts[numContacts].point = 0.5f*(xi+xj); 
					contacts[numContacts].normal = Normalize(xij);
					contacts[numContacts].e = 0.05f;

					if (contacts[numContacts].j != particleIndex)
						contacts[numContacts].lambda = 0.0f;

					contacts[numContacts].j = particleIndex; 

					++numContacts;
				}	
			}
		}
	}
	
	// collide planes
	for (int i=0; i < numPlanes; ++i)
	{
		float3 p = planes[i];

		// distance to plane
		float d = xi.x*p.x + xi.y*p.y - p.z;
			
		if (d <= ri)
		{
			if (numContacts == maxContacts)
			{
				//printf("Dropping contacts!\n");
				return numContacts;
			}

			float2 n(p.x, p.y);

			contacts[numContacts].point = xi - n*d; 
			contacts[numContacts].normal = n; 
			contacts[numContacts].e = 1.0f;
				
			if (contacts[numContacts].j != -1)
				contacts[numContacts].lambda = 0.0f;

			contacts[numContacts].j = -1;
			contacts[numContacts].lambda = 0.0f;
			
			//printf("adding contact %f %f at %f %f \n", xi.x, xi.y, contacts[numContacts].point.x, contacts[numContacts].point.y);

			++numContacts;
		}
	}

	return numContacts;
}

inline float2 SolvePositions(
		int index,
		const float2* positions,
		const float2* velocities,
		const float* radii,
		Contact* contacts, 
		int numContacts,
		float& pressure,
		float dt)

{
	float2 xi = positions[index];
	float2 vi = velocities[index];
	float ri = radii[index];

	// collide particles
	float2 impulse;

	float scale = 0.95f;

	ri *= scale;

	for (int i=0; i < numContacts; ++i)
	{
		Contact& c = contacts[i];
			
		Vec2 vj = (c.j != -1)?velocities[c.j]:0.0f;

		Vec2 delta = xi-c.point;

		float dn = Dot(delta, c.normal);
		
		if (dn < ri)
		{
			Vec2 vij = vj-vi;
			float vn = Dot(vij, c.normal);
		
			//Vec2 vt = vij - vn;
			//impulse -= c.e*SafeNormalize(vt, Vec2(0.0f))*min(Length(vn), Length(vt));

			// calculate velocity correction
			float l = max((ri-dn)/dt, vn*0.5f);
			Vec2 dv = l*c.normal;
			impulse += dv;

			c.lambda = l;


			pressure += 1.0f;
		}
	}

	return impulse/max(1.0f, pressure);
}

inline float2 SolveVelocities(
		int index,
		const float2* positions,
		const float2* velocities,
		const float* radii,
		const Contact* contacts, 
		int numContacts,
		float& pressure,
		float dt)

{
	float2 xi = positions[index];
	float2 vi = velocities[index];
	float ri = radii[index];

	// collide particles
	float2 impulse;

	float scale = 0.95f;

	ri *= scale;

	for (int i=0; i < numContacts; ++i)
	{
		const Contact& c = contacts[i];
			
		Vec2 vj = (c.j != -1)?velocities[c.j]:0.0f;

		Vec2 delta = xi-c.point;

		float dn = Dot(delta, c.normal);
		
		if (dn < ri)
		{
			Vec2 vij = vi-vj;
			Vec2 vn = Dot(vij, c.normal)*c.normal;
			Vec2 vt = vij - vn;
			impulse -= c.e*SafeNormalize(vt, Vec2(0.0f))*min(c.lambda, Length(vt));

			//impulse -= vn*c.e;

			pressure += c.e;
		}
	}

	return impulse/max(1.0f, pressure);
}
void Integrate(int index, const float2* positions, float2* candidatePositions, float2* velocities, const float* mass, float2 gravity, float damp, float dt)
{
	// v += f*dt
	velocities[index] += (gravity - damp*velocities[index])*dt;

	// x += v*dt
	candidatePositions[index] = positions[index] + velocities[index]*dt;
}

void Update(GrainSystem s, float dt, float invdt)
{		
	for (int i=0; i < s.mNumGrains; ++i)
		Integrate(i, s.mPositions, s.mCandidatePositions, s.mVelocities, s.mMass, s.mParams.mGravity, s.mParams.mDamp, dt);
	

	const int kNumIterations = 5;

	memset(s.mCellStarts, 0, sizeof(unsigned int)*128*128);
	memset(s.mCellEnds, 0, sizeof(unsigned int)*128*128);

	ConstructGrid(invCellEdge, s.mNumGrains, s.mCandidatePositions, s.mIndices, s.mCellStarts, s.mCellEnds); 
		
	for (int k=0; k < kNumIterations; ++k)
	{
		// find neighbours
		for (int i=0; i < s.mNumGrains; ++i)
		{
			float2 xi = s.mCandidatePositions[i];

			s.mContactCounts[i] = Collide(s.mCandidatePositions[i],
										  s.mRadii[i],
										  s.mCellStarts,
										  s.mCellEnds,
										  s.mIndices,
										  s.mCandidatePositions,
										  s.mRadii,
										  s.mParams.mPlanes,
										  s.mParams.mNumPlanes,
										  &s.mContacts[i*kMaxContactsPerSphere],
										  kMaxContactsPerSphere);

		}

		for (int i=0; i < s.mNumGrains; ++i)
		{
			s.mForces[i] = 0.0f;

			float pressure = 0.0f;

			float2 j = SolvePositions(
					i,
				   	s.mCandidatePositions,
					s.mVelocities,
				   	s.mRadii,
					&s.mContacts[i*kMaxContactsPerSphere],
					s.mContactCounts[i],
					pressure,
					dt);
			
				s.mMass[i] = pressure;
				s.mForces[i] += j;//*float(k+1)/kNumIterations;
		}

		for (int i=0; i < s.mNumGrains; ++i)
		{
			float pressure = 0.0f;

			s.mVelocities[i] += s.mForces[i];
			s.mForces[i] = 0.0f;

			float2 j = SolveVelocities(
					i,
				   	s.mCandidatePositions,
					s.mVelocities,
				   	s.mRadii,
					&s.mContacts[i*kMaxContactsPerSphere],
					s.mContactCounts[i],
					pressure,
					dt);
			
				s.mMass[i] = pressure;
				s.mForces[i] += j;//*float(k+1)/kNumIterations;
		}

		// apply impulses	
		for (int i=0; i < s.mNumGrains; ++i)
		{
			s.mVelocities[i] += s.mForces[i];
			s.mCandidatePositions[i] = s.mPositions[i] + dt*s.mVelocities[i];	
		}
	}

	// apply impulses	
	for (int i=0; i < s.mNumGrains; ++i)
	{
		//s.mVelocities[i] /= max(1.0f, s.mContactCounts[i]*0.3f); 
		s.mPositions[i] = s.mPositions[i] + s.mVelocities[i]*dt; 
	}
}

//------------------------------------------------------------------


GrainSystem* grainCreateSystem(int numGrains)
{
	GrainSystem* s = new GrainSystem();
	
	s->mNumGrains = numGrains;
	
	s->mPositions = (float2*)malloc(numGrains*sizeof(float2));
	s->mVelocities = (float2*)malloc(numGrains*sizeof(float2));
	s->mRadii = (float*)malloc(numGrains*sizeof(float));
	s->mMass = (float*)malloc(numGrains*sizeof(float));
	
	s->mForces = (float2*)malloc(numGrains*sizeof(float2));

	s->mContacts = (Contact*)malloc(numGrains*kMaxContactsPerSphere*sizeof(Contact));
	s->mContactCounts = (int*)malloc(numGrains*sizeof(int));

	s->mCandidatePositions = (float2*)malloc(numGrains*sizeof(float2));
	s->mCandidateVelocities = (float2*)malloc(numGrains*sizeof(float2));

	s->mNewMass = (float*)malloc(numGrains*sizeof(float));

	for (int i=0; i < s->mNumGrains; ++i)
	{
		s->mMass[i] = 1.0f;
		s->mNewMass[i] = 1.0f;
	}

	s->mCellStarts = (unsigned int*)malloc(128*128*sizeof(unsigned int));
	s->mCellEnds = (unsigned int*)malloc(128*128*sizeof(unsigned int));
	s->mIndices = (unsigned int*)malloc(numGrains*sizeof(unsigned int));

	s->mSpringIndices = NULL;
	s->mSpringLengths = NULL;
	s->mNumSprings = 0;

	return s;
}

void grainDestroySystem(GrainSystem* s)
{
	free(s->mPositions);
	free(s->mVelocities);
	free(s->mRadii);
	free(s->mMass);
	
	free(s->mForces);

	free(s->mContacts);
	free(s->mContactCounts);
	
	free(s->mCandidateVelocities);
	free(s->mCandidatePositions);

	free(s->mNewMass);

	free(s->mCellStarts);
	free(s->mCellEnds);
	free(s->mIndices);

	delete s;
}

void grainSetPositions(GrainSystem* s, float* p, int n)
{
	memcpy(&s->mPositions[0], p, sizeof(float2)*n);
}

void grainSetVelocities(GrainSystem* s, float* v, int n)
{
	memcpy(&s->mVelocities[0], v, sizeof(float2)*n);	
}

void grainSetRadii(GrainSystem* s, float* r)
{
	memcpy(&s->mRadii[0], r, sizeof(float)*s->mNumGrains);
}

void grainSetSprings(GrainSystem* s, const uint32_t* springIndices, const float* springLengths, uint32_t numSprings)
{
	s->mSpringIndices = (uint32_t*)malloc(numSprings*2*sizeof(uint32_t));
	s->mSpringLengths = (float*)malloc(numSprings*sizeof(float));

	memcpy(s->mSpringIndices, springIndices, numSprings*2*sizeof(uint32_t));
	memcpy(s->mSpringLengths, springLengths, numSprings*sizeof(float));
	
	s->mNumSprings = numSprings;
}

void grainGetPositions(GrainSystem* s, float* p)
{
	memcpy(p, &s->mPositions[0], sizeof(float2)*s->mNumGrains);
}

void grainGetVelocities(GrainSystem* s, float* v)
{
	memcpy(v, &s->mVelocities[0], sizeof(float2)*s->mNumGrains);
}

void grainGetRadii(GrainSystem* s, float* r)
{
	memcpy(r, &s->mRadii[0], sizeof(float)*s->mNumGrains);
}

void grainGetMass(GrainSystem* s, float* r)
{
	memcpy(r, &s->mMass[0], sizeof(float)*s->mNumGrains);
}


void grainSetParams(GrainSystem* s, GrainParams* params)
{
	//cudaMemcpy(s->mParams, params, sizeof(GrainParams), cudaMemcpyHostToDevice);
	s->mParams = *params;
}

void grainUpdateSystem(GrainSystem* s, float dt, int iterations, GrainTimers* timers)
{
	dt /= iterations;
			
	const float invdt = 1.0f / dt;
	
	for (int i=0; i < iterations; ++i)
	{
		Update(*s, dt, invdt);
	}	
}

#endif

