#if 1

#include <core/maths.h>
#include <core/shader.h>

#include "solve.h"

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

const int kMaxContactsPerSphere = 27; 

struct GrainSystem
{
public:
	
	Vec3* mPositions;
	Vec3* mVelocities;
	float* mRadii;
	float* mMass;

	Vec3* mCandidatePositions;
	Vec3* mNewPositions;
	float*  mNewMass;

	Vec3* mForces;

	Matrix22* mStress;

	unsigned int* mCellStarts;
	unsigned int* mCellEnds;
	unsigned int* mIndices;

	uint32_t* mSpringIndices;
	float* mSpringLengths;
	int mNumSprings;

	int* mContacts;
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


unsigned int GridHash(int x, int y, int z)
{
	const unsigned int kDim = 128;
	
	unsigned int cx = x & (kDim-1);
	unsigned int cy = y & (kDim-1);
	unsigned int cz = z & (kDim-1);
	
	return cz*kDim*kDim + cy*kDim + cx;
}

struct IndexPair
{
	IndexPair(unsigned int c, unsigned int p) : mCellId(c), mParticleId(p) {}
	
	unsigned int mCellId;
	unsigned int mParticleId;
	
	bool operator<(const IndexPair& a) const { return mCellId < a.mCellId; }
};

void ConstructGrid(float invCellEdge, int numGrains, const Vec3* positions, 
				   unsigned int* indices, unsigned int* cellStarts, unsigned int* cellEnds)
{	
	std::vector<IndexPair> indexPairs;
	indexPairs.reserve(numGrains);
	
	for (int i=0; i < numGrains; ++i)
	{
		indexPairs.push_back(IndexPair(GridHash(GridCoord(positions[i].x, invCellEdge),
												GridCoord(positions[i].y, invCellEdge),
												GridCoord(positions[i].z, invCellEdge)), i));
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

inline float sqr(float x) { return x*x; }

inline int Collide(
		Vec3 xi, float ri,
		const unsigned int* cellStarts, 
		const unsigned int* cellEnds, 
		const unsigned int* indices,
		const Vec3* positions,
		const float* radii,
		int* contacts,
		int maxContacts) 
{

	// collide particles
	int cx = GridCoord(xi.x, invCellEdge);
	int cy = GridCoord(xi.y, invCellEdge);
	int cz = GridCoord(xi.z, invCellEdge);
	
	int numContacts = 0;

	for (int i=cx-1; i <= cx+1; ++i)
	{
		for (int j=cy-1; j <= cy+1; ++j)
		{
			for (int k=cz-1; k <= cz+1; ++k)
			{
				const unsigned int cellIndex = GridHash(i, j, k);
				const unsigned int cellStart = cellStarts[cellIndex];
				const unsigned int cellEnd = cellEnds[cellIndex];
						
				// iterate over cell
				for (unsigned int i=cellStart; i < cellEnd; ++i)
				{
					const unsigned int particleIndex = indices[i];
					
					const Vec3 xj = positions[particleIndex];
					const float rj = radii[particleIndex];

					// distance to sphere
					const Vec3 xij = xi - xj; 
					
					const float dSq = LengthSq(xij);
					const float rsum = ri + rj;
				
					if (dSq < sqr(rsum) && dSq > 0.001f)
					{	
						contacts[numContacts++] = particleIndex;	

						if (numContacts == maxContacts)
							return numContacts;
					}	
				}
			}
		}
	}

	return numContacts;
}

inline Vec3 SolvePositions(
		int index,
		const Vec3* prevPositions,
		const Vec3* positions,
		const float* radii,
		const float* mass,
		float* newMass,
		const int* contacts, 
		int numContacts,
		const Vec4* planes,
		int numPlanes,
		float& pressure)

{
	Vec3 xi = positions[index];
	float ri = radii[index];

	// collide particles
	Vec3 impulse;
	float weight = 0.0f;

	float scale = 0.95f;

	ri *= scale;

	for (int i=0; i < numContacts; ++i)
	{
		const int particleIndex = contacts[i];

		const Vec3 xj = positions[particleIndex];
		const float rj = scale*radii[particleIndex];

		// distance to sphere
		const Vec3 xij = xi - xj; 
		
		const float dSq = LengthSq(xij);
		const float rsum = (ri + rj);

		//if (wsum == 0.0f)
		//	continue;
	
		if (dSq < sqr(rsum) && dSq > 0.001f)
		{
			const float d = sqrtf(dSq);
			const Vec3 n = xij / d;
	
			// project out of sphere
			impulse += 0.5f*(rsum-d)*n;
		
			//weight += wi/wsum;//1.0f; 
			weight += 1.0f;
		}
	}
	
	// collide planes
	for (int i=0; i < numPlanes; ++i)
	{
		Vec4 p = planes[i];
						
		// distance to plane
		float d = xi.x*p.x + xi.y*p.y + xi.z*p.z + p.w;
		float mtd = d - ri;
			
		if (mtd <= 0.0f)
		{
			// lerp between both to get varying friction coefficients 
			impulse += -mtd*Vec3(p.x, p.y, p.z);

			newMass[index] = 0.0f;

			// weight
			weight += 1.0f;	
		}
	}

	pressure = weight;

	return impulse/max(1.0f, weight);
}


Vec3 SolveVelocities(
		int index,
		const Vec3* positions,
		const Vec3* velocities,
		const float* radii,
		const int* contacts,
		int numContacts,
		const Vec4* planes,
		int numPlanes,
		float& weight)
{

	Vec3 impulse;
	float scale = 1.0f;

	Vec3 xi = positions[index];
	float ri = scale*radii[index];
	
	for (int i=0; i < numContacts; ++i)
	{
		const int particleIndex = contacts[i];

		const Vec3 xj = positions[particleIndex];
		const float rj = scale*radii[particleIndex];

		// vector to sphere
		const Vec3 xij = xj - xi; 
		
		const float dSq = LengthSq(xij);
		const float rsum = (ri + rj);

		if (dSq < sqr(rsum) && dSq > 0.001f)
		{
			const float d = sqrtf(dSq);
			const Vec3 n = xij / d;

			Vec3 vr = (velocities[particleIndex] - velocities[index]);

			const float kFriction = 0.0f;
			const float kViscosity = 0.0f;

			Vec3 vn = Dot(vr, n)*n; 
			Vec3 vt = vr - vn;

			// don't apply friction if separating
			if (Dot(vr, n) < 0.0f)
				impulse += kFriction*vt*0.5f;

			// apply viscosity if separating 
			if (Dot(vr, n) >= 0.0f)
				impulse += kViscosity*vn;

			//weight += 1.0f;
		}
	}

	// collide planes
	for (int i=0; i < numPlanes; ++i)
	{
		Vec4 p = planes[i];
						
		// distance to plane
		float d = xi.x*p.x + xi.y*p.y + xi.z*p.z + p.w;
		float mtd = d - ri;
	
		const float kFriction = 0.2f;
		
		if (mtd <= 0.0f)
		{
			Vec3 n(p.x, p.y, p.z);

			// friction
			Vec3 vn = Dot(velocities[index], n)*n;
			Vec3 vt = velocities[index] - vn;

			// friction
			impulse -= kFriction*vt;

			weight += 1.0f;
		}
	}

	return impulse;
}

void SolveSprings(const Vec3* positions, const uint32_t* indices, float* lengths, int numSprings, Vec3* deltas, float* weights, bool doBreak)
{
	const uint32_t* sIt = indices;
	const uint32_t* sEnd = indices+numSprings*2;
	float* lIt = lengths; 

	for (; sIt < sEnd; sIt+=2, ++lIt )
	{
		Vec3 xi = positions[sIt[0]];
		Vec3 xj = positions[sIt[1]];
		Vec3 xij = xi-xj;

		float rij = Dot(xij, xij);
		float l = *lIt;
		
		if (l < 0.0f)
			continue;

		//if (rij > sqr(l))
		{
			float k = 1.0f;	
			float d = sqrtf(rij);
		
			// handle breaking spring	
			if (d > l*1.05f)
				*lIt = -1.0f;

			float e = d-l;	
			xij / d;

			Vec3 j = 0.5f*k*e*xij;
			
			deltas[sIt[0]] -= j;
			deltas[sIt[1]] += j;

			weights[sIt[0]] += 1.0f;
			weights[sIt[1]] += 1.0f;
		}
	}	
}

void SolveSpringsGauss(Vec3* positions, const uint32_t* indices, float* lengths, int numSprings, bool dobreak)
{
	const uint32_t* sIt = indices;
	const uint32_t* sEnd = indices+numSprings*2;
	float* lIt = lengths; 

	for (; sIt < sEnd; sIt+=2, ++lIt )
	{
		Vec3 xi = positions[sIt[0]];
		Vec3 xj = positions[sIt[1]];
		Vec3 xij = xi-xj;

		float rij = Dot(xij, xij);
		float l = *lIt;
		
		if (l < 0.0f)
			continue;

		//if (rij > sqr(l))
		{
			float k = 1.0f;	
			float d = sqrtf(rij);
		
			// handle breaking spring	
			if (dobreak && d > l*1.03f)
				*lIt = -1.0f;

			float e = d-l;	
			xij /= d;

			xi -= 0.5f*k*e*xij;
			xj += 0.5f*k*e*xij;
			
			positions[sIt[0]] = xi;
			positions[sIt[1]] = xj;
		}
	}	
}

void SolveSpringDamping(const Vec3* velocities, const uint32_t* indices, float* lengths, int numSprings, Vec3* forces, float* weights)
{
	const uint32_t* sIt = indices;
	const uint32_t* sEnd = indices+numSprings*2;
	float* lIt = lengths; 

	for (; sIt < sEnd; sIt+=2, ++lIt )
	{
		float l = *lIt;
		
		if (l < 0.0f)
			continue;

		Vec3 vi = velocities[sIt[0]];
		Vec3 vj = velocities[sIt[1]];
		Vec3 vij = vi-vj;

		float k = 0.5f;	

		Vec3 j = 0.5f*k*vij;
		
		forces[sIt[0]] -= j;
		forces[sIt[1]] += j;
	
		weights[sIt[0]] += 1.0f;
		weights[sIt[1]] += 1.0f;
	}	
}

void Integrate(int index, const Vec3* positions, Vec3* candidatePositions, Vec3* velocities, const float* mass, Vec3 gravity, float damp, float dt)
{
	// v += f*dt
	velocities[index] += (gravity - damp*velocities[index])*dt;

	// x += v*dt
	candidatePositions[index] = positions[index] + velocities[index]*dt;
}

void Update(GrainSystem s, float dt, float invdt)
{		
	#pragma omp parallel for
	for (int i=0; i < s.mNumGrains; ++i)
		Integrate(i, s.mPositions, s.mCandidatePositions, s.mVelocities, s.mMass, s.mParams.mGravity, s.mParams.mDamp, dt);
	
	memset(s.mCellStarts, 0, sizeof(unsigned int)*128*128*128);
	memset(s.mCellEnds, 0, sizeof(unsigned int)*128*128*128);
	
	ConstructGrid(invCellEdge, s.mNumGrains, s.mCandidatePositions, s.mIndices, s.mCellStarts, s.mCellEnds); 

	// find neighbours
	for (int i=0; i < s.mNumGrains; ++i)
	{
		Vec3 xi = s.mCandidatePositions[i];

		s.mContactCounts[i] = Collide(s.mCandidatePositions[i],
									  s.mRadii[i],
									  s.mCellStarts,
									  s.mCellEnds,
									  s.mIndices,
									  s.mCandidatePositions,
									  s.mRadii,
									  &s.mContacts[i*kMaxContactsPerSphere],
									  kMaxContactsPerSphere);

	}

	int kNumPosIters = 3;

	for (int k=0; k < kNumPosIters; ++k)
	{
		// zero forces
		for (int i=0; i < s.mNumGrains; ++i)
		{
			s.mForces[i] = 0.0f;
			s.mNewMass[i] = 0.0f;
		}

		// solve springs 		
		//SolveSprings(s.mCandidatePositions, s.mSpringIndices, s.mSpringLengths, s.mNumSprings, s.mForces, s.mNewMass, k==kNumPosIters-1); 
		SolveSpringsGauss(s.mCandidatePositions, s.mSpringIndices, s.mSpringLengths, s.mNumSprings, k==kNumPosIters-1);
	
		// solve position constraints
		for (int i=0; i < s.mNumGrains; ++i)
		{
			float p = 0.0f;
			s.mForces[i] += SolvePositions(
					i,
					s.mPositions,
				   	s.mCandidatePositions,
				   	s.mRadii,
					s.mMass,
					s.mNewMass,
					&s.mContacts[i*kMaxContactsPerSphere],
					s.mContactCounts[i],
				   	s.mParams.mPlanes,
				   	s.mParams.mNumPlanes,
					p);

			s.mNewMass[i] += p;
			s.mMass[i] = p;
		}

		for (int i=0; i < s.mNumGrains; ++i)
			s.mCandidatePositions[i] += s.mForces[i];///max(1.0f, s.mNewMass[i]);
	}
	 

	// update velocities 
	//
	for (int i=0; i < s.mNumGrains; ++i)
		s.mVelocities[i] = (s.mCandidatePositions[i]-s.mPositions[i])*invdt; 
	
	for (int k=0; k < 1; ++k)
	{	
		for (int i=0; i < s.mNumGrains; ++i)
		{
			s.mForces[i] = 0.0f;
			s.mNewMass[i] = 0.0f;
		}
		
		SolveSpringDamping(s.mVelocities, s.mSpringIndices, s.mSpringLengths, s.mNumSprings, s.mForces, s.mNewMass);

		// solve velocity constraints
		//
		for (int i=0; i < s.mNumGrains; ++i)
		{
			s.mForces[i] += SolveVelocities(
					i,
					s.mCandidatePositions,
					s.mVelocities,
					s.mRadii,
					&s.mContacts[i*kMaxContactsPerSphere],
					s.mContactCounts[i],
					s.mParams.mPlanes,
					s.mParams.mNumPlanes,
					s.mNewMass[i]);
		}	

		//  apply velocities to positions
		//
		for (int i=0; i < s.mNumGrains; ++i)
		{
			s.mCandidatePositions[i] += s.mForces[i]*dt/max(1.0f, s.mNewMass[i]);
			//s.mVelocities[i] += s.mForces[i];///max(1.0f, s.mNewMass[i]);// (s.mNewPositions[i]-s.mPositions[i])*invdt;
		} 

		for (int i=0; i < s.mNumGrains; ++i)
			s.mVelocities[i] = (s.mCandidatePositions[i]-s.mPositions[i])*invdt;
	}

	for (int i=0; i < s.mNumGrains; ++i)
	{
		s.mVelocities[i] /= max(1.0f, s.mMass[i]*s.mParams.mDissipation); 
		//s.mVelocities[i] /= max(1.0f, s.mContactCounts[i]*0.3f); 

		s.mMass[i] = 1.0f;//
		s.mPositions[i] = s.mCandidatePositions[i];
	}
}

//------------------------------------------------------------------


GrainSystem* grainCreateSystem(int numGrains)
{
	GrainSystem* s = new GrainSystem();
	
	s->mNumGrains = numGrains;
	
	s->mPositions = (Vec3*)malloc(numGrains*sizeof(Vec3));
	s->mVelocities = (Vec3*)malloc(numGrains*sizeof(Vec3));
	s->mRadii = (float*)malloc(numGrains*sizeof(float));
	s->mMass = (float*)malloc(numGrains*sizeof(float));
	
	s->mForces = (Vec3*)malloc(numGrains*sizeof(Vec3));

	s->mContacts = (int*)malloc(numGrains*kMaxContactsPerSphere*sizeof(int));
	s->mContactCounts = (int*)malloc(numGrains*sizeof(int));

	s->mCandidatePositions = (Vec3*)malloc(numGrains*sizeof(Vec3));
	s->mNewPositions = (Vec3*)malloc(numGrains*sizeof(Vec3));

	s->mNewMass = (float*)malloc(numGrains*sizeof(float));

	for (int i=0; i < s->mNumGrains; ++i)
	{
		s->mMass[i] = 0.0f;
		s->mNewMass[i] = 0.0f;
	}

	s->mCellStarts = (unsigned int*)malloc(128*128*128*sizeof(unsigned int));
	s->mCellEnds = (unsigned int*)malloc(128*128*128*sizeof(unsigned int));
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
	
	free(s->mNewPositions);
	free(s->mCandidatePositions);

	free(s->mNewMass);

	free(s->mCellStarts);
	free(s->mCellEnds);
	free(s->mIndices);

	delete s;
}

void grainSetPositions(GrainSystem* s, float* p, int n)
{
	memcpy(&s->mPositions[0], p, sizeof(Vec3)*n);
}

void grainSetVelocities(GrainSystem* s, float* v, int n)
{
	memcpy(&s->mVelocities[0], v, sizeof(Vec3)*n);	
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
	memcpy(p, &s->mPositions[0], sizeof(Vec3)*s->mNumGrains);
}

void grainGetVelocities(GrainSystem* s, float* v)
{
	memcpy(v, &s->mVelocities[0], sizeof(Vec3)*s->mNumGrains);
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

