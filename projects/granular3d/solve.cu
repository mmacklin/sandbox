#if 0

#include "solve.h"

#include <cuda.h>
#include <cuda_runtime_api.h>

#include <iostream>

#ifdef _WIN32
typedef unsigned int uint32_t;
//typedef unsigned short uint32_t;
#endif

using namespace std;

#define PROFILE 0
#define USE_GRID 1
#define USE_BOX_PRUNING 0

#define kRadius 0.1f
#define kMaxRadius (kRadius)
#define kInvCellEdge (0.5f/kMaxRadius)

#if USE_GRID
typedef uint32_t CellId;
#else
typedef float CellId;
#endif

struct GrainSystem
{
public:
	
	Vec3* mPositions;
	Vec3* mVelocities;
	float* mRadii;
	
	Vec3* mSortedPositions;
	Vec3* mSortedVelocities;
	float* mSortedRadii;

	Vec3* mNewVelocities;

	uint32_t* mCellStarts;
	uint32_t* mCellEnds;
	CellId* mCellIds;
	uint32_t* mIndices;

	uint32_t mNumGrains;
	GrainParams mParams;
};

#if PROFILE

struct CudaTimer
{
	CudaTimer(const char* name, cudaEvent_t start, cudaEvent_t stop, float& timer) : mTimer(timer), mName(name), mStart(start), mStop(stop)
	{
		cudaEventRecord(mStart, 0);
	}
	
	~CudaTimer()
	{
		cudaEventRecord(mStop, 0);
		cudaEventSynchronize(mStop);
		
		float elapsedTime;
		cudaEventElapsedTime(&elapsedTime, mStart, mStop);
		
		mTimer += elapsedTime;

		//cout << mName << " took: " << elapsedTime << endl;
	}
	
	float& mTimer;
	cudaEvent_t mStart;
	cudaEvent_t mStop;
	const char* mName;
};

#else
struct CudaTimer
{
	CudaTimer(const char*, cudaEvent_t, cudaEvent_t, float& ) {}
};
#endif

void SortCellIndices(uint32_t* cellIds, uint32_t* particleIndices, uint32_t numGrains);
void SortCellIndices(float* cellIds, uint32_t* particleIndices, uint32_t numGrains);

__device__ inline float sqr(float x) { return x*x; }


// calculate collision impulse
__device__ inline Vec3 CollisionImpulse(Vec3 va, Vec3 vb, float ma, float mb, Vec3 n, float d, float baumgarte, float friction, float overlap)
{
	// calculate relative velocity
	Vec3 vd = vb-va;
	
	// calculate relative normal velocity
	float vn = Dot(vd, n);
	
	Vec3 j = Vec3(0.0f, 0.0f, 0.0f);
	
	//if (vn < 0.0f)
	vn = min(vn, 0.0f);

	{
		// calculate relative tangential velocity
		Vec3 vt = vd - n*vn;	
		float vtsq = Dot(vt, vt);
		float rcpvt = rsqrtf(vtsq);// + 0.001f);
		
		// position bias
		float bias = baumgarte*min(d+overlap, 0.0f);

		Vec3 jn = -(vn + bias)*n;
		Vec3 jt = max(friction*vn*rcpvt, -1.0f)*vt;
		
		// crappy static friction
		if (fabsf(vtsq*rcpvt) < fabsf(friction*vn*2.0f) && vn < 0.0f)
			jt = -vt;				

		// total mass 
		float msum = ma + mb;
	
		// normal impulse
		j = (jn + jt)*mb/msum;
	}
	
	return j;
}

#if USE_GRID

const uint32_t kGridDim = 128;

// transform a world space coordinate into cell coordinate
__device__ inline uint32_t GridCoord(float x, float invCellEdge)
{
	// offset to handle negative numbers
	float l = x+1000.0f;
	
	uint32_t c = (uint32_t)(floorf(l*invCellEdge));
	return c;
}

__device__ inline uint32_t GridHash(int x, int y, int z)
{	
	uint32_t cx = x & (kGridDim-1);
	uint32_t cy = y & (kGridDim-1);
	uint32_t cz = z & (kGridDim-1);
	
	return cy*(kGridDim*kGridDim) + cx*kGridDim + cz;
}

/*
__device__ inline uint32_t GridHash(int x, int y, int z)
{
	const uint32_t p1 = 73856093; 
	const uint32_t p2 = 19349663;
	const uint32_t p3 = 53471161;
		
	uint32_t n = x*p1 ^ y*p2 ^ z*p3;
	return n&(kGridDim*kGridDim*kGridDim-1);
}
*/

__global__ void CreateCellIndices(const Vec3* positions, uint32_t* cellIds, uint32_t* particleIndices)
{
	uint32_t i = blockIdx.x*blockDim.x + threadIdx.x;

	Vec3 p = positions[i];
	
	cellIds[i] = GridHash(GridCoord(p.x, kInvCellEdge), GridCoord(p.y, kInvCellEdge), GridCoord(p.z, kInvCellEdge));
	particleIndices[i] = i;	
}

__global__ void CreateGrid(const uint32_t* cellIds, uint32_t* cellStarts, uint32_t* cellEnds, uint32_t numGrains)
{	
	uint32_t i = blockIdx.x*blockDim.x + threadIdx.x;
	
	// scan the particle-cell array to find the start and end
	uint32_t c = cellIds[i];
	
	if (i == 0)
	{
		cellStarts[c] = i;
	}
	else
	{
		uint32_t p = cellIds[i-1];

		if (c != p)
		{
			cellStarts[c] = i;
			cellEnds[p] = i;
		}
	}
	
	if (i == numGrains-1)
	{
		cellEnds[c] = i+1;
	}
}

__device__ inline Vec3 CollideSphere(Vec3 xa, Vec3 xb, Vec3 va, Vec3 vb, float ra, float rb, float baumgarte, float friction, float overlap)
{
	// distance to sphere
	Vec3 t = xa - xb;
	Vec3 j = Vec3(0.0f, 0.0f, 0.0f);

	float d = Dot(t, t);
	float rsum = ra + rb;
	float mtd = d - sqr(rsum);
			
	if (mtd < 0.0f)
	{
		Vec3 n = Vec3(0.0f, 1.0f, 0.0f);
				
		if (d > 0.0f)
		{
			float rcpDist = rsqrtf(d);

			n = t * rcpDist;
			d = d * rcpDist;
		}
				
		j = CollisionImpulse(vb, va, 1.0f, 1.0f, n, d-rsum, baumgarte, friction, overlap);
	}

	return j;
}

__device__ inline Vec3 CollideCell(int index, int cx, int cy, int cz, const uint32_t* cellStarts, const uint32_t* cellEnds, const uint32_t* indices,
				 const Vec3* positions, const Vec3* velocities, const float* radii, Vec3 x, Vec3 v, float r, float baumgarte, float friction, float overlap)
{
	Vec3 j = Vec3(0.0f, 0.0f, 0.0f);
	
	uint32_t cellIndex = GridHash(cx, cy, cz);
	uint32_t cellStart = cellStarts[cellIndex];
	uint32_t cellEnd = cellEnds[cellIndex];
	
	for (int i=cellStart; i < cellEnd; ++i)
	{
		uint32_t particleIndex = i;//indices[i];
		
		if (particleIndex != index)
		{		
			j += CollideSphere(x, positions[particleIndex], v, velocities[particleIndex], r, radii[particleIndex], baumgarte, friction, overlap);
		}		
	}
	
	return j;
}


#endif


__global__ void ReorderParticles(const Vec3* positions, const Vec3* velocities, const float* radii, Vec3* sortedPositions, Vec3* sortedVelocities, float* sortedRadii, const uint32_t* indices)
{
	uint32_t i = blockIdx.x*blockDim.x + threadIdx.x;
	
	int originalIndex = indices[i];

	sortedPositions[i] = positions[originalIndex];
	sortedVelocities[i] = velocities[originalIndex];
	sortedRadii[i] = radii[originalIndex];
}


__global__ void Collide(const Vec3* positions, const Vec3* velocities, const float* radii, const uint32_t* cellStarts, const uint32_t* cellEnds, const uint32_t* indices,
						Vec3* newVelocities, int numGrains, GrainParams params, float dt, float scale)
{
	const int index = blockIdx.x*blockDim.x + threadIdx.x;
		
	const Vec3 x = positions[index];
	const Vec3 v = velocities[index];
	const float  r = radii[index];

	Vec3 vd = Vec3(0.0f, 0.0f, 0.0f);

#if USE_GRID

	// collide particles
	int cx = GridCoord(x.x, kInvCellEdge);
	int cy = GridCoord(x.y, kInvCellEdge);
	int cz = GridCoord(x.z, kInvCellEdge);

	for (int k=cz-1; k <= cz+1; ++k)
	{
		for (int j=cy-1; j <= cy+1; ++j)
		{
			for (int i=cx-1; i <= cx+1; ++i)
			{
				vd += CollideCell(index, i, j, k, cellStarts, cellEnds, indices, positions, velocities, radii, x, v, r, params.mBaumgarte, params.mFriction, params.mOverlap);
			}
		}
	}
#endif

	// collide planes
	for (int i=0; i < params.mNumPlanes; ++i)
	{
		Vec4 p = params.mPlanes[i];
						
		// distance to plane
		float d = x.x*p.x + x.y*p.y + x.z*p.z + p.w;
			
		float mtd = d - r;
			
		if (mtd < 0.0f)
		{
			vd += CollisionImpulse(Vec3(0.0f, 0.0f, 0.0f), v, 0.0f, 1.0f, Vec3(p.x, p.y, p.z), mtd, params.mBaumgarte, 0.8f, params.mOverlap);
		}
	}
	
	// write back velocity
	newVelocities[index] = v + vd * scale;
}

__global__ void IntegrateForce(Vec3* velocities, Vec3 gravity, float damp, float dt)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;

	velocities[index] += (gravity - damp*velocities[index])*dt;
}


__global__ void IntegrateVelocity(Vec3* positions, Vec3* velocities, const Vec3* newVelocities, float dt)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;

	// x += v*dt
	velocities[index] = newVelocities[index];
	positions[index] += velocities[index]*dt;
}

/*
__global__ void PrintCellCounts(uint32_t* cellStarts, uint32_t* cellEnds)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;

	printf("%d\n", cellEnds[index]-cellStarts[index]);

}
*/

//------------------------------------------------------------------


GrainSystem* grainCreateSystem(int numGrains)
{
	GrainSystem* s = new GrainSystem();
	
	s->mNumGrains = numGrains;
	
	cudaMalloc(&s->mPositions, numGrains*sizeof(Vec3));
	cudaMalloc(&s->mVelocities, numGrains*sizeof(Vec3));
	cudaMalloc(&s->mNewVelocities, numGrains*sizeof(Vec3));
	cudaMalloc(&s->mRadii, numGrains*sizeof(float));
	
	cudaMalloc(&s->mSortedPositions, numGrains*sizeof(Vec3));
	cudaMalloc(&s->mSortedVelocities, numGrains*sizeof(Vec3));
	cudaMalloc(&s->mSortedRadii, numGrains*sizeof(float));

	// grid
#if USE_GRID
	cudaMalloc(&s->mCellStarts, kGridDim*kGridDim*kGridDim*sizeof(uint32_t));
	cudaMalloc(&s->mCellEnds, kGridDim*kGridDim*kGridDim*sizeof(uint32_t));
#endif

	cudaMalloc(&s->mCellIds, numGrains*sizeof(uint32_t));
	cudaMalloc(&s->mIndices, numGrains*sizeof(uint32_t));
	
	return s;
}

void grainDestroySystem(GrainSystem* s)
{
	cudaFree(s->mPositions);
	cudaFree(s->mVelocities);
	cudaFree(s->mNewVelocities);
	cudaFree(s->mRadii);	
	
	cudaFree(s->mSortedPositions);
	cudaFree(s->mSortedVelocities);
	cudaFree(s->mSortedRadii);	
	
#if USE_GRID
	cudaFree(s->mCellStarts);
	cudaFree(s->mCellEnds);
#endif
	cudaFree(s->mCellIds);
	cudaFree(s->mIndices);

	delete s;
}
void grainSetSprings(GrainSystem* s, const uint32_t* springIndices, const float* springLengths, uint32_t numSprings)
{
	/*
	s->mSpringIndices = (uint32_t*)malloc(numSprings*2*sizeof(uint32_t));
	s->mSpringLengths = (float*)malloc(numSprings*sizeof(float));

	memcpy(s->mSpringIndices, springIndices, numSprings*2*sizeof(uint32_t));
	memcpy(s->mSpringLengths, springLengths, numSprings*sizeof(float));
	
	s->mNumSprings = numSprings;
	*/
}


void grainSetPositions(GrainSystem* s, float* p, int n)
{
	cudaMemcpy(&s->mPositions[0], p, sizeof(Vec3)*n, cudaMemcpyHostToDevice);
}

void grainSetVelocities(GrainSystem* s, float* v, int n)
{
	cudaMemcpy(&s->mVelocities[0], v, sizeof(Vec3)*n, cudaMemcpyHostToDevice);	
}

void grainSetRadii(GrainSystem* s, float* r)
{
	cudaMemcpy(&s->mRadii[0], r, sizeof(float)*s->mNumGrains, cudaMemcpyHostToDevice);
}

void grainGetPositions(GrainSystem* s, float* p)
{
	cudaMemcpy(p, &s->mPositions[0], sizeof(Vec3)*s->mNumGrains, cudaMemcpyDeviceToHost);
}

void grainGetVelocities(GrainSystem* s, float* v)
{
	cudaMemcpy(v, &s->mVelocities[0], sizeof(Vec3)*s->mNumGrains, cudaMemcpyDeviceToHost);
}

void grainGetRadii(GrainSystem* s, float* r)
{
	cudaMemcpy(r, &s->mRadii[0], sizeof(float)*s->mNumGrains, cudaMemcpyDeviceToHost);
}

void grainSetParams(GrainSystem* s, GrainParams* params)
{
	//cudaMemcpy(s->mParams, params, sizeof(GrainParams), cudaMemcpyHostToDevice);
	s->mParams = *params;
}

void grainUpdateSystem(GrainSystem* s, float dt, int iterations, GrainTimers* timers)
{
	//iterations = 10;

	dt /= iterations;

	const int kNumThreadsPerBlock = 128;
	const int kNumBlocks = s->mNumGrains / kNumThreadsPerBlock;

	GrainParams params = s->mParams;
	params.mBaumgarte /= dt;
	
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaFuncSetCacheConfig(CreateCellIndices, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(CreateGrid, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(ReorderParticles, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(IntegrateForce, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(IntegrateVelocity, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(Collide, cudaFuncCachePreferL1);

	for (int i=0; i < iterations; ++i)
	{
		{
			CudaTimer timer("CreateCellIndices", start, stop, timers->mCreateCellIndices);
			
			CreateCellIndices<<<kNumBlocks, kNumThreadsPerBlock>>>(s->mPositions, s->mCellIds, s->mIndices);
		}

		{ 
			CudaTimer timer("SortCellIndices", start, stop, timers->mSortCellIndices);
			
			SortCellIndices(s->mCellIds, s->mIndices, s->mNumGrains);
		}

#if USE_GRID
		{
			CudaTimer timer("CreateGrid", start, stop, timers->mCreateGrid);
			
			cudaMemset(s->mCellStarts, 0, sizeof(uint32_t)*kGridDim*kGridDim*kGridDim);
			cudaMemset(s->mCellEnds, 0, sizeof(uint32_t)*kGridDim*kGridDim*kGridDim);

			CreateGrid<<<kNumBlocks, kNumThreadsPerBlock>>>(s->mCellIds, s->mCellStarts, s->mCellEnds, s->mNumGrains);
		}
#endif

		{
			CudaTimer timer("ReorderParticles", start, stop, timers->mReorder);

			ReorderParticles<<<kNumBlocks, kNumThreadsPerBlock>>>(s->mPositions, s->mVelocities, s->mRadii, s->mSortedPositions, s->mSortedVelocities, s->mSortedRadii, s->mIndices);
		}
		
		//PrintCellCounts<<<kGridDim*kGridDim/kNumThreadsPerBlock, kNumThreadsPerBlock>>>(s->mCellStarts, s->mCellEnds);

		{
			float t;
			CudaTimer timer("Integrate Force", start, stop, t);

			IntegrateForce<<<kNumBlocks, kNumThreadsPerBlock>>>(s->mSortedVelocities, s->mParams.mGravity, s->mParams.mDamp, dt);
		}

		{
			CudaTimer timer("Collide", start, stop, timers->mCollide);
			
			float scale = 1;//float(i+1)/(iterations);

			Collide<<<kNumBlocks, kNumThreadsPerBlock>>>(s->mSortedPositions, s->mSortedVelocities, s->mSortedRadii, s->mCellStarts, s->mCellEnds, s->mIndices, s->mNewVelocities, s->mNumGrains, params, dt, scale);
		}

		{
			CudaTimer timer("Integrate", start, stop, timers->mIntegrate);
	
			IntegrateVelocity<<<kNumBlocks, kNumThreadsPerBlock>>>(s->mSortedPositions, s->mSortedVelocities, s->mNewVelocities, dt); 
		}
	
		swap(s->mSortedPositions, s->mPositions);
		swap(s->mSortedVelocities, s->mVelocities);
		swap(s->mSortedRadii, s->mRadii);
		
	}		


	cudaEventDestroy(start);
	cudaEventDestroy(stop);
}

#endif
