#pragma once

#include "core/maths.h"

struct GrainSystem;

struct GrainParams
{
	Vec3 mGravity;
	float mDamp;
	
	float mBaumgarte;
	float mFriction;
	float mRestitution;
	float mOverlap;
	float mDissipation;

	Vec4 mPlanes[8];
	int mNumPlanes;
};

struct GrainTimers
{
	GrainTimers() 
		: mCreateCellIndices(0.0f)
		, mSortCellIndices(0.0f)
		, mCreateGrid(0.0f)
		, mCollide(0.0f)
		, mIntegrate(0.0f)
		, mReorder(0.0f)
	{
	}

	float mCreateCellIndices;
	float mSortCellIndices;
	float mCreateGrid;
	float mCollide;
	float mIntegrate;
	float mReorder;
};

GrainSystem* grainCreateSystem(int numGrains);
void grainDestroySystem(GrainSystem* s);

void grainSetSprings(GrainSystem* s, const uint32_t* springIndices, const float* springLengths, uint32_t numSprings);

void grainSetPositions(GrainSystem* s, float* p, int n);
void grainSetVelocities(GrainSystem* s, float* v, int n);

void grainGetPositions(GrainSystem* s, float* p);
void grainGetVelocities(GrainSystem* s, float* v);

void grainSetRadii(GrainSystem* s, float* r);
void grainGetRadii(GrainSystem* s, float* r);

void grainSetParams(GrainSystem* s, GrainParams* params);

void grainUpdateSystem(GrainSystem* s, float dt, int iterations, GrainTimers* timers);
