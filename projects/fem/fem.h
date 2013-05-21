#pragma once

#include "core/types.h"
#include "core/maths.h"
#include "core/vec2.h"
#include "core/mat22.h"

namespace fem
{

struct Scene;

struct Particle
{
	Particle() {};
	Particle(Vec2 x, float im) : p(x), invMass(im) {}

	Vec2 p;
	Vec2 v;
	Vec2 f;
    Vec2 c;
	float invMass;
	uint32_t index;
};

struct Triangle 
{
	Triangle() {};
	Triangle(int a, int b, int c) : i(a), j(b), k(c) {}

	uint32_t i, j, k;
};

struct SceneParams
{
	Vec2  mGravity;
	float mLameLambda; 
	float mLameMu;
	float mYield;
	float mCreep;
	float mDamping;
	float mDrag;
	float mFriction;
	float mToughness;

	SceneParams() :
		mGravity(0.0f, -9.8f),
		mLameLambda(1000.0f),
		mLameMu(1000.0f),
		mYield(0.5f),
		mCreep(25.0f),
		mDamping(10.0f),
		mDrag(0.0f),
		mFriction(0.5f),
		mToughness(500.0f) {}
	
};

Scene* CreateScene(
	const Particle* particles, uint32_t numParticles,
	const Triangle* elements, uint32_t numTriangles);

void DestroyScene(Scene* sim);

void SetParams(Scene* scene, SceneParams params);
SceneParams GetParams(const Scene* scene);

void SetParticles(Scene* sim, const Particle* src);
void GetParticles(const Scene* sim, Particle* dest);
uint32_t NumParticles(const Scene* sim);

void GetTriangles(const Scene* sim, Triangle* dest);
uint32_t NumTriangles(const Scene* sim);

void SetPlanes(Scene* scene, const Vec3* planes, uint32_t n);

void Update(Scene* sim, float dt);

void DrawDebug();

} // namespace fem

