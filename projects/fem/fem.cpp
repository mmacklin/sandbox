#include "fem.h"
#include "mesher.h"

#include "core/vec2.h"
#include "core/mat22.h"
#include "core/shader.h"

#include <vector>
#include <algorithm>
#include <set>

// Todo List
// --------- 
//
// * Implicit solver 
// * Amortized polar decomposition / SVD if possible 
// * Plasticity model
// * Stress averaging, temporal and spatial
// * CUDA implementation
// * Self-collision

using namespace std;

namespace fem
{

const uint32_t kInvalidIndex = uint32_t(-1);

void Print(const Matrix22& m)
{
	printf("%f, %f\n", m(0,0), m(0, 1));
	printf("%f, %f\n", m(1,0), m(1, 1));
}	

struct Element
{
	Element(const Vec2 x[3])
	{
		Vec2 e1 = x[1]-x[0]; 
		Vec2 e2 = x[2]-x[0]; 
		Vec2 e3 = x[2]-x[1];

		Matrix22 m(e1, e2);

		float det;	
		mInvDm = Inverse(m, det);

		assert(det > 0.0f);

		mB[0] = PerpCCW(e3);
		mB[1] = PerpCW(e2);
		mB[2] = PerpCCW(e1);

		printf("mInvDm:\n");
		Print(mInvDm);

		printf("mB:\n");
		printf("%f, %f\n", mB[0].x, mB[0].y);
		printf("%f, %f\n", mB[1].x, mB[1].y);
		printf("%f, %f\n", mB[2].x, mB[2].y);
	}

	Matrix22 mInvDm; // inverse rest configuration
	Matrix22 mEp;	 // plastic strain

	Vec2 mB[3];		 // area weighted normals in material space
};

struct FractureEvent
{
	uint32_t mTri;
	uint32_t mNode;
	Vec3 mPlane;
};

struct Scene
{
	vector<Particle> mParticles;

	vector<Triangle> mTriangles;
	vector<Element> mElements;
	vector<Vec3> mPlanes;

	vector<FractureEvent> mFractures;

	SceneParams mParams;
};

namespace 
{

float FrobeniusNorm(const Matrix22& m)
{
	float f = 0.0f;

	for (uint32_t i=0; i < 2; ++i)
		for (uint32_t j=0; j < 2; ++j)
			f += m(i, j)*m(i, j);

	return sqrtf(f);
}

// deformation gradient
Matrix22 CalcDeformation(const Vec2 x[3], const Matrix22& invM)
{	
	Vec2 e1 = x[1]-x[0]; 
	Vec2 e2 = x[2]-x[0]; 

	Matrix22 m(e1, e2);

	// mapping from material coordinates to world coordinates	
	Matrix22 f = m*invM; 	
	return f;
}

// calculate Green's non-linear strain tensor
Matrix22 CalcGreenStrainTensor(const Matrix22& f)
{
	Matrix22 e = 0.5f*(f*Transpose(f) - Matrix22::Identity());
	return e;
}

// calculate time derivative of Green's strain
Matrix22 CalcGreenStrainTensorDt(const Matrix22& f, const Matrix22& dfdt)
{
	Matrix22 e = 0.5f*(f*Transpose(dfdt) + dfdt*Transpose(f));
	return e;
}

// calculate Cauchy's linear strain tensor
Matrix22 CalcCauchyStrainTensor(const Matrix22& f)
{
	Matrix22 e = 0.5f*(f + Transpose(f)) - Matrix22::Identity();
	return e;
}

// calculate time derivative of Cauchy's strain tensor
Matrix22 CalcCauchyStrainTensorDt(const Matrix22& dfdt)
{
	Matrix22 e = 0.5f*(dfdt + Transpose(dfdt));
	return e;
}

// calculate isotropic Hookean stress tensor, lambda and mu are the Lame parameters
Matrix22 CalcStressTensor(const Matrix22& e, float lambda, float mu)
{
	Matrix22 s = lambda*Trace(e)*Matrix22::Identity() + mu*2.0f*e;
	return s;
}

void CollidePlanes(Particle* particles, uint32_t numParticles, const Vec3* planes, uint32_t numPlanes, float friction)
{
	for (uint32_t i=0; i < numParticles; ++i)
	{
		for (uint32_t p=0; p < numPlanes; ++p)
		{
			float d = Dot(particles[i].p, Vec2(planes[p])) + planes[p].z; 

			if (d < 0.0f)
			{
				Vec2 n(planes[p]);

				// push out of halfspace
				particles[i].p -= d*n;

				// make relative velocity separating
				float rv = Dot(particles[i].v, n);

				if (rv < 0.0f)
				{
					// zero normal velocity, material simulation will take care of restitution
					Vec2 nv = -rv*n;

					// friction
					Vec2 tv = (particles[i].v + nv)*friction;

					// update velocity
					particles[i].v = tv;
				}	
			}		
		}
	}
}

void EigenDecompose(const Matrix22& m, float& e1, float& e2)
{
	// solve the characteristic polynomial
	float a = 1.0f;
	float b = -(m(0,0) + m(1,1));
	float c = m(0,0)*m(1,1)-m(0,1)*m(1,0);

	SolveQuadratic(a, b, c, e1, e2);
}

void ShowMatrix(const Matrix22& p, const Matrix22& q, Vec2 c)
{
	float e1, e2;
	EigenDecompose(p, e1, e2);

	DrawString(0, -1, "%f %f", e1, e2);

	// calculate Eigenvectors
	Vec2 ev1 = q*Normalize(Vec2(p(0,1), e1-p(0,0)));
	Vec2 ev2 = q*Normalize(Vec2(e2-p(1,1), p(0,1)));

	//DrawString(1, 1, "%f %f %f %f", p(0,0), p(0,1), p(1,0), p(1,1));
	//DrawString(0, 1.5, "%f", Dot(p*ev2, ev2));

	glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex2fv(c);
	glVertex2fv(c+ev1);
	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex2fv(c);
	glVertex2fv(c+ev2);
	glEnd();

}

uint32_t GetVertex(const Triangle& tri, uint32_t i)
{
	assert(i < 3);
	return ((const uint32_t*)(&tri))[i];
}
uint32_t& GetVertex(Triangle& tri, uint32_t i)
{
	assert(i < 3);
	return ((uint32_t*)(&tri))[i];
}

uint32_t FindVertex(const Triangle& tri, uint32_t x)
{
	if (tri.i == x)
		return 0;
	else if (tri.j == x)
		return 1;
	else if (tri.k == x)
		return 2;
	else
		return kInvalidIndex;
}

bool FindSharedEdge(const Triangle& t1, const Triangle& t2)
{
	// count shared verts
	uint32_t sharedVerts = 0;

	for (uint32_t i=0; i < 3; ++i)
		if (FindVertex(t2, GetVertex(t1, i)) != kInvalidIndex)
			++sharedVerts;
	
	return sharedVerts > 1;	
}

struct AdjFace
{
	AdjFace(uint32_t triIndex) :
		tri(triIndex), component(-1) {}

	uint32_t tri;
	uint32_t component;
};


// detects singular vertices and separates by duplicating particles
uint32_t SeparateSingular(Particle* particles, uint32_t numParticles, Triangle* triangles, uint32_t numTriangles)
{
	const uint32_t n = numParticles;

	for (uint32_t p=0; p < n; ++p)
	{
		// find all faces connected to particle
		vector<AdjFace> adjFaces;
		for (uint32_t i=0; i < numTriangles; ++i)
		   if (FindVertex(triangles[i], p) != kInvalidIndex)
			   adjFaces.push_back(AdjFace(i));

		uint32_t componentIndex = 0;
		vector<uint32_t> stack;

		// identify connected components in the incident faces group
		for (uint32_t i=0; i < adjFaces.size(); ++i)
		{
			// dfs on each unexplored face, this is the standard
			// algorithm for identifying connected components in a graph
			if (adjFaces[i].component == kInvalidIndex)
			{
				adjFaces[i].component = componentIndex;
				stack.push_back(i);

				while (!stack.empty())
				{
					const uint32_t t = stack.back();
					stack.pop_back();

					// try and find a face connected to t
					// that is also incident on the particle
					// and doesn't belong to any other component
					for (uint32_t j=0; j < adjFaces.size(); ++j)
					{
						if (adjFaces[j].component == kInvalidIndex && FindSharedEdge(triangles[adjFaces[t].tri], triangles[adjFaces[j].tri]))
						{
							adjFaces[j].component = componentIndex;
							stack.push_back(j);
						}
					}
				}

				// have found all reachable components from i at this point
				++componentIndex;
			}
		}

		if (componentIndex > 1)
		{
			particles[p].invMass = 1.0f;

			// more than one component, create n-1 duplicates
			// of the particle and assign to each component
			for (uint32_t i=0; i < componentIndex-1; ++i)
				particles[numParticles+i] = particles[p];

			// assign the new particle to each face
			for (uint32_t i=0; i < adjFaces.size(); ++i)
			{
				assert(adjFaces[i].component != kInvalidIndex);

				// leave first component mapped to original paricle
				if (adjFaces[i].component > 0)
				{
					Triangle& tri = triangles[adjFaces[i].tri];

					uint32_t s = FindVertex(tri, p); 

					GetVertex(tri, s) = numParticles+adjFaces[i].component-1; 
				}
			}

			numParticles += componentIndex-1; 
		}
	}

	return numParticles;
}


uint32_t Fracture(Particle* particles, uint32_t numParticles, 
		Triangle* triangles, uint32_t numTriangles, 
		const FractureEvent* fractures, uint32_t numFractures)
{
	for (uint32_t f=0; f < numFractures; ++f)
	{
		const FractureEvent& event = fractures[f];

		uint32_t splitNode = GetVertex(triangles[event.mTri], event.mNode);
		
		particles[splitNode].invMass = 1.0f;

		// duplicate vertex
		particles[numParticles] = particles[splitNode];

		vector<uint32_t> left;
		vector<uint32_t> right;		

		for (uint32_t t=0; t < numTriangles; ++t)
		{
			Triangle& tri = triangles[t];

			uint32_t s = FindVertex(tri, splitNode); 

			if (s == kInvalidIndex)
				continue;

			// classify all tris to either side of the 
			// fracture plane according to triangle centroid
			Vec2 x[3] = { particles[tri.i].p, particles[tri.j].p, particles[tri.k].p };
			Vec2 c = (x[0]+x[1]+x[2])/3.0f;

			float d = Dot(c, Vec2(event.mPlane)) + event.mPlane.z;
			
			if (d < 0.0f)
				left.push_back(t);
			else
				right.push_back(t);
		}
	
		if (left.size() > 0 && right.size() > 0)
		{
			// assign left side to new particle
			for (uint32_t i=0; i < left.size(); ++i)
			{
				Triangle& tri = triangles[left[i]];

				uint32_t s = FindVertex(tri, splitNode);
				GetVertex(tri, s) = numParticles;
			}

			++numParticles;
		}	
	}	

	return numParticles;
}


uint32_t UpdateForces(Particle* particles, uint32_t numParticles, 
	const Triangle* triangles, Element* elements, uint32_t numTriangles,
   	Vec2 gravity, float lameLambda, float lameMu, float damp, float drag, float dt, 
	FractureEvent* fractures, uint32_t maxFractures, float toughness, float yield, float creep)
{

	for (uint32_t i=0; i < numParticles; ++i)
	{
		particles[i].f += particles[i].invMass>0.0f?(gravity/particles[i].invMass):Vec2(0.0f) - drag*particles[i].v;
	}

	uint32_t numFractures = 0;

	for (uint32_t i=0; i < numTriangles; ++i)
	{
		const Triangle& tri = triangles[i];
		Element& elem = elements[i];

		// read particles into a local array
		Vec2 x[3] = { particles[tri.i].p, particles[tri.j].p, particles[tri.k].p };
		Vec2 v[3] = { particles[tri.i].v, particles[tri.j].v, particles[tri.k].v };

		if (1)
		{
			Matrix22 f = CalcDeformation(x, elem.mInvDm);
			Matrix22 q = QRDecomposition(f);

			// strain 
			Matrix22 e = CalcCauchyStrainTensor(Transpose(q)*f);
			//if (FrobeniusNorm(e) > 0.2f)
			//	printf("%f\n", FrobeniusNorm(e));

			// update plastic strain
			float ef = FrobeniusNorm(e);
		
			//if (ef > yield)
			//	printf("%f\n", ef);

			if (ef > yield)
				elem.mEp += dt*creep*e;
			
			const float epmax = 0.6f;	
			if (ef > epmax)	
				elem.mEp *= epmax / ef;  

			// adjust strain
			e -= elem.mEp;

			Matrix22 s = CalcStressTensor(e, lameLambda, lameMu);

			// damping forces	
			Matrix22 dfdt = CalcDeformation(v, elem.mInvDm);
			Matrix22 dedt = CalcCauchyStrainTensorDt(Transpose(q)*dfdt);
			Matrix22 dsdt = CalcStressTensor(dedt, damp, damp);

			Matrix22 p = s + dsdt;
		
			/*	
			static int z = 0;
		   	if (1)
			{	
				Vec2 c = (x[0]+x[1]+x[2])/3.0f;
				ShowMatrix(e, q, c);
			}
			++z;	
			*/

			float e1, e2;
			EigenDecompose(p, e1, e2);

			float me = max(e1, e2);

			if (me > toughness && numFractures < maxFractures)
			{
				// calculate Eigenvector corresponding to max Eigenvalue
				Vec2 ev = q*Normalize(Vec2(p(0,1), me-p(0,0)));

				// pick a random vertex to split on
				uint32_t splitNode = rand()%3;

				// don't fracture immovable nodes
				if (particles[GetVertex(tri, splitNode)].invMass == 0.0f)
					break;

				// fracture plane perpendicular to ev
				Vec3 p(ev.x, ev.y, -Dot(ev, particles[GetVertex(tri, splitNode)].p));

				FractureEvent f = { i, splitNode, p };

				fractures[numFractures++] = f;
			}

			// calculate force on each edge due to stress and distribute to the nodes
			Vec2 f1 = q*p*elem.mB[0];
			Vec2 f2 = q*p*elem.mB[1];
			Vec2 f3 = q*p*elem.mB[2];
		
			particles[tri.i].f -= f1/3.0f;
			particles[tri.j].f -= f2/3.0f;
			particles[tri.k].f -= f3/3.0f;
		}
		else
		{
			Matrix22 f = CalcDeformation(x, elem.mInvDm);

			// elastic forces
			Matrix22 e = CalcGreenStrainTensor(f);	
			Matrix22 s = CalcStressTensor(e, lameLambda, lameMu);
	
			// damping forces	
			Matrix22 dfdt = CalcDeformation(v, elem.mInvDm);
			Matrix22 dedt = CalcGreenStrainTensorDt(f, dfdt);
			Matrix22 dsdt = CalcStressTensor(dedt, damp, damp);

			Matrix22 p = s + dsdt;

			float det;	
			Matrix22 finv = Inverse(Transpose(f), det);

			Vec2 f1 = p*(finv*elem.mB[0]);
			Vec2 f2 = p*(finv*elem.mB[1]);
			Vec2 f3 = p*(finv*elem.mB[2]);

			particles[tri.i].f -= f1/3.0f;
			particles[tri.j].f -= f2/3.0f;
			particles[tri.k].f -= f3/3.0f;
		}
	}

	return numFractures;
}

void IntegrateForces(Particle* particles, uint32_t numParticles, float dt)
{
	// integrate particles forward in time, symplectic Euler step 	
	for (uint32_t i=0; i < numParticles; ++i)
	{
		particles[i].v += particles[i].f*particles[i].invMass*dt;
		particles[i].p += particles[i].v*dt;// + 0.5f*gParticles[i].f*dt*dt*gParticles[i].invMass;

		particles[i].f = 0.0f;
	}
}

} // anonymous namespace


Scene* CreateScene(
	const Particle* particles, uint32_t numParticles,
   	const Triangle* triangles, uint32_t numTriangles)
{
	Scene* scene = new Scene();

	scene->mElements.reserve(numTriangles);

	// calculate inverse of the initial configuration
	for (uint32_t i=0; i < numTriangles; ++i)
	{
		const Triangle& t = triangles[i];

		assert(t.i < numParticles);
		assert(t.j < numParticles);
		assert(t.k < numParticles);

		// read particles into a local array
		Vec2 x[3] = { particles[t.i].p, particles[t.j].p, particles[t.k].p }; 

		scene->mElements.push_back(Element(x));	
	}


	scene->mParticles.assign(particles, particles+numParticles);
	scene->mTriangles.assign(triangles, triangles+numTriangles);
		
	// space for fractures
	scene->mFractures.resize(scene->mParticles.size());

	return scene;
}

void DestroyScene(Scene* scene)
{
	delete scene;
}

void SetParams(Scene* scene, SceneParams params)
{
	scene->mParams = params;
}

void SetPlanes(Scene* scene, const Vec3* planes, uint32_t n)
{
	scene->mPlanes.assign(planes, planes+n);
}

void SetParticles(Scene* scene, const Particle* src)
{
	scene->mParticles.assign(src, src+NumParticles(scene)); 
}

void GetParticles(const Scene* scene, Particle* dest)
{
	if (!scene->mParticles.empty())
		memcpy(dest, &scene->mParticles[0], sizeof(Particle)*NumParticles(scene));
}

uint32_t NumParticles(const Scene* scene)
{
	return uint32_t(scene->mParticles.size());
}

void GetTriangles(const Scene* scene, Triangle* dest)
{
	if (!scene->mTriangles.empty())
		memcpy(dest, &scene->mTriangles[0], sizeof(Triangle)*NumTriangles(scene));
}

uint32_t NumTriangles(const Scene* scene)
{
	return uint32_t(scene->mTriangles.size());
}

void Update(Scene* scene, float dt)
{
	Particle* particles = &scene->mParticles[0];
	Triangle* triangles = &scene->mTriangles[0];
	Element* elements = &scene->mElements[0];
	const Vec3* planes = &scene->mPlanes[0];

	uint32_t numParticles = NumParticles(scene);
	uint32_t numTriangles = NumTriangles(scene);
	uint32_t numPlanes = scene->mPlanes.size();

	SceneParams params = scene->mParams;

	FractureEvent* fractures = &scene->mFractures[0];
	uint32_t maxFractures = scene->mFractures.size();

	uint32_t numFractures = UpdateForces(particles, numParticles, triangles, elements, numTriangles,
		params.mGravity, params.mLameLambda, params.mLameMu, params.mDamping, params.mDrag, dt, 
		fractures, maxFractures, params.mToughness, params.mYield, params.mCreep);

	CollidePlanes(particles, numParticles, planes, numPlanes, params.mFriction);

	IntegrateForces(particles, numParticles, dt);

	//if (numFractures)
	//	cout << "numFractures: " << numFractures << endl;

	// todo. ugh figure out a better way to manage this
	scene->mParticles.resize(2*numParticles);
	particles = &scene->mParticles[0];

	/*
	glColor3f(1.0f, 0.0f, 0.0f);
	glBegin(GL_LINES);
	for (size_t i=0; i < numFractures; ++i)
	{
		Vec2 p = Vec2(fractures[i].mPlane)*-fractures[i].mPlane.z;
		Vec2 n = PerpCCW(Vec2(fractures[i].mPlane));

		glVertex2fv(p+n*100.0f);
		glVertex2fv(p-n*100.0f);

		glVertex2fv(p);
		glVertex2fv(p + Vec2(fractures[i].mPlane));

	}
	glEnd();
	*/

	if (params.mToughness > 0.0f)
	{
		numParticles = Fracture(particles, numParticles, triangles, numTriangles, fractures, numFractures); 

		numParticles = SeparateSingular(particles, numParticles, triangles, numTriangles);
	}

	scene->mParticles.resize(numParticles);
	particles = &scene->mParticles[0];	
}

void DrawDebug();

/*
float CalcTotalEnergy()
{
	// gravitational potential
	float gp = 0.0f;
	float ep = 0.0f;
	float ke = 0.0f;

	for (size_t i=0; i < gParticles.size(); ++i)
	{
		const Particle& p = gParticles[i];

		if (p.invMass > 0.0f)
		{
			float m = 1.0f/p.invMass;

			gp += (20.0f + p.p.y)*(-gGravity.y)*m;
			ke += 0.5f*m*Dot(p.v,p.v);
		}
	}

	for (size_t i=0; i < gTris.size(); ++i)
	{
		const Triangle& t = gTris[i];

		Matrix22 f = CalcDeformation(t);
		Matrix22 e = CalcStrainTensor(f);
		Matrix22 s = CalcStressTensor(e, gLameLambda, gLameMu);

		float a = 0.5f*Cross(gParticles[t.j].p-gParticles[t.i].p, gParticles[t.k].p-gParticles[t.i].p);

		ep += a*0.5f*(Dot(e.cols[0], s.cols[0]) + Dot(e.cols[1], s.cols[1]));
	}

	printf("%f + %f + %f = %f\n", gp, ep, ke, gp+ep+ke); 
	return 0.0f;
}
*/

} // namespace fem
