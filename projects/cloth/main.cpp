#include <core/types.h>
#include <core/maths.h>
#include <core/platform.h>
#include <core/shader.h>

#include <iostream>

using namespace std;

typedef float real;
typedef XVector2<real> Vec2r;
typedef XVector3<real> Vec3r;

const uint32_t kHeight = 600;
const uint32_t kWidth = 600;
const real kWorldSize = 4.0f;
const real kZoom = kWorldSize + 1.0f;

// to allow switching between floats / doubles
#define glVertex2x(x) glVertex2fv(x)


struct Bounds
{
	Bounds(const Vec2r& p) : lower(p), upper(p) {}
	
	Vec2r lower;
	Vec2r upper;
	
	void Expand(const Vec2r& p)
	{
		lower = Min(lower, p);
		upper = Max(upper, p);
	}
};

struct Particle
{
	Particle(const Vec2r& x, real invMass=1.0f)
		: mX(x), mInvMass(invMass)
	{
	}
	
	Vec2r mX;
	Vec2r mV;
	Vec2r mF;
	
	real mInvMass;
	
	void ApplyForce(const Vec2r& f)
	{
		/*
		float l = Length(f);
		if(l > 3000.0f);
		{
			glBegin(GL_POINTS);
			glVertex2x(mX);
			glEnd();
			
			cout << l << endl;
			
			glutSwapBuffers();
			//assert(0);
		}
		*/
		mF += f;
	}
	
	void ApplyImpulse(const Vec2r& i)
	{
		//assert(Length(i) < 10.0f);
		
		mV += i;
	}
};

struct Force
{
	virtual void ApplyForce()=0;
	virtual void Draw()=0;
};

struct MouseSpring : public Force
{
	MouseSpring(Particle* p, const Vec2r& m) : mP(p), mMouse(m), mRestLength(1.0f)
	{
	}
	
	virtual void ApplyForce()
	{
		const real kStiffness = 3.0;
		const real kDamping = 5.0;
		
		// position delta
		Vec2r dx = mP->mX - mMouse;
		real dl = Length(dx);
		dx /= dl;
		
		// velocity delta
		Vec2r dv = mP->mV - mMouse;
		
		Vec2r f = (kStiffness*(dl-mRestLength) + kDamping*(Dot(dv, dx)))*dx;
		
		mP->ApplyForce(-f);
	}
	
	virtual void Draw()
	{
		glBegin(GL_LINES);
		glColor3f(1.0f, 1.0f, 0.0f);
		glVertex2x(mP->mX);
		glVertex2x(mMouse);
		glEnd();
	}
	
	Particle* mP;
	Vec2r mMouse;
	
	real mRestLength;
};


struct LinearSpring : public Force
{
	LinearSpring(Particle* p, Particle* q) : mP(p), mQ(q), mRestLength(Length(p->mX-q->mX)), mVisited(false)
	{
	}

	LinearSpring(Particle* p, Particle* q, real l) : mP(p), mQ(q), mRestLength(l)
	{
	}
	
	
	virtual void ApplyForce()
	{
		const real kStiffness = 10000.0f;
		const real kDamping = 50.0f;
		
		// position delta
		Vec2r dx = mP->mX - mQ->mX;
		real dl = Length(dx);
		dx /= dl;
		
		// velocity delta
		Vec2r dv = mP->mV - mQ->mV;
		
		Vec2r f = (kStiffness*(dl-mRestLength) + kDamping*(Dot(dv, dx)))*dx;
		
		mP->ApplyForce(-f);
		mQ->ApplyForce(f);
	}
	
	virtual void Draw()
	{
		glBegin(GL_LINES);
		glColor3f(0.9f, 0.85f, 0.6f);
		glVertex2x(mP->mX);
		glVertex2x(mQ->mX);
		glEnd();
	}
	
	Bounds CalculateBounds(real dt) const
	{
		Bounds b(mP->mX);
		b.Expand(mP->mX + mP->mV*dt);
		b.Expand(mQ->mX);
		b.Expand(mQ->mX + mQ->mV*dt);
		
		return b;
	}
	
	Particle* mP;
	Particle* mQ;
	
	real mRestLength;
	bool mVisited;
};

struct AngularSpring : public Force
{
	AngularSpring(Particle* p, Particle* q, Particle* r) : mP(p), mQ(q), mR(r)
	{
	}

	virtual void ApplyForce()
	{
		const real kStiffness = 100.0;
		const real kDamping = 0.5;
	
		Vec2r a = mQ->mX - mP->mX;
		Vec2r b = mR->mX - mQ->mX;

		real lenA = Length(a);
		real lenB = Length(b);
		
		real adotb = Dot(a, b);
		
		real s = -0.5f / (lenA * lenB);
		
		Vec2r cdp = s*(-b + a*(adotb/(lenA*lenA)));
		Vec2r cdr = s*( a - b*(adotb/(lenB*lenB)));
		
		Vec2r fp = -kStiffness*cdp;
		Vec2r fr = -kStiffness*cdr;
		Vec2r fq = -(fp+fr);
		
		mP->ApplyForce(fp);
		mQ->ApplyForce(fq);
		mR->ApplyForce(fr);
	
	}
	
	virtual void Draw()
	{
		/*
		glBegin(GL_LINE_STRIP);
		glColor3f(1.0f, 1.0f, 0.0f);
		glVertex2x(mP->mX);
		glVertex2x(mQ->mX);
		glVertex2x(mR->mX);
		glEnd();
		 */
	}
	
	Particle* mP;
	Particle* mQ;
	Particle* mR;
	
};

Vec2r g_gravity(0.0f, -9.8f);
real g_drag = 2.5f;

std::vector<Particle*> g_particles;
std::vector<Force*> g_forces;
std::vector<Vec3r> g_planes;
MouseSpring* g_mouseSpring;


template <typename T>
class FixedGrid
{
public:

	typedef int Index;
	static const Index kNull = -1;
	
	struct Entry
	{
		T mValue;
		int mNext;
	};
	
	FixedGrid(float left, float right, float bottom, float top, float dx, uint32_t maxEntries) :
		mLeft(left),
		mRight(right),
		mBottom(bottom),
		mTop(top)
	{
		assert(right > left);
		assert(top > bottom);
		
		mRcpDx = 1.0f/dx;
		
		mWidth = int(ceilf((right-left) / dx));
		mHeight = int(ceilf((top-bottom) / dx));
		
		uint32_t numCells = mWidth*mHeight;
		
		mGrid.resize(numCells, kNull);
		mEntries.reserve(maxEntries);
		
	}
	
	template <bool unique>
	void Insert(const T& value, float left, float right, float bottom, float top)
	{		
		// rasterize onto the grid
		int startx, starty, endx, endy;
		GetCoord(left, bottom, startx, starty);
		GetCoord(right, top, endx, endy);
		
		startx = max(startx-1, 0);
		starty = max(starty-1, 0);
		endx = min(endx+1, mWidth-1);
		endy = min(endy+1, mHeight-1);
		
		for (int y=starty; y <= endy; ++y)
		{
			for (int x=startx; x <= endx; ++x)
			{
				Index& head = GetCell(x, y);
				
				if (unique)
				{
					// check if value already exists, if so don't re-add
					bool found = false;
					Index iter = head;
					while (iter != kNull)
					{
						if (mEntries[iter].mValue == value)
						{
							found = true;
							break;
						}
						
						iter = mEntries[iter].mNext;
					}
					
					if (found)
						break;
				}
				
				Entry e;
				e.mValue = value;
				e.mNext = head;

				head = Index(mEntries.size());
				mEntries.push_back(e);
			}
		}
	}
	

	void Query(float left, float right, float bottom, float top, std::vector<T>& out)
	{
		if (0)
		{
			for (size_t  i=0; i < g_forces.size(); ++i)
			{
				LinearSpring* s = dynamic_cast<LinearSpring*>(g_forces[i]);
				if (s)
					out.push_back(s);
			}
			return;
		}
		
		int startx, starty, endx, endy;
		GetCoord(left, bottom, startx, starty);
		GetCoord(right, top, endx, endy);
		
		startx = max(startx-1, 0);
		starty = max(starty-1, 0);
		endx = min(endx+1, mWidth-1);
		endy = min(endy+1, mHeight-1);
		
		for (int y=starty; y <= endy; ++y)
		{
			for (int x=startx; x <= endx; ++x)
			{
				Index iter = GetCell(x, y);
				
				while (iter != kNull)
				{
					const Entry& e = mEntries[iter];
					
					/*
					if (find(out.begin(), out.end(), e.mValue) == out.end())
					{
						// invoke user callback
						out.push_back(e.mValue);	
					}
					*/
					if(!e.mValue->mVisited)
					{
						e.mValue->mVisited = true;
						out.push_back(e.mValue);
					}
					
					// move forward
					iter = e.mNext;
				}
			}
		}
		
		for (size_t i=0; i < out.size(); ++i)
			out[i]->mVisited = false;
	}
	
	void GetCoord(float x, float y, int& ix, int& iy)
	{
		ix = int((x-mLeft)*mRcpDx);
		iy = int((y-mBottom)*mRcpDx);
	}	
			
	Index& GetCell(int x, int y)
	{
		return mGrid[y*mWidth + x];
	}	
		
	float mLeft;
	float mRight;
	float mBottom;
	float mTop;
	
	int mWidth;
	int mHeight;
	
	float mRcpDx;
	
	std::vector<Entry> mEntries;
	std::vector<Index> mGrid;
	
};

bool PointOnLineSegment(const Vec2r& p, const Vec2r& q, const Vec2r& r, real& s)
{
	Vec2r pq = p-q;
	Vec2r qr = r-q;
	real lSq = LengthSq(qr);
	
	s = Dot(pq, qr) / lSq;
	
	return s >= 0.0f && s <= 1.0f;
}

bool IntersectPointLineSegmentContinuous(const Vec2r& p, const Vec2r& pv, const Vec2r& q, const Vec2r& qv, const Vec2r& r, const Vec2r& rv, real& t, real& s, real maxt, bool draw=false)
{
	Vec2r rq = PerpCCW(r-q);
	Vec2r pq = p-q;
	Vec2r pvqv = pv-qv;
	Vec2r rvqv = PerpCCW(rv-qv);
	
	real c = Dot(rq, pq);
	real b = Dot(pvqv, rq) + Dot(rvqv, pq);
	real a = Dot(rvqv, pvqv);
	
	real minT, maxT;
		
	bool hit = SolveQuadratic(a, b, c, minT, maxT);
	
	if (!hit)
		return false;

	t = minT;
	
	if (minT < 0.0)
		t = maxT;
	
	if (t < 0.0 || t > maxt)
		return false;
	
	// test if point lies within line segment at t
	if (!PointOnLineSegment(p+t*pv, q+t*qv, r+t*rv, s))
		return false;
	
	if (draw)
	{
		Vec2r c = p + t*pv;
		glMatrixMode(GL_PROJECTION);
		glLoadMatrixf(OrthographicMatrix(c.x-0.01f, c.x+0.01f, c.y-0.01f, c.y+0.01f, 0.0f, 1.0f));

		glBegin(GL_QUADS);
		
		if (hit)
			glColor3f(0.0f, 1.0f, 0.0f);
		else
			glColor3f(1.0f, 0.0f, 0.0f);

		glVertex2x(q);
		glVertex2x(r);
		glVertex2x(r+rv*maxt);
		glVertex2x(q+qv*maxt);
		glEnd();
		
		glBegin(GL_LINES);
		glColor3f(1.0f, 1.0f, 1.0f);
		glVertex2x(p);
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex2x(p+pv*maxt);
		
		glColor3f(0.0f, 0.0f, 1.0f);
		glVertex2x(q+t*qv);
		glVertex2x(r+t*rv);
		
		glEnd();
		
		glBegin(GL_POINTS);
		glColor3f(1.0f, 1.0f, 0.0f);
		glVertex2x(p+t*pv);
		glColor3f(1.0f, 0.0f, 1.0f);
		glVertex2x(p+maxT*pv);
		glEnd();
		
		//glFinish();
		glutSwapBuffers();
		int i;
		cin >> i;
	}
	
	return true;
}


void AddPlane(const Vec2r& p, const Vec2r& n)
{
	g_planes.push_back(Vec3r(n.x, n.y, -Dot(p, n)));
}
					   
void DrawPlanes()
{
	for (size_t i=0; i < g_planes.size(); ++i)
	{
		Vec2r p = -g_planes[i].z * Vec2r(g_planes[i].x, g_planes[i].y);
		Vec2r d = Vec2r(-g_planes[i].y, g_planes[i].x);
			
		glBegin(GL_LINES);
		glColor3f(1.0f, 1.0f, 1.0f);
		glVertex2x(p - d*1000.0);
		glVertex2x(p + d*1000.0);
		glEnd();
	}
}

void AddSpaghetti(float dt)
{
	static Particle* p0=NULL;
	static Particle* p1=NULL;
	static Particle* p2=NULL;
		
	static int length = 0;
	static double rate = 0.05f;
	static double t = 0.0f;
	static double lastT = 0.0f;
	
	t += dt;
	
	if (t-lastT > rate)
	{	
		// shuffle down
		p0 = p1;
		p1 = p2;
		
		const Vec2r initialX(3.0f, 4.0f);
		const Vec2r initialV(-10.0f, -2.0f);
		
		p2 = new Particle(initialX);
		p2->mV = Vec2r(-10.0f, -2.0f);

		g_particles.push_back(p2);
		
		if (p1)
		{
			g_forces.push_back(new LinearSpring(p2, p1));//, 0.2f));
			
			if (p0)	
			{
				g_forces.push_back(new AngularSpring(p0, p1, p2));
			}
		}
		
		++length;
		lastT = t;
	}
	
	if (length == 13)
	{
		length = 0;
		p0 = p1 = p2 = NULL;
		lastT += 2.0f;
	}
}

void AddLine(const Vec2r& p, const Vec2r& pv, const Vec2r& q, const Vec2r& qv)
{
	Particle* a = new Particle(p);
	a->mV = pv;
	
	Particle* b = new Particle(q);
	b->mV = qv;
	
	g_forces.push_back(new LinearSpring(a, b));
	g_particles.push_back(a);
	g_particles.push_back(b);
}

void AddTri(const Vec2r& centroid, real edgeLength)
{
	real halfHeight = edgeLength*Sin(DegToRad(60.0f))*0.5f;
	real halfWidth = 0.5f*edgeLength;
	
	int i = g_particles.size();
	
	g_particles.push_back(new Particle(Vec2r(0.0f, halfHeight), 1.0f));
	g_particles.push_back(new Particle(Vec2r(-halfWidth, -halfHeight), 1.0f));
	g_particles.push_back(new Particle(Vec2r(halfWidth, -halfHeight), 1.0f));
	
	g_forces.push_back(new LinearSpring(g_particles[i], g_particles[i+1]));
	g_forces.push_back(new LinearSpring(g_particles[i+1], g_particles[i+2]));
	g_forces.push_back(new LinearSpring(g_particles[i+2], g_particles[i]));

	g_particles[i]->mV += Vec2r(8.0f, 15.0f);

}


void Init()
{
	real x = -Sin(DegToRad(-0.0f));
	real y = Cos(DegToRad(-0.0f));

	AddPlane(Vec2r(0.0f, -2.0f), Vec2r(x, y));
	AddPlane(Vec2r(4.0f, 0.0f), Vec2r(-1.0f, 0.0f));
	AddPlane(Vec2r(-4.0f, 0.0f), Vec2r(1.0f, 0.0f));
	
}

void Reset()
{
	g_particles.resize(0);
	g_forces.resize(0);
	
	Init();
}

void ApplyForces()
{
	// apply gravity and damping
	for (size_t  i=0; i < g_particles.size(); ++i)
	{
		Particle& p = *g_particles[i];

		p.mF = Vec2r(0.0f, 0.0f);
		p.ApplyForce(-g_drag*p.mV);

		if (p.mInvMass > 0.0f)
			p.ApplyForce(g_gravity/p.mInvMass);
	}
	
	for (size_t  i=0; i < g_forces.size(); ++i)
	{
		g_forces[i]->ApplyForce();
	}
}

bool TestIntersections(real dt)
{
	for (size_t  j=0; j < g_particles.size(); ++j)
	{
		Particle* particle = g_particles[j];
		
		for (size_t  i=0; i < g_forces.size(); ++i)
		{
			LinearSpring* spring = dynamic_cast<LinearSpring*>(g_forces[i]);
			
			if (spring)
			{
				if (spring->mP == particle || spring->mQ == particle)
					continue;
				
				real t, s;
				if (IntersectPointLineSegmentContinuous(particle->mX, particle->mV, spring->mP->mX, spring->mP->mV, spring->mQ->mX, spring->mQ->mV, t, s, dt))
				{			
					return true;
				}
			}
		}
	}
	
	return false;
					
}

typedef FixedGrid<LinearSpring*> SpringGrid;
SpringGrid* g_grid;

void UpdateGrid(real dt)
{
	//return;
	
	delete g_grid;
	g_grid = new SpringGrid(-kWorldSize, kWorldSize, -kWorldSize, kWorldSize*2.0f, 0.1f, g_forces.size()*10);
	
	for (size_t  i=0; i < g_forces.size(); ++i)
	{
		LinearSpring* spring = dynamic_cast<LinearSpring*>(g_forces[i]);
		
		if (spring)
		{
			Bounds b = spring->CalculateBounds(dt);
		
			g_grid->Insert<false>(spring, b.lower.x, b.upper.x, b.lower.y, b.upper.y);
		}
	}
}


bool Collide(real dt);
void ApplyPenaltyForces(real dt);

void IntegrateSymEuler(real dt)
{
	// apply gravity and damping
	for (size_t  i=0; i < g_particles.size(); ++i)
	{
		Particle& p = *g_particles[i];
		
		p.ApplyImpulse(dt * (p.mF*p.mInvMass));
	}
	
	UpdateGrid(dt);
	
	ApplyPenaltyForces(dt);
	
	int c = 0;
	static int maxc = 0;
	
	while (Collide(dt) && c < 200)
		++c;
	
	maxc = max(c, maxc);
		
	assert(c < 200);
	//assert(TestIntersections(dt) == false);
		

	for (size_t  i=0; i < g_particles.size(); ++i)
	{
		Particle& p = *g_particles[i];
		
		p.mX += dt * (p.mV);
	}
}


Vec2r VectorToLineSegment(const Vec2r& p, const Vec2r& a, const Vec2r& b, real& t)
{
	Vec2r ap = p-a;
	Vec2r ab = b-a;
	t = Dot(ap, ab)/LengthSq(ab);
	Vec2r v = ap - ab*t;
	
	return v;
}

bool TestOverlap(const Vec2r& a, const Vec2r& b, const Vec2r& c, const Vec2r& d, Vec2r& out, real& t, real& u)
{
	Vec2r m = a-c;
	Vec2r cdp = PerpCCW(d-c);
	
	t = -Dot(m, cdp) / Dot(b-a, cdp);

	out = a + t*(b-a);

	
	u = Dot(out-c, (d-c));
	u /= LengthSq(d-c);
	
	return (t > 0.0 && t < 1.0 && u > 0 && u < 1.0);
}

void TestOverlaps()
{
	
	for (size_t  i=0; i < g_forces.size(); ++i)
	{
		LinearSpring* a = dynamic_cast<LinearSpring*>(g_forces[i]);
		
		if (a)
		{
			for (size_t  j=i+1; j < g_forces.size(); ++j)
			{
				LinearSpring* b = dynamic_cast<LinearSpring*>(g_forces[j]);

				if (b)
				{
					//assert(TestOverlap(a->mP->mX, a->mQ->mX, b->mP->mX, b->mQ->mX) == false);
					
					Vec2r i;
					real t, u;
					if (TestOverlap(a->mP->mX, a->mQ->mX, b->mP->mX, b->mQ->mX, i, t, u))
					{
						assert(false);
						
						glBegin(GL_POINTS);
						
						glVertex2x(i);
						glEnd();
					}
				}
			}
		}
	}
}

void ApplyPenaltyForces(real dt)
{	
	const real h = 0.1;
	const real maxf = 0.1f;
	const real k = 10000.0;

	
	// check each particle against the collision planes
	for (size_t  i=0; i < g_particles.size(); ++i)
	{
		Particle& p = *g_particles[i];
		
		for (size_t  j=0; j < g_planes.size(); ++j)
		{
			Vec2r n(g_planes[j].x, g_planes[j].y);
			real d = Dot(p.mX, n) + g_planes[j].z;
			
			if (d < h)
			{				
				float e = h-d;
				
				// normal component of velocity
				real rvn = Dot(p.mV, n);

				if (rvn >= maxf*e/dt)
					continue;
				
				Vec2r i = -min(dt*k*e, maxf*e/dt - rvn)*n;//-std::max(dt*k*d, real(0.0))*n;
				
				// distribute the impulse to the particles
				Vec2r id = i;// / (particle->mInvMass + spring->mP->mInvMass*(s*s) + spring->mQ->mInvMass*(1.0-s)*(1.0-s));
				
				p.ApplyImpulse(-id);

			}
		}
	}
	
	
	std::vector<LinearSpring*> overlaps;
	
	for (size_t  j=0; j < g_particles.size(); ++j)
	{
		Particle* particle = g_particles[j];
	
		Bounds b(particle->mX);
		
		overlaps.resize(0);
		g_grid->Query(b.lower.x-h, b.upper.x+h, b.lower.y-h, b.upper.y+h, overlaps);
		
		for (size_t  i=0; i < overlaps.size(); ++i)
		{
			LinearSpring* spring = (overlaps[i]);
			
			if (spring)
			{
				// calculate vector to spring, s is the poing along 
				// the line segment closest to the particle
				real s;
				Vec2r v = VectorToLineSegment(particle->mX, spring->mP->mX, spring->mQ->mX, s);
				
				
				if (s < 0.0f || s > 1.0f)
					continue;
								
				const real l = Length(v);
				
				if (l > 0.0)
				{
					const real d = h - l;
					if (d <= 0.0)
						continue;

					Vec2r cn = v/l;

					// calculate line segment velocity at point of collision
					Vec2r lv = Lerp(spring->mP->mV, spring->mQ->mV, s);
					
					// relative relocity
					Vec2r rv = particle->mV-lv;
					
					// relative velocity in the normal direction
					real rvn = Dot(rv, cn);
					
					if (rvn >= maxf*d/dt)
						continue;

					// relative tangential velocity
					Vec2r rvt = (rv - rvn*cn)*0.1;
					
					Vec2r i = -min(dt*k*d, maxf*d/dt - rvn)*cn;//-std::max(dt*k*d, real(0.0))*n;
					
					// distribute the impulse to the particles
					Vec2r id = Vec2r(1.0, 1.0) / (particle->mInvMass + spring->mP->mInvMass*(s*s) + spring->mQ->mInvMass*(1.0-s)*(1.0-s));
					
					particle->ApplyImpulse(-i*id);
					
					spring->mP->ApplyImpulse(i*id*(1.0f-s));
					spring->mQ->ApplyImpulse(i*id*s);
					
					particle->ApplyImpulse(-rvt*id);
					
					spring->mP->ApplyImpulse(rvt*id*(1.0f-s));
					spring->mQ->ApplyImpulse(rvt*id*s);
				}
			}
		}
	}
}


bool Collide(real dt)
{	
	bool collided = false;

	const real kRestitution = 0.01;
	
	// check each particle against the collision planes
	for (size_t  i=0; i < g_particles.size(); ++i)
	{
		Particle& p = *g_particles[i];

		for (size_t  j=0; j < g_planes.size(); ++j)
		{
			Vec2r n(g_planes[j].x, g_planes[j].y);
			real d = Dot(p.mX, n);
			real e = d + g_planes[j].z;
			
			if (e < 0.0f)
			{
				// project out
				//p.mX -= n*e;
				
				// normal component of velocity
				real rvn = Dot(p.mV, n);
				if (rvn < 0.0f)
				{
					Vec2r vn = rvn*n;
					
					p.ApplyImpulse(-(1.0f + kRestitution - e)*vn);			
					
					collided = true;
				}
			}
		}
	}
		
	UpdateGrid(dt);
		
	std::vector<LinearSpring*> overlaps;
	
	for (size_t  j=0; j < g_particles.size(); ++j)
	{
		Particle* particle = g_particles[j];

		Bounds b(particle->mX);
		b.Expand(particle->mX + dt*particle->mV);
		
		overlaps.resize(0);
		g_grid->Query(b.lower.x, b.upper.x, b.lower.y, b.upper.y, overlaps);
				
		for (size_t  i=0; i < overlaps.size(); ++i)
		{						
			//LinearSpring* spring = dynamic_cast<LinearSpring*>(overlaps[i]);
			LinearSpring* spring = (overlaps[i]);
			
			if (spring)
			{
				if (spring->mP == particle || spring->mQ == particle)
					continue;
				
				real t, s;
				if (IntersectPointLineSegmentContinuous(particle->mX, particle->mV, spring->mP->mX, spring->mP->mV, spring->mQ->mX, spring->mQ->mV, t, s, dt, false))
				{			
					collided = true;
					
					assert(s >= 0.0f && s <= 1.0f);
										
					// calculate point of collision
					Vec2r cp = particle->mX + t*particle->mV;
					
					Vec2r q = spring->mP->mX+t*spring->mP->mV;
					Vec2r r = spring->mQ->mX+t*spring->mQ->mV;
										
					// calculate collision normal
					Vec2r cn = PerpCCW(Normalize(r-q));
					
					//if (s == 0.0f || s == 1.0f)
					//	cn = SafeNormalize(particle->mX-Lerp(q, r, s), Normalize(q-r));
					
					// calculate line segment velocity at point of collision
					Vec2r lv = Lerp(spring->mP->mV, spring->mQ->mV, s);
					
					// relative relocity
					Vec2r rv = particle->mV-lv;
					
					// relative velocity in the normal direction
					real rvn = Dot(rv, cn);
					
					//assert(rvn <= 0.0f);
					if (rvn > 0.0f)
					{
						cn *= -1.0f;
						rvn *= -1.0f;
					}
					
					//rvn = min(rvn, -0.1);
					
					// impulse in the normal direction
					Vec2r i = real(1.0 + kRestitution)*rvn*cn;
					Vec2r id = i / (particle->mInvMass + spring->mP->mInvMass*(s*s) + spring->mQ->mInvMass*(1.0-s)*(1.0-s));
					
					particle->ApplyImpulse(-id*particle->mInvMass);
					
					spring->mP->ApplyImpulse(id*(1.0f-s)*spring->mP->mInvMass);
					spring->mQ->ApplyImpulse(id*s*spring->mQ->mInvMass);
					
					// calculate new relative velocity
					lv = Lerp(spring->mP->mV, spring->mQ->mV, s);
					rv = particle->mV-lv;
					real newrvn = Dot(rv, cn);

					assert(newrvn > 0.0);
					
					
					//bool draw = (tp == 2066564);
					/*
					bool check = IntersectPointLineSegmentContinuous(particle->mX, particle->mV, spring->mP->mX, spring->mP->mV, spring->mQ->mX, spring->mQ->mV, t, s, dt, draw);
					assert(draw == false);
					if (check)
					{
						lv = Lerp(spring->mP->mV, spring->mQ->mV, s);
						rv = particle->mV-lv;
						newrvn = Dot(rv, cn);
						
						assert(newrvn > 0.0);
						assert(false);
					}
					*/
					
					//Bounds b = spring->CalculateBounds(dt);
					g_grid->Insert<true>(spring, b.lower.x, b.upper.x, b.lower.y, b.upper.y);		
					
					
					//UpdateGrid(dt);
				}
			}
		}
	}
	
	return collided;
}

Particle* GetClosestParticle(const Vec2r& p)
{
	Particle* closest = NULL;
	real closestDistSq = std::numeric_limits<real>::max();
	
	for (size_t  i=0; i < g_particles.size(); ++i)
	{
		real d = LengthSq(g_particles[i]->mX-p);
		
		if (d < closestDistSq)
		{
			closest = g_particles[i];
			closestDistSq = d;
		}
	}
	
	return closest;
}

void TickSim(real dt)
{
	
	dt = 1.0f / 60.0f;	
	//	cout << dt << endl;
	
	const int kSubsteps = 10;
	dt /= kSubsteps;
	
	//random_shuffle(g_forces.begin(), g_forces.end());
	
	for (int i=0; i < kSubsteps; ++i)
	{	
		ApplyForces();
		
		IntegrateSymEuler(dt);
	}
	
	AddSpaghetti(1.0f / 60.0f);
}

void DrawGrid()
{
	glBegin(GL_QUADS);

	float dx = (g_grid->mRight-g_grid->mLeft) / float(g_grid->mWidth);
	float dy = dx;
	
	float left = g_grid->mLeft;
	float bottom = g_grid->mBottom;
		
	for (int y=0; y < g_grid->mHeight; ++y)
	{
		left = g_grid->mLeft;
		
		for (int x=0; x < g_grid->mWidth; ++x)
		{
			if (g_grid->GetCell(x, y) != SpringGrid::kNull)
			{
				glVertex2f(left, bottom);
				glVertex2f(left+dx, bottom);
				glVertex2f(left+dx, bottom+dy);
				glVertex2f(left, bottom+dy);
			}
			
			left += dx;
		}

		bottom += dy;
	}

	glEnd();
}

Vec2r g_pend;
Vec2r g_pstart(0.0f, 1.0f);

bool g_step = false;

void GLUTUpdate()
{
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND); 
	glLineWidth(4.0f);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);	
	glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
	
	glDisable(GL_CULL_FACE);
	glClear(GL_COLOR_BUFFER_BIT);
	glDisable(GL_DEPTH_TEST);
	
	
	glPointSize(5.0f);

	
	real aspect = real(kWidth)/kHeight;
	
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(OrthographicMatrix(-kZoom*aspect, kZoom*aspect, -kZoom, kZoom, 0.0f, 1.0f));

	double t = GetSeconds();
	
	if (1 || g_step)
	{
		TickSim(1.0f/60.0f);
		g_step = false;
		
		TestOverlaps();

	}
	
	real dt = (GetSeconds()-t);
	t = GetSeconds();

	//Vec2r c = g_particles[0]->mX;
	//float zoom = 0.01f;
	//glLoadMatrixf(OrthographicMatrix(c.x-zoom, c.x+zoom, c.y-zoom, c.y+zoom, 0.0f, 1.0f));
				
	glBegin(GL_POINTS);
	glColor3f(1.0f, 1.0f, 1.0f);
	
	for (size_t  i=0; i < g_particles.size(); ++i)
	{
	//	glVertex2x(g_particles[i]->mX);
	}
	
	glEnd();
	
	DrawPlanes();

	//DrawGrid();
	
	// draw forces
	for (size_t  i=0; i < g_forces.size(); ++i)
	{
		g_forces[i]->Draw();
	}
	
	//Vec2r p(0.0f, 2.0f);
	//Vec2r pv(1.0f, -4.0f);
	Vec2r p = g_pstart;
	Vec2r pv = g_pend-p;
	Vec2r q(-1.0f, 0.0f);
	Vec2r qv(0.0f, 1.0f);
	Vec2r r(1.0f, 0.0f);
	Vec2r rv(0.0f, 0.0f);
	
	//real toi, s;
	//IntersectPointLineSegmentContinuous(p, pv, q, qv, r, rv, toi, s);
	
	int x = 10;
	int y = 15;
	
	char line[1024];

	glColor3f(1.0f, 1.0f, 1.0f);
	sprintf(line, "Frame: %.2fms", dt*1000.0f);
	DrawString(x, y, line); y += 13;
	
	sprintf(line, "Particles: %d", int(g_particles.size()));
	DrawString(x, y, line); y += 13;
			   
	sprintf(line, "Forces: %d", int(g_forces.size()));
	DrawString(x, y, line); y += 13;
			   
	glutSwapBuffers();
	
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

void GLUTKeyboardDown(unsigned char key, int x, int y)
{
 	switch (key)
	{
		case 'e':
		{
			break;
		}
		case 'r':
		{
			Reset();
			break;
		}
		case 't':
		{
			AddTri(Vec2r(-1.0f, 1.0f), Randf(0.5f, 2.0f));
			break;
		}
		case 's':
		{
			g_step = true;
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
	}
}

static int lastx;
static int lasty;

Vec2r ScreenToScene(int x, int y)
{
	return Vec2r(-kZoom, -kZoom) + 2.0f*kZoom*Vec2r(real(x)/kWidth, 1.0f-real(y)/kHeight);
}

void GLUTMouseFunc(int b, int state, int x, int y)
{	
	switch (state)
	{
		case GLUT_UP:
		{
			lastx = x;
			lasty = y;
			
			g_pstart = ScreenToScene(x, y);
			
			if (b == GLUT_RIGHT_BUTTON && g_mouseSpring)
			{
				g_forces.erase(find(g_forces.begin(), g_forces.end(), g_mouseSpring));
				delete g_mouseSpring;
			}
			break;
		}
		case GLUT_DOWN:
		{
			lastx = x;
			lasty = y;
			
			if (b == GLUT_RIGHT_BUTTON)
			{
				Particle* p = GetClosestParticle(ScreenToScene(x, y));
				if (p)
				{
					g_mouseSpring = new MouseSpring(p, ScreenToScene(x, y));
					g_forces.push_back(g_mouseSpring);
				}
			}
			else
			{
				break;
				int s = g_particles.size();
				g_particles.push_back(new Particle(g_particles.back()->mX + Vec2r(0.0f, 1.f), 1.0f));

				g_forces.push_back(new LinearSpring(g_particles[s], g_particles[s-1]));
				g_forces.push_back(new AngularSpring(g_particles[s], g_particles[s-1], g_particles[s-2]));
			}
			break;
		}
	}
}

void GLUTMotionFunc(int x, int y)
{
	g_pend = ScreenToScene(x, y);
	
    int dx = x-lastx;
    int dy = y-lasty;
	
	lastx = x;
	lasty = y;
	
	if (g_mouseSpring)
	{
		g_mouseSpring->mMouse = ScreenToScene(x, y);
	}
}

int* topk(int* begin, int* end, int k)
{
	nth_element(begin, end-k, end);
	return end-k;
}

int main(int argc, char* argv[])
{		
	RandInit();
	Init();
	
    // init gl
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH);
	
    glutInitWindowSize(kWidth, kHeight);
    glutCreateWindow("Spaghetti");
    glutPositionWindow(200, 200);
		
    glutMouseFunc(GLUTMouseFunc);
    glutReshapeFunc(GLUTReshape);
    glutDisplayFunc(GLUTUpdate);
    glutKeyboardFunc(GLUTKeyboardDown);
    glutKeyboardUpFunc(GLUTKeyboardUp);
    glutIdleFunc(GLUTUpdate);	
    glutSpecialFunc(GLUTArrowKeys);
    glutSpecialUpFunc(GLUTArrowKeysUp);
    glutMotionFunc(GLUTMotionFunc);

#ifndef WIN32
	int swap_interval = 1;
	CGLContextObj cgl_context = CGLGetCurrentContext();
	CGLSetParameter(cgl_context, kCGLCPSwapInterval, &swap_interval);
#endif

    glutMainLoop();

}

