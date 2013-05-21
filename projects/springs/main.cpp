#include <core/maths.h>
#include <core/shader.h>
#include <core/platform.h>
#include <core/mat22.h>
#include "cg.h"

#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <vector>
#include <stdint.h>

#include <xmmintrin.h>

using namespace std;

const int kWidth = 800;
const int kHeight = 600;

float gViewLeft = -1.0f;
float gViewBottom = -2.0f;
float gViewWidth = 5.0f;
float gViewAspect = kHeight/float(kWidth);

Vec2 gGravity(0.0f, -9.8f);
float gStiffness = 100000.0f; 
float gDamping = 0.0f;
const int gSubsteps = 1;

Vec2 gMousePos;
int gMouseIndex=-1;
float gMouseStrength = 100.0f;

const Matrix22 kIdentity(Vec2(1.0f, 0.0f), Vec2(0.0f, 1.0f));

struct Particle
{
	Particle(Vec2 x, float im) : p(x), invMass(im)  {}

	Vec2 p;
	Vec2 v;
	Vec2 f;
	float invMass;
};

struct Spring
{
	uint32_t i, j;
	float r;
	float ks;
	float kd;
};

std::vector<Particle> gParticles;
std::vector<Spring> gSprings;

int FindClosestParticle(Vec2 p)
{
	float minDistSq = FLT_MAX;
	int minIndex = -1;

	for (size_t i=0; i < gParticles.size(); ++i)
	{
		float d = LengthSq(p-gParticles[i].p);
		
		if (d < minDistSq)
		{
			minDistSq = d;
			minIndex = i;
		}		
	}
	return minIndex;
}

void Init()
{
	gParticles.clear();
	gSprings.clear();

	/* Create chain */
	for (int i=0; i < 10; ++i)
	{
		const float kLength = 0.1f;

		gParticles.push_back(Particle(Vec2(i*kLength, 0.0f), i?1.0f:0.0f));
		
		if (i)
		{
			Spring s = { i-1, i, kLength, gStiffness, gDamping };
			gSprings.push_back(s);
		}
	}
}

void SolveExplicit(float dt)
{
	// integrate particles forward in time 	
	for (size_t i=0; i < gParticles.size(); ++i)
	{
		gParticles[i].v += gParticles[i].invMass*(gParticles[i].f - gDamping*gParticles[i].v)*dt;
		gParticles[i].p += gParticles[i].v*dt;// + 0.5f*gParticles[i].f*dt*dt*gParticles[i].invMass;
	}
}

Matrix22 SpringJdx(uint32_t i, uint32_t j, float r, float ks)
{
	Vec2 p = gParticles[i].p;
	Vec2 q = gParticles[j].p;

	Vec2 u = p-q;
	float m = Length(u);

	// normalize
	u /= m;

	Matrix22 l = Outer(u, u); 

	Matrix22 c = -ks*(l + (1.0f - r/m)*(kIdentity-l));
	return c;
}

Matrix22 SpringJdv(uint32_t i, uint32_t j, float kd)
{
	Vec2 p = gParticles[i].p;
	Vec2 q = gParticles[j].p;

	Vec2 u = Normalize(p-q);

	return -kd*Outer(u, u);
}

// returns the force on i from spring
Vec2 SpringF(uint32_t i, uint32_t j, float r, float ks, float kd, Vec2 dp, Vec2 dq, Vec2 dvi=Vec2(), Vec2 dvj=Vec2())
{
	Vec2 p = gParticles[i].p+dp;
	Vec2 q = gParticles[j].p+dq;

	Vec2 u = p-q;
	Vec2 v = gParticles[i].v+dvi-gParticles[j].v+dvj; 

	float l = Length(u);
	u /= l;

	float e = l-r;

	return -(ks *e + kd*Dot(v, u))*u;
}


// finite difference approximation to J
Matrix22 SpringJn(uint32_t i, uint32_t j, float r, float ks, float kd, float dp, float dq)
{
	Vec2 fdx = SpringF(i, j, r, ks, kd, Vec2(dp, 0.0f), Vec2(dq, 0.0f)) - SpringF(i, j, r, ks, kd, Vec2(), Vec2());
	Vec2 fdy = SpringF(i, j, r, ks, kd, Vec2(0.0f, dp), Vec2(0.0f, dq)) - SpringF(i, j, r, ks, kd, Vec2(), Vec2());

	if (dp > 0.0f)	
		return Matrix22(fdx/dp, fdy/dp);	
	else
		return Matrix22(fdx/dq, fdy/dq);
}

void SolveImplicit(float dt)
{
	size_t n = gParticles.size();

	// construct n by n block matrix of force Jacobians
	std::vector<NodeRow> A(n, NodeRow(n));
	std::vector<NodeRow> dJdx(n, NodeRow(n));

	for (size_t i=0; i < gSprings.size(); ++i)
	{
		const Spring& s = gSprings[i];

		// calculate Jacobian of force on i
		Matrix22 Jp = -dt*dt*SpringJdx(s.i, s.j, s.r, s.ks);

		Matrix22 Jpdp = Jp;
		Matrix22 Jpdq = -1.0f*Jp;
		Matrix22 Jqdq = Jp;
		Matrix22 Jqdp = -1.0f*Jp;
	
		Matrix22 Jv = -dt*SpringJdv(s.i, s.j, s.kd);

		// force derivative wrt first particle
		A[s.i][s.i] += Jpdp+Jv;

		// force on i due to other particle is in the opposite direction
		A[s.i][s.j] += Jpdq-Jv;

		// force on j due to j is the same as i
		A[s.j][s.j] += Jqdq+Jv;

		// force on j due to i is same as force on i due to j 
		A[s.j][s.i] += Jqdp-Jv;

	//	-----------------------------------
	//
		// force derivative wrt first particle
		dJdx[s.i][s.i] += Jpdp;

		// force on i due to other particle is in the opposite direction
		dJdx[s.i][s.j] += Jpdq;

		// force on j due to j is the same as i
		dJdx[s.j][s.j] += Jqdq;

		// force on j due to i is same as force on i due to j 
		dJdx[s.j][s.i] += Jqdp;
	}	

	std::vector<Vec2> v(n);
	std::vector<Vec2> f(n);

	for (size_t i=0; i < n; ++i)
	{
		f[i] = dt*gParticles[i].f;
		v[i] = gParticles[i].v;
	}

	// calculate b
	std::vector<Vec2> b = CgSub(f, CgMul(dJdx, v));

	// multiply through by mass
	for (size_t i=0; i < n; ++i)
		A[i][i] = 1.0f/max(gParticles[i].invMass, 0.01f)*kIdentity + A[i][i];

	// solve for dv
	std::vector<Vec2> dv = CgSolve(A, b, 20, 0.001f);

	// update v
	for (size_t i=0; i < n; ++i)
	{
		gParticles[i].v += dv[i]*gParticles[i].invMass;
		gParticles[i].p += gParticles[i].v*dt;
	}	
}

void Advance(float dt)
{
	static float et = 0.0f;
	et += dt;

	for (size_t i=0; i < gParticles.size(); ++i)
	{
		Particle& p = gParticles[i];

		//if (p.invMass > 0.0f)
			p.f = gGravity;//*(1.0f/p.invMass);
	}

	for (size_t i=0; i < gSprings.size(); ++i)
	{
		Spring& s = gSprings[i];

		Vec2 f = SpringF(s.i, s.j, s.r, s.ks, s.kd, Vec2(), Vec2());

		
		float msum = gParticles[s.i].invMass + gParticles[s.j].invMass;
		float ma = gParticles[s.i].invMass / msum;
	    float mb = gParticles[s.j].invMass / msum;

		gParticles[s.i].f += f;
		gParticles[s.j].f -= f;
	}

	if (gMouseIndex != -1)
		gParticles[gMouseIndex].f += gMouseStrength*(gMousePos-gParticles[gMouseIndex].p);

	//SolveExplicit(dt);
	SolveImplicit(dt);
}


void Update()
{
	const float dt = 1.0f/60.0f;

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDisable(GL_DEPTH_TEST);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(gViewLeft, gViewLeft+gViewWidth, gViewBottom, gViewBottom+gViewWidth*gViewAspect);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	double startTime = GetSeconds();

	for (int i=0; i < gSubsteps; ++i)
	{	
		Advance(dt/gSubsteps);
	}
	
	double elapsedTime = GetSeconds()-startTime;

	glBegin(GL_LINES);

	for (size_t i=0; i < gSprings.size(); ++i)
	{
		Vec2 a = gParticles[gSprings[i].i].p;
		Vec2 b = gParticles[gSprings[i].j].p;

		glVertex2fv(a);
		glVertex2fv(b);	
	}

	glEnd();

	glPointSize(4.0f);
	glBegin(GL_POINTS);
	glColor3f(0.0f, 1.0f, 0.0f);

	for (size_t i=0; i < gParticles.size(); ++i)
	{
		glVertex2fv(gParticles[i].p);
	}

	glEnd();

	if (gMouseIndex != -1)
	{
		glBegin(GL_LINES);
		glColor3f(0.0f, 1.0f, 0.0f);
		glVertex2fv(gMousePos);
		glVertex2fv(gParticles[gMouseIndex].p);
		glEnd();
	}
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, kWidth, kHeight, 0);

	int x = 10;
	int y = 15;
	
	char line[1024];

	glColor3f(1.0f, 1.0f, 1.0f);
	sprintf(line, "Time: %.2fms", float(elapsedTime)*1000.0f);
	DrawString(x, y, line); y += 13;

	glutSwapBuffers();
}

Vec2 RasterToScene(int x, int y)
{
	float vx = gViewLeft + gViewWidth*x/float(kWidth);
	float vy = gViewBottom + gViewWidth*gViewAspect*(1.0f-y/float(kHeight));
	
	return Vec2(vx, vy);
}

int lastx = 0;
int lasty = 0;

void GLUTMouseFunc(int b, int state, int x, int y)
{	
	switch (state)
	{
		case GLUT_UP:
		{
			lastx = x;
			lasty = y;

			gMouseIndex = -1;
			
			break;
		}	
		case GLUT_DOWN:
		{
			lastx = x;
			lasty = y;

			gMousePos = RasterToScene(x, y);
			gMouseIndex = FindClosestParticle(gMousePos);
		}
	};
}

void GLUTKeyboardDown(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'w':
		{
			break;
		}
		case 's':
		{
			break;
		}
		case 'a':
		{
			break;
		}
		case 'd':
		{
			break;
		}
		case 'u':
		{
			break;
		}
		case 'j':
		{
			break;
		}
		case 'r':
		{
			Init();
			break;
		}
		case 27:
		case 'q':
		{
			exit(0);
			break;
		}
	};
}

void GLUTKeyboardUp(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'w':
		{
			break;
		}
		case 's':
		{
			break;
		}
		case 'a':
		{
			break;
		}
		case 'd':
		{
			break;
		}
	};
}

void GLUTMotionFunc(int x, int y)
{	
//	int dx = x-lastx;
//	int dy = y-lasty;
	
	lastx = x;
	lasty = y;
	
	gMousePos = RasterToScene(x, y);
}

int main(int argc, char* argv[])
{
	//_MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);		
		
	// init gl
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH);
	
	glutInitWindowSize(kWidth, kHeight);
	glutCreateWindow("FEM");
	glutPositionWindow(200, 100);

#if _WIN32
	glewInit();
#endif

	Init();

	glutIdleFunc(Update);	
	glutDisplayFunc(Update);
	glutMouseFunc(GLUTMouseFunc);
	glutMotionFunc(GLUTMotionFunc);
	glutKeyboardFunc(GLUTKeyboardDown);
	glutKeyboardUpFunc(GLUTKeyboardUp);
/*	

	glutReshapeFunc(GLUTReshape);
	glutDisplayFunc(GLUTUpdate);
   
	glutSpecialFunc(GLUTArrowKeys);
	glutSpecialUpFunc(GLUTArrowKeysUp);

*/
#if __APPLE__
	int swap_interval = 1;
	CGLContextObj cgl_context = CGLGetCurrentContext();
	CGLSetParameter(cgl_context, kCGLCPSwapInterval, &swap_interval);
#endif

	glutMainLoop();
	return 0;
}





