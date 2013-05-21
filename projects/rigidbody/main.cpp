#include <core/types.h>
#include <core/maths.h>
#include <core/platform.h>

#include <graphics/rendergl/glutil.h>

#ifndef WIN32
#include <opengl/opengl.h>
#endif

#include <iostream>

using namespace std;

const uint32_t kHeight = 600;
const uint32_t kWidth = 600;
const float kWorldSize = 4.0f;
const float kZoom = kWorldSize + 0.0f;

struct Particle
{
	Particle(Vec2 p, Vec2 v, float r) : position(p), velocity(v), radius(r) {}
	
	Vec2 position;
	Vec2 velocity;
	Vec2 force;
	float radius;
};

std::vector<Vec3> g_planes;
std::vector<Particle> g_particles;

Vec2 g_gravity(0.0f, -9.8f);


void Init()
{
	g_planes.push_back(Vec3(0.0f, 1.0f, 0.0f));
	g_planes.push_back(Vec3(-1.0f, 0.0f, -2.0f));
	g_planes.push_back(Vec3(1.0f, 0.0f, -2.0f));	

	// fixed world particle
	g_particles.push_back(Particle(Vec2(), Vec2(), 0.0f));
	
	float radius = 0.2f;

	for (int i=0; i < 50; ++i)
		 g_particles.push_back(Particle(Vec2(0.0f, radius + i*2.0f*radius), Vec2(0.0f, 0.0f), radius));
	
	//g_particles[1].velocity += Vec2(1.0f, 0.0f);
}

void ApplyForces(float dt)
{
	for (int i=1; i < g_particles.size(); ++i)
	{
		Particle& p = g_particles[i];
		
		p.force = g_gravity;//f*(p.velocity);
	}
}

struct Contact
{
	Contact() : x(0.0f, 0.0f), n(0.0f, 0.0f), j(0.0f), tj(0.0f), d(0.0f), a(0), b(0) {}
	Vec2 x;
	Vec2 n;
	
	float j;
	float tj;
	float d;
	
	int a;
	int b;
};


void CollidePlanes(const std::vector<Vec3>& planes, const std::vector<Particle>& particles, std::vector<Contact>& contacts)
{
	// collide with planes
	for (int i=0; i < planes.size(); ++i)
	{
		Vec3 p = planes[i];
		
		for (int r=1; r < particles.size(); ++r)
		{
			Vec2 x = particles[r].position;			
			float radius = particles[r].radius;
			
			// distance to plane
			float d = x.x*p.x + x.y*p.y - p.z;
			
			float mtd = d - radius;
						
			if (mtd < 0.0f)
			{
				Contact c;
				c.n = Vec2(p.x, p.y);
				c.d = mtd;
				c.a = r;
				c.b = 0;
				
				contacts.push_back(c);
			}
		}
	}			
}

// takes an array of particles and planes and returns a list of contacts
void CollideParticles(const std::vector<Particle>& particles, std::vector<Contact>& contacts)
{
	// collide with particles
	for (int i=1; i < particles.size(); ++i)
	{
		const Particle& a = particles[i];
		
		for (int j=i+1; j < particles.size(); ++j)
		{
			const Particle& b = particles[j];
			
			Vec2 pd = a.position - b.position;
			
			float l = Length(pd);
			float mtd = l - a.radius - b.radius;

/*
void SolveImpulse(const std::vector<Contact>& contacts, float dt)
{
	// apply penalty forces for each contact
	for (int i=0; i < contacts.size(); ++i)
	{
		const Contact& c = contacts[i];
		
		assert(c.d < 0.0f);
		
		float me = 1.0f / (c.a->invmass + c.b->invmass);
		
		Vec2 va = c.a->GetVelocityAtPoint(c.x);
		Vec2 vb = c.b->GetVelocityAtPoint(c.x);
		
		float rvn = Dot(c.n, va-vb);
		
		float fdamp = 0.004f*k*rvn;
		
		// apply a penalty force to bodies
		c.a->ApplyForce(-c.n*(f+fdamp), c.x);
		c.b->ApplyForce( c.n*(f+fdamp), c.x);
	}
}
*/

class Matrix
{
public:

			if (mtd < 0.0f)
			{
				Contact c;
				c.n = n;
				c.d = mtd;
				c.a = i;
				c.b = j;
				
				contacts.push_back(c);
			}
		}
	}
}

void Solve(std::vector<Contact>& contacts, std::vector<Contact>& prevContacts, std::vector<Particle>& particles, float dt)
{	
	const float kOverlap = 0.01f;
	const float kBaumgarte = 0.2f / dt;
	const int   kIterations = 10;
	
	if (1)
	{
		for (int i=0; i < contacts.size(); ++i)
		{
			Contact& c = contacts[i];

			for (int j=0; j < prevContacts.size(); ++j)
			{
				Contact& pc = prevContacts[j];
				
				if (pc.a == c.a && pc.b == c.b)
				{
					// apply previous frames impulse
					
					Particle& a = particles[c.a];
					Particle& b = particles[c.b];
					
					// treat radius as 1.0/mass
					float ma = c.a?1.0f:0.0f;
					float mb = c.b?1.0f:0.0f;				
					float msum = ma + mb;
					
					assert(msum > 0.0f);
					
					// impulse
					a.velocity -= pc.tj*ma*c.n/msum;
					b.velocity += pc.tj*mb*c.n/msum;
					
					c.tj = pc.tj;
				}
			}
		}
	}
	
	for (int s=0; s < kIterations; ++s)
	{		
		for (int i=0; i < contacts.size(); ++i)
		{
			Contact& c = contacts[i];
				
			Particle& a = particles[c.a];
			Particle& b = particles[c.b];
			
			Vec2 vd = a.velocity-b.velocity;
						
			// calculate relative normal velocity
			float vn = Dot(vd, c.n);

			// Baumgarte stabilisation
			float bias = kBaumgarte*min(c.d+kOverlap, 0.0f);
			
			float l = vn + bias;
			
			float j = min(c.tj + l, 0.0f);
			l = j - c.tj;
			c.tj = j;
			
			c.j = l; 
			//c.tj += c.j;
			/*
			
		}

		// update all particles velocities
		for (int i=0; i < contacts.size(); ++i)
		{
			Contact& c = contacts[i];
			
			Particle& a = particles[c.a];
			Particle& b = particles[c.b];
*/
			
			// treat radius as 1.0/mass
			float ma = c.a?1.0f:0.0f;
			float mb = c.b?1.0f:0.0f;				
			float msum = ma + mb;
			
			assert(msum > 0.0f);
			
			// impulse
			a.velocity -= c.j*ma*c.n/msum;
			b.velocity += c.j*mb*c.n/msum;
			
			//c.j = 0.0f;
		}	
	}
}


void Integrate(float dt)
{
	// v += a*dt
	for (int i=1; i < g_particles.size(); ++i)
	{
		Particle& p = g_particles[i];
		
		p.velocity += p.force*dt;
	}
	
	std::vector<Contact> contacts;
	
	CollideParticles(g_particles, contacts);
	CollidePlanes(g_planes, g_particles, contacts);

	static std::vector<Contact> prevContacts;
	
	// solve
	Solve(contacts, prevContacts, g_particles, dt);
	
	swap(contacts, prevContacts);
	
	// p += v*dt
	for (int i=1; i < g_particles.size(); ++i)
	{
		Particle& p = g_particles[i];
		
		p.position += p.velocity*dt;
	}
	
	static int counter = 0;
	counter++;
	
	// calc error
	float error = fabsf(g_particles.back().position.y - (1 + 2*(g_particles.size()-1))*g_particles.back().radius);
	static float maxerror;
	maxerror = max(error, maxerror);

	cout << error << " max=" << maxerror << endl;
	
	//if (counter == 300)
	//	exit(0);
}

//	printf("\nA:\n");
//	A.print();

	// solve
	Matrix x = GaussSeidel(A, b, 30, false);
	
	CollideParticles(g_particles, contacts);
	CollidePlanes(g_planes, g_particles, contacts);

	for (int i=0; i < contacts.size(); ++i)
	{
		Contact& c = contacts[i];
		
		const float kStiff = 20000.0f;
		const float kDamp = 100.0f;
		
		float rv = Dot(g_particles[c.a].velocity - g_particles[c.b].velocity, c.n);
	
		if (rv < 0.0f)
		{
			g_particles[c.a].force -= (kStiff*c.d + kDamp*rv)*c.n;
			g_particles[c.b].force += (kStiff*c.d + kDamp*rv)*c.n;
		}
		
	}
	
	for (int i=1; i < g_particles.size(); ++i)
	{
		Particle& p = g_particles[i];
		
		p.velocity += p.force*dt;
		p.position += p.velocity*dt;
	}
	
	static int counter = 0;
	counter++;
	
	// calc error
	float error = fabsf(g_particles.back().position.y - (1 + 2*(g_particles.size()-1))*g_particles.back().radius);
	static float maxerror;
	maxerror = max(error, maxerror);
	
	RigidBody* b = new RigidBody(Vec2(x, y), radius, 1.0f, &points[0], &points[0] + count);
	g_bodies.push_back(b);
}

void CreateBalls()
{
	 const float kRadius = 0.05f;
	 const int kBalls = 40;
	 
	 Vec2 points(0.0f);
	 
	 for (int i=0; i < kBalls; ++i)
	 {
		g_bodies.push_back(new RigidBody(Vec2(0.0f, kRadius + kRadius * (2.0f*i)), kRadius, 1.0f, &points, &(points)+1));
	 }
}

void Init()
{
	g_planes.resize(0);
	g_forces.resize(0);
	g_bodies.resize(0);
	g_mouseSpring = NULL;
	
}

void TickSim(float dt)
{
	dt = 1.0f / 30.0f;	
	//	cout << dt << endl;
	
	const int kSubsteps = 1;
	dt /= kSubsteps;
	
	for (int i=0; i < kSubsteps; ++i)
	{		
		ApplyForces(dt);
		Integrate(dt);
	}
}

void DrawCircle(const Vec2& c, float r)
{
	glBegin(GL_TRIANGLE_FAN);
	glVertex2fv(c);
	
	const int kSegments = 40;
	for (int i=0; i < kSegments+1; ++i)
	{
		float theta = k2Pi*float(i)/kSegments;
		
		float y = c.y + r*Cos(theta);
		float x = c.x + r*Sin(theta);
		
		glVertex2f(x, y);		
	}
	
	glEnd();
	
}

void DrawString(int x, int y, const char* s)
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, kWidth, kHeight, 0);
	
	glRasterPos2d(x, y);
	while (*s)
	{
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *s);
		++s;
	}
}

bool g_step = false;

void GLUTUpdate()
{
	//glEnable(GL_LINE_SMOOTH);
	//glEnable(GL_BLEND); 
	//glLineWidth(1.0f);
	//glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);	
	//glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
	
	glDisable(GL_CULL_FACE);
	glClear(GL_COLOR_BUFFER_BIT);
	glDisable(GL_DEPTH_TEST);
	
	glPointSize(5.0f);

	float aspect = float(kWidth)/kHeight;
	
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(OrthographicMatrix(-kZoom*aspect, kZoom*aspect, -kZoom, kZoom, 0.0f, 1.0f));

	double t = GetSeconds();
	
	if (1 || g_step)
	{
		TickSim(1.0f/60.0f);
		g_step = false;
	}
	
	float dt = GetSeconds()-t;
		
	for (int i=0; i < g_planes.size(); ++i)
	{	
		Vec2 p = g_planes[i].z * Vec2(g_planes[i].x, g_planes[i].y);
		Vec2 d = Vec2(-g_planes[i].y, g_planes[i].x);
		
		glBegin(GL_LINES);
		glColor3f(1.0f, 1.0f, 1.0f);
		glVertex2fv(p - d*1000.0);
		glVertex2fv(p + d*1000.0);
		glEnd();		
	}
	
	for (int i=0; i < g_particles.size(); ++i)
	{
		DrawCircle(g_particles[i].position, g_particles[i].radius);
	}
		
	int x = 10;
	int y = 15;
	
	char line[1024];

	glColor3f(1.0f, 1.0f, 1.0f);
	sprintf(line, "Frame: %.2fms", dt*1000.0f);
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
			Init();
			break;
		}
		case 't':
		{
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

Vec2 ScreenToScene(int x, int y)
{
	return Vec2(-kZoom, -kZoom) + 2.0f*kZoom*Vec2(float(x)/kWidth, 1.0f-float(y)/kHeight);
}

void GLUTMouseFunc(int b, int state, int x, int y)
{	
	switch (state)
	{
		case GLUT_UP:
		{
			lastx = x;
			lasty = y;
			
			break;
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
}


int main(int argc, char* argv[])
{	
	RandInit();
	Init();
	
    // init gl
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH);
	
    glutInitWindowSize(kWidth, kHeight);
    glutCreateWindow("Granular");
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

