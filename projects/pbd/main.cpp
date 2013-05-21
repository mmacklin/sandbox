#include <iostream>

#include "core/types.h"
#include "core/shader.h"

using namespace std;

uint32_t g_screenWidth = 600;
uint32_t g_screenHeight = 600;

struct Particle
{
	Vec2 x0,x1;
	Vec2 v0;
};

Particle g_particle;
FILE* g_file;
int g_frame;

bool g_step = false;
bool g_pause = false;

float sqr(float x) { return x*x; }

void Init()
{
	g_file = fopen("dump.txt", "w");
	fprintf(g_file, "{");

	g_particle.x0 = Vec2(1.0f, 0.0f);
	g_particle.x1 = Vec2(1.0f, 0.0f);
	g_particle.v0 = Vec2(0.0f);
}

void GLUTUpdate()
{	
	glViewport(0, 0, g_screenWidth, g_screenHeight);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluPerspective(45.0f, g_screenWidth/float(g_screenHeight), 0.01f, 1000.0f);
	gluOrtho2D(-2.0, 2.0, -2.0f, 2.0f);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	//gluLookAt(sin(g_angle)*2.0f, 0.4f, cos(g_angle)*2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
	float dt = 1.0f/60.0f;

	Vec2 v = g_particle.v0 + 0.5f*Vec2(0.0f, -9.8f)*dt + 0.5f*Vec2(0.0f, -9.8f)*sqr(0.5f*dt);
	
	Vec2 newx = g_particle.x0 + v*dt; 

	// constrain
	Vec2 j = newx - Vec2(0.0f);
	float e = Length(j)-1.0f;
	newx -= Normalize(j)*e;

	
/*
	float k0 = 1.0f;
	float k1 = -1.0f;
	float k2 = 0.0f;
	*/

	//float k0 = 3.0f/2.0f;
	//float k1 = -2.0f;
	//float k2 = 1.0f/2.0f;
	//float k0 = 2.5f + 2.0f*sqrtf(2.0f);
    //float k2 = 1.5f + sqrtf(2.0f);	
	//float k1 = -(k0 + k2);  
	
	float alpha = 0.5f;//2.0f - sqrtf(2.0f);
	float k0 = 2.0f-alpha;
	float k1 = 1.0f/alpha;
	float k2 = sqr(1.0f-alpha)/alpha;


	glPointSize(5.0f);
	glBegin(GL_POINTS);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex2fv(g_particle.x0);	
	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex2fv(g_particle.x1);	
	glColor3f(1.0f, 1.0f, 1.0f);
	glVertex2fv(newx);	
	glEnd();	


	if (!g_pause || g_step)
	{
		float ke0 = (g_particle.x0.y+1.0f)*9.8f + Dot(g_particle.v0, g_particle.v0)*0.5f;

		//g_particle.v0 = (newx-g_particle.x0)/dt;
		g_particle.v0 = (newx*k0 - g_particle.x0*k1 + g_particle.x1*k2)/dt;
		//g_particle.v0 = (newx - g_particle.x1)/(2.0f*dt);
		//g_particle.v0 = (g_particle.x0-g_particle.x1)/dt;
		//g_particle.v0 = (newx-g_particle.x0)/dt + (newx-2.0f*g_particle.x0+g_particle.x1)/(2.0f*dt);

		//Vec2 dx = (newx-g_particle.x0)/dt;
		//Vec2 ddx = (-2.0f*g_particle.x0 + newx + g_particle.x1)/(dt*dt);
		//g_particle.v0 = dx + ddx*dt*0.3f;
	
		g_particle.x1 = g_particle.x0;
		g_particle.x0 = newx;

		/*	
		// make velocity energy conserving
		float ke1 = (g_particle.x0.y+1.0f)*9.8f + Dot(g_particle.v0, g_particle.v0)*0.5f;

		float a = 0.5f*Dot(g_particle.v0, g_particle.v0);
	    float b = Dot(g_particle.v0, g_particle.v0);
		float c = -(ke0-ke1);

		float e1, e2;
		SolveQuadratic(a, b, c, e1, e2);

		g_particle.v0 *= (1.0f + e2);
		*/
		
	}

	float ke = (g_particle.x0.y+1.0f)*9.8f + Dot(g_particle.v0, g_particle.v0)*0.5f;
	DrawString(0.0, 0.0, "%f", ke); 

	fprintf(g_file, "%f,", g_particle.x0.y);
	
	if (g_frame++ == 1000)
	{
		fprintf(g_file, "0.0}\n");
		fclose(g_file);

		exit(0);
	}

	g_step = false;


	// flip
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
		case ' ':
		{
			g_step = true;
			break;
		}
		case 'p':
		{
			g_pause = !g_pause;
			break;
		}
		case 'r':
		{
			Init();
			break;
		}
		case 'q':
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
			break;
		}
		case ' ':
		{
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
	
}

void GLUTPassiveMotionFunc(int x, int y)
{
    int dx = x-lastx;
    int dy = y-lasty;
	
	lastx = x;
	lasty = y;

}


int main(int argc, char* argv[])
{	
    // init gl
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH);
	
    glutInitWindowSize(g_screenWidth, g_screenHeight);
    glutCreateWindow("Empty");
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

