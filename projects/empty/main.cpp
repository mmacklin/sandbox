#include <iostream>

#include "core/types.h"
#include "core/shader.h"

using namespace std;

uint32_t g_screenWidth = 1280;
uint32_t g_screenHeight = 720;

void Init()
{

}

void GLUTUpdate()
{	
	glViewport(0, 0, g_screenWidth, g_screenHeight);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0f, g_screenWidth/float(g_screenHeight), 0.01f, 1000.0f);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	//gluLookAt(sin(g_angle)*2.0f, 0.4f, cos(g_angle)*2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
	
	
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

