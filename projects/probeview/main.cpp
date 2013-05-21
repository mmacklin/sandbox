#include <iostream>

#include "core/types.h"
#include "core/aabbtree.h"
#include "core/mesh.h"
#include "core/pfm.h"

#include "graphics/rendergl/glutil.h"

using namespace std;

const uint32_t g_screenWidth = 1024;
const uint32_t g_screenHeight = 512;

Point3 g_camPos;
Rotation g_camDir;

PfmImage g_probe;

uint32_t g_buffer[g_screenWidth*g_screenHeight];

Colour SampleProbe(const Vector3& dir, const PfmImage& image)
{
	// convert world space dir to probe space
	float c = (1.0f / kPi) * acosf(dir.z)/sqrt(dir.x*dir.x + dir.y*dir.y);
	
	uint32_t px = (0.5f + 0.5f*(dir.x*c))*g_probe.m_width;
	uint32_t py = (0.5f + 0.5f*(-dir.y*c))*g_probe.m_height;
	
	float r = g_probe.m_data[py*g_probe.m_width*3 + (px*3)+0];
	float g = g_probe.m_data[py*g_probe.m_width*3 + (px*3)+1];
	float b = g_probe.m_data[py*g_probe.m_width*3 + (px*3)+2];

	return Colour(r, g, b, 1.0f);
}

void Init()
{
	//const char* probeFile = "../../probes/grace_probe.pfm.pfm";
	const char* probeFile = "../../probes/uffizi_probe.pfm.pfm";
	
	if (!PfmLoad(probeFile, g_probe))
	{
		cout << "Couldn't load probe\n" << endl;
		exit(-1);
	}
}


void GLUTUpdate()
{	
	GlVerify(glEnable(GL_CULL_FACE));
	GlVerify(glEnable(GL_DEPTH_TEST));
	GlVerify(glDisable(GL_LIGHTING));
	GlVerify(glDisable(GL_BLEND));
	
	GlVerify(glViewport(0, 0, g_screenWidth, g_screenHeight));
	GlVerify(glClearColor(0.5f, 0.5f, 0.5f, 1.0f));
	GlVerify(glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT));
	
		
	Matrix44 rasterToScreen( 2.0f / g_screenWidth, 0.0f, 0.0f, 0.0f,
							0.0f, -2.0f / g_screenHeight, 0.0f, 0.0f,
							0.0f,  0.0f, 1.0f, 0.0f,
							-1.0f,  1.0f, 1.0f, 1.0f);
	
	float f = tanf(DegToRad(45.0f)*0.5f);
	float aspect = float(g_screenWidth) / g_screenHeight;
	
	Matrix44 screenToCamera(f*aspect, 0.0f, 0.0f, 0.0f,
							0.0f, f, 0.0f, 0.0f, 
							0.0f, 0.0f, -1.0f, 0.0f,
							0.0f, 0.0f, 0.0f, 1.0f);

	Matrix44 cameraToWorld = TransformMatrix(g_camDir, Point3(0.0f));
	
	Matrix44 rasterToWorld = cameraToWorld*screenToCamera*rasterToScreen;
	
	for (uint32_t y=0; y < g_screenHeight; ++y)
	{
		for (uint32_t x=0; x < g_screenWidth; ++x)
		{
			Point3 p = rasterToWorld * Point3(float(x) + 0.5f, float(y) + 0.5f, 0.0f);		
			Vector3 dir = Normalize(p-Point3(0.0f));
						
			g_buffer[y*g_screenWidth + x] = ColourToRGBA8(LinearToSrgb(ToneMap(SampleProbe(dir, g_probe))));
		}
	}
	    
    GlVerify(glDrawPixels(g_screenWidth, g_screenHeight, GL_RGBA, GL_UNSIGNED_BYTE, g_buffer));
	
	// flip
	GlVerify(glutSwapBuffers());
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
	
    const float sensitivity = 0.1f;
	
    g_camDir.yaw -= dx*sensitivity;
    g_camDir.roll += dy*sensitivity;
	
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
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	
    glutInitWindowSize(g_screenWidth, g_screenHeight);
    glutCreateWindow("ProbeView");
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

