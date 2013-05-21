// Main.cpp : Defines the entry point for the console application.
//
#include <vector>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>

#include "stdlib.h" // for atexit()

#include "maths.h"
#include "scene.h"
#include "particlecontainer.h"
#include "particleemitter.h"
#include "camera.h"
#include "view.h"
#include "texture.h"
#include "illuminationsim.h"

using namespace std;

int g_width = 512;
int g_height = 512;

Scene* g_scene;
Camera* g_camera;
View* g_view;
Mesh* g_mesh;
IlluminationSim* g_sim;
IlluminationSim::DrawMode g_drawmode = IlluminationSim::kDrawMesh;
int g_drawdepth = 0;
Surfel* g_light = NULL;

bool g_enable = false;

void DrawText(int x, int y, const char *string);

void DeInit()
{
}

void Init(int argc, char* argv[])
{
	g_scene = new Scene();

	// tracking camera
	g_camera = new Camera();
	g_camera->SetPosition(Point3(0, 60, 250));
	g_camera->SetRotation(Rotation(0, 0, -DegToRad(10.0f)));

	// main scene view
	g_view = new View();
	g_view->SetCamera(g_camera);
	g_view->m_x = 0;
	g_view->m_y = 0;
	g_view->m_width = g_width;
	g_view->m_height = g_height;

	// load test mesh
	g_mesh = new Mesh("cornell.dae");

	// create illumination sim
	g_sim = new IlluminationSim(g_mesh);
/*
	g_light = g_sim->AddLight();
	g_light->position = Vec3(0,0,0);
	g_light->area = 10.0f;
	g_light->normal = Vec3(0,-1,0);
	g_light->colour = Colour(1,1,1,1);
	g_light->emissive = Colour(50, 50, 50, 1.0f);
*/
}

void Reshape(int width, int height)
{
	g_view->m_width = width;
	g_view->m_height = height;

	g_width = width;
	g_height = height;
}

void SimulationLoop()
{
	// update the simulation
	glDepthMask(true);
	glClearColor(0.2,0.2,0.2,0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	static double s = GetSeconds();
	double delta = GetSeconds()-s;
	s = GetSeconds();
	
	g_sim->Step();

	// draw output
	g_view->Apply();
	static float r=0.0f;
	r += delta*10.0f;

	glRotatef(r, 0.0f, 1.0f, 0.0f);
	g_sim->DrawDebug(g_drawmode, g_drawdepth);

	char buf[255];
	sprintf(buf, "Frame: %.2fms", delta*1000.0f);
	DrawText(10, 20, buf);

	// flip
	glutSwapBuffers();
}

void DrawText(int x, int y, const char *string)
{
	int len, i;

	glDisable(GL_TEXTURE_2D);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	int w = glutGet(GLUT_WINDOW_WIDTH);
	int h = glutGet(GLUT_WINDOW_HEIGHT);
	gluOrtho2D(0, w, h, 0);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glColor3f(1.0f, 1.0f, 1.0f);
	glRasterPos2i(x, y);
	len = (int) strlen(string);
	for (i = 0; i < len; i++)
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[i]);

	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}

void ArrowKeys(int key, int /*x*/, int /*y*/)
{
	//if (key == GLUT_KEY_A)
	{
	}
}

void ArrowKeysUp(int key, int /*x*/, int /*y*/)
{
}

void Keyboard(unsigned char key, int /*x*/, int /*y*/)
{
	switch (key)
	{
	case '1':
		g_drawmode = IlluminationSim::kDrawSurfels;
		break;
	case '2':
		g_drawmode = IlluminationSim::kDrawIrradiance;
		break;
	case '3':
		g_drawmode = IlluminationSim::kDrawRadiance;
		break;
	case '4':
		g_drawmode = IlluminationSim::kDrawEmissive;
		break;
	case '5':
		g_drawmode = IlluminationSim::kDrawMesh;
		break;
	case 'u':
		g_drawdepth++;
		break;
	case 'i':
		g_drawdepth--;
		break;
	case ' ':
		g_sim->Step();
		break;
	case 27:
		exit(0);
		break;
 	}
}

static int lastx;
static int lasty;
static int button;

void MouseFunc(int b, int state, int x, int y)
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
			button = b;
			lastx = x;
			lasty = y;
		}
	}
}

void MotionFunc(int x, int y)
{

	lastx = x;
	lasty = y;
}

int main(int argc, char* argv[])
{
	// init gl
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH );

	// fullscreen
	//glutGameModeString( "1280x720:32@75" );
	//glutEnterGameMode();

	glutInitWindowSize(g_width, g_height);
	glutCreateWindow("Indirect Illumination");

	glewInit();
	if (glewIsSupported("GL_VERSION_2_0"))
		printf("Ready for OpenGL 2.0\n");
	else {
		printf("OpenGL 2.0 not supported\n");
		exit(1);
	}

	Init(argc, argv);

	glutMouseFunc(MouseFunc);
	glutReshapeFunc(Reshape);
	glutDisplayFunc(SimulationLoop);
	glutKeyboardFunc(Keyboard);
	glutIdleFunc(SimulationLoop);	
	glutSpecialFunc(ArrowKeys);
	glutSpecialUpFunc(ArrowKeysUp);
	glutMotionFunc(MotionFunc);
	//glutPassiveMotionFunc(MotionFunc);

	glutMainLoop();

	return 0;
}







