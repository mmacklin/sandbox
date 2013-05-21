// Main.cpp : Defines the entry point for the console application.
//
#include <vector>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>

//#include "glut.h"
#include "stdlib.h" // for atexit()

#include "kernel.h"
#include "maths.h"
#include "scene.h"
#include "particlecontainer.h"
#include "particleemitter.h"
#include "camera.h"
#include "view.h"
#include "texture.h"

void FluidInit();
void FluidTick(float);
void FluidDraw();

using namespace std;

int g_width = 512;
int g_height = 512;

Scene* g_scene;
ParticleContainer* g_particleContainer;
ParticleEmitter* g_particleEmitter;
Camera* g_camera;
View* g_view;
Texture* g_metaTex;

float g_alphaRef = 0.75f;
bool g_enable = false;

void DrawText(int x, int y, const char *string);
void DrawBackground();

void DeInit()
{
	Kernel::Destroy();
}

void Init(int argc, char* argv[])
{
	atexit(DeInit);

	Kernel::Create();

	// runs init script and such
	Kernel::Get().ParseCommandLine(argc, argv);
	Kernel::Get().PostCreate();

	g_scene = new Scene();

	Kernel::Get().SetScene(g_scene);

	// tracking camera
	g_camera = new Camera();
	g_camera->SetPosition(Vec3(0, -1, 10));

	// main scene view
	g_view = new View();
	g_view->SetCamera(g_camera);
	g_view->m_x = 0;
	g_view->m_y = 0;
	g_view->m_width = g_width;
	g_view->m_height = g_height;

	g_particleContainer = new ParticleContainer();
	g_particleContainer->m_additive = true;
	g_particleContainer->m_acceleration = Vec4(0.0f, -9.8f, 0.0, 0.0f);
	g_particleContainer->m_collisionEnabled = true;
	g_particleContainer->m_collisionPlane = Plane(Vec4(0.0f, 1.0f, 0.0f, 2.0f));
	g_particleContainer->m_scale.LoadFromString("3, 0, 1, 0.8, 1, 1, 0");
	//g_particleContainer->m_billboardMode = ParticleContainer::kBillboardStretch;
	g_particleContainer->m_stretchAmount = 0.1f;

	g_particleEmitter = new ParticleEmitter();
	g_particleEmitter->SetContainer(g_particleContainer);
	g_particleEmitter->SetRotation(Rotation(0.0f, DegToRad(45.0f), 0.0f));
	g_particleEmitter->m_spawnSize = 0.5f;
	g_particleEmitter->m_spawnLifetime = 10.0f;
	g_particleEmitter->m_emitOuterRadius = 10.0f;
	g_particleEmitter->m_spawnSpeed = 6.0f;
	g_particleEmitter->m_emitRate = 50;
	
	//g_particleEmitter->m_

	g_scene->AddActor(g_particleContainer);
	g_scene->AddActor(g_particleEmitter);
	g_scene->AddActor(g_camera);
	g_scene->AddView(g_view);

	g_metaTex = new Texture(g_width, g_height);

	// FluidInit
	//FluidInit();
}

void Reshape(int width, int height)
{
	Kernel::Get().ResizeWindow(width, height);

	// recreate meta texture
	delete g_metaTex;
	g_metaTex = new Texture(width, height);

	g_view->m_width = width;
	g_view->m_height = height;

	g_width = width;
	g_height = height;
}

void SimulationLoop()
{
	Kernel::Get().Tick();

//	if (g_enable)
//		FluidTick(1/60.0f);//Kernel::Get().GetFrameDelta());

}

void OnPostRender()
{
	if (!g_enable)
		return;

	g_metaTex->BindToSampler(0);

	// copy the frame buffer into our texture then blit to the screen
	glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, g_width, g_height);
	

	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glDisable(GL_CULL_FACE);
	glDepthMask(false);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
	//glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER, g_alphaRef);

	glDisable(GL_STENCIL_TEST);

	g_metaTex->BindToSampler(0);

	glBegin(GL_QUADS);
		
	glColor4f(0.2,0.5,1,1);
	glTexCoord2f(0,1);
	glVertex4f(-1, 1, 0,1);
	
	glTexCoord2f(1,1);
	glVertex4f(1, 1, 0,1);
	
	glTexCoord2f(1,0);
	glVertex4f(1, -1, 0,1);
	
	glTexCoord2f(0,0);
	glVertex4f(-1, -1, 0,1);

	glEnd();

	glDisable(GL_ALPHA_TEST);

	static double f = 0.0f;
	double dt = GetSeconds() - f;
	f = GetSeconds();

	char buf[255];
	sprintf(buf, "Frame time: %.2f ms", dt*1000.0f);
	DrawText(10, 20, buf);

	// controls
	sprintf(buf, "Edge threshold %.2f (a/z to change)", g_alphaRef);
	DrawText(10, 50, buf);

	sprintf(buf, "Emit speed %.2f (s/x to change)", g_particleEmitter->m_spawnSpeed);
	DrawText(10, 70, buf);

	//FluidDraw();
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
	case 'a':
		g_alphaRef += 0.01f;
		break;
	case 'z':
		g_alphaRef -= 0.01f;
		break;
	case 's':
		g_particleEmitter->m_spawnSpeed += 0.1f;
		break;
	case 'x':
		g_particleEmitter->m_spawnSpeed -= 0.1f;
		break;
	case 'd':
		g_particleEmitter->Rotate(Rotation(0.0f, DegToRad(1.0f), 0.0f));
		break;
	case 'c':
		g_particleEmitter->Rotate(Rotation(0.0f, DegToRad(-1.0f), 0.0f));
		break;
	case ' ':		
		FluidTick(1/60.0f);
		break;
	case 'u':
		FluidInit();
		break;
	case 'p':
		g_enable = !g_enable;
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

	if (button = GLUT_LEFT_BUTTON)
		g_particleEmitter->Rotate(Rotation(0.0f, DegToRad(y-lasty), 0.0f));
	else
		g_particleEmitter->Translate(Vec3(x-lastx, y-lasty, 0.0f)*0.05f);

	lastx = x;
	lasty = y;
}

void JoystickFunc(int x, int y, int z, unsigned long buttons)
{
}

void SqRegisterGameClasses()
{

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
	glutCreateWindow("Metaballs");

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







