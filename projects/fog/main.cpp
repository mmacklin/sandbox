// FogVolumes.cpp : 
//

#include "core/perlin.h"
#include "core/shader.h"
#include "core/types.h"
#include "core/maths.h"
#include "core/platform.h"

#include <vector>

using namespace std;

uint32_t g_width = 1280;
uint32_t g_height = 720;

GLuint g_staticPrograms[2];
GLuint g_depthBuffer;

enum LightType 
{
	kLightSpot=0,
	kLightPoint
};

struct Light
{	
	LightType m_type;
	
	Matrix44 m_lightToWorld;
	Vector4 m_lightColour;	
};

Light g_lights[2];

int g_selectedLightIndex = 0;
int g_activeLights = 2;
float g_scatteringCoefficient = 0.3f;
bool g_showHelp =  true;

Light& GetSelectedLight() { return g_lights[g_selectedLightIndex]; }

void Init()
{ 
	g_staticPrograms[kLightSpot] = CompileProgramFromFile("Data/SimpleVertex.glsl", "Data/SpotLight.glsl");
	g_staticPrograms[kLightPoint]= CompileProgramFromFile("Data/SimpleVertex.glsl", "Data/PointLight.glsl");

	// create some default lights
	g_lights[0].m_lightToWorld = TransformMatrix(Rotation(), Point3(5.0f, 10.0f, 0.0f));
	g_lights[0].m_lightColour = Vector4(0.3f, 0.3f, 0.8f, 1.0f)*4.0f;
	g_lights[0].m_type = kLightSpot ;

	g_lights[1].m_lightToWorld = TransformMatrix(Rotation(), Point3(-5.0f, 10.0f, 0.0f));
	g_lights[1].m_lightColour = Vector4(0.8f, 0.3f, 0.3f, 1.0f)*4.0f;
	g_lights[1].m_type = kLightPoint;

	RandInit();
}

void Shutdown()
{

}

void RenderStaticGeometry()
{	
	glDisable(GL_TEXTURE_2D);

	glBegin(GL_QUADS);
	
	glNormal3f(0.0f, 1.0f, 0.0f);
	glVertex4f(-200.0f, -2.0f, 200.0f, 1.0f);
	glVertex4f(200.0f, -2.0f, 200.0f, 1.0f);
	glVertex4f(200.0f, -2.0f, -200.0f, 1.0f);
	glVertex4f(-200.0f, -2.0f, -200.0f, 1.0f);

	glNormal3f(0.0f, 0.0f, 1.0f);
	glVertex4f(-200.0f, -200.0f, -2.0f, 1.0f);
	glVertex4f(200.0f, -200.0f, -2.0f, 1.0f);
	glVertex4f(200.0f, 200.0f, -2.0f, 1.0f);
	glVertex4f(-200.0f, 200.0f, -2.0f, 1.0f);
	
	glEnd();
}


void SetLightParams(Light& l, GLuint program)
{
	// disable depth test for full screen quad
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glDepthMask(true);
	glBlendFunc(GL_ONE, GL_ONE);
	glEnable(GL_CULL_FACE);
	glDisable(GL_LIGHTING);
	glDisable(GL_ALPHA_TEST);
	glDepthFunc(GL_LEQUAL);

	glUseProgram(program);
	
	GLuint param = glGetUniformLocation(program, "g_lightToWorld");
	glUniformMatrix4fv(param, 1, false, l.m_lightToWorld);

	param = glGetUniformLocation(program, "g_worldToLight");
	glUniformMatrix4fv(param, 1, false, AffineInverse(l.m_lightToWorld));

	param = glGetUniformLocation(program, "g_lightCol");
	glUniform4fv(param, 1, l.m_lightColour);

	param = glGetUniformLocation(program, "g_scatteringCoefficient");
	glUniform1f(param, g_scatteringCoefficient);

	/*
	param = glGetUniformLocation(program, "g_noiseTexture");
	if (param != -1)
	{
		glVerify(glActiveTexture((GLenum)GL_TEXTURE1));	
		glVerify(glEnable(GL_TEXTURE_3D));
		glVerify(glBindTexture(GL_TEXTURE_3D, g_noiseTexture));
	}
	*/
	
}


void RenderHelp()
{	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0f, g_width, g_height, 0.0f);

	glUseProgram(0);

	int x = 50;
	int y = 50;
	int dy = 12;

	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

	DrawString(x, y, "space - switch light"); y += dy;
	DrawString(x, y, "left drag - translate"); y += dy;
	DrawString(x, y, "right drag - rotate"); y += dy;
	DrawString(x, y, "t - change light type"); y += dy;
	DrawString(x, y, "u,j - scattering coefficient (%.3f)", g_scatteringCoefficient); y += dy;
	DrawString(x, y, "i,k - light intensity (%.2f, %.2f, %.2f)", GetSelectedLight().m_lightColour.x, GetSelectedLight().m_lightColour.y, GetSelectedLight().m_lightColour.z); y += dy;
	DrawString(x, y, "h - toggle help"); y += dy;

}

void Render()
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0f, 20.0f, 25.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0f, float(g_width) / g_height, 0.1f, 1000.0f);
	
	// z-pre pass
	glColorMask(false, false, false, false);
	
	RenderStaticGeometry();

	glColorMask(true, true, true, true);

	for (int i=0; i < g_activeLights; ++i)
	{
		Light& l = g_lights[i];

		
		// select shader program and render
		GLuint staticProgram = g_staticPrograms[l.m_type];
		SetLightParams(l, staticProgram);

		RenderStaticGeometry();
	}

	// debug
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
	
	if (g_showHelp)
		RenderHelp();
}



//--------------------------------------------------------------------------
// Glut callbacks

void GlutUpdate()
{
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	static double s_lastUpdate = GetSeconds();
	
	double curTime = GetSeconds();
	float dt = float(curTime - s_lastUpdate);
	s_lastUpdate = curTime;

	Render();

	// flip
	glutSwapBuffers();
}

void GlutReshape(int width, int height)
{
	g_width = width;
	g_height = height;

	glViewport(0, 0, width, height);
}

void GlutArrowKeys(int key, int x, int y)
{
}

void GlutArrowKeysUp(int key, int x, int y)
{
}

void GlutKeyboardDown(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'u':
		{
			g_scatteringCoefficient += 0.001f;
			break;
		}
	case 'j':
		{
			g_scatteringCoefficient = max(0.0f, g_scatteringCoefficient-0.001f);
			break;
		}
	case 't':
		{
			// swap light type
			GetSelectedLight().m_type = (GetSelectedLight().m_type == kLightPoint) ? kLightSpot : kLightPoint;
			break;
		}
	case ' ':
		{
			g_selectedLightIndex = (g_selectedLightIndex+1) % g_activeLights;
			break;
		}
	case 'i':
		{
			GetSelectedLight().m_lightColour *= 1.1f;
			break;
		}
	case 'k':
		{
			GetSelectedLight().m_lightColour *= 0.9f;
			break;
		}
	case 'h':
		{
			g_showHelp = !g_showHelp;
			break;
		}
	case 27:
		exit(0);
		break;
	};

}

void GlutKeyboardUp(unsigned char key, int x, int y)
{
}

static int lastx;
static int lasty;
static int button;

void GlutMouseFunc(int b, int state, int x, int y)
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
			button = b;
		}
	}
}

void GlutMotionFunc(int x, int y)
{
	int dx = x - lastx;
	int dy = y - lasty;

	lastx = x;
	lasty = y;
	
	Light& selectedLight = g_lights[g_selectedLightIndex];
	
	if (button == GLUT_LEFT_BUTTON)
	{
		Point3 t = selectedLight.m_lightToWorld.GetTranslation();
		
		t.x += dx * 0.01f;
		t.z += dy * 0.01f;

		selectedLight.m_lightToWorld.SetTranslation(t);
	}
	else
	{
		selectedLight.m_lightToWorld *= RotationMatrix(dx * 0.01f, Vec3(1.0f, 0.0f, 0.0f));
		selectedLight.m_lightToWorld *= RotationMatrix(dy * 0.01f, Vec3(0.0f, 0.0f, 1.0f));
	}
}

int main(int argc, char* argv[])
{
	// init glc
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH );

	glutInitWindowSize(g_width, g_height);
	glutCreateWindow("Fog Volumes");
	glutPositionWindow(200, 200);

	glewInit();
	if (glewIsSupported("GL_VERSION_2_0"))
	{
		Log::Info << "Ready for OpenGL 2.0" << endl;
	}
	else 
	{
		Log::Warn << "OpenGL 2.0 not supported" << endl;
		return 1;
	}

	float one = 1.0f;

	Init();

	glutMouseFunc(GlutMouseFunc);
	glutReshapeFunc(GlutReshape);
	glutDisplayFunc(GlutUpdate);
	glutKeyboardFunc(GlutKeyboardDown);
	glutKeyboardUpFunc(GlutKeyboardUp);
	glutIdleFunc(GlutUpdate);	
	glutSpecialFunc(GlutArrowKeys);
	glutSpecialUpFunc(GlutArrowKeysUp);
	glutMotionFunc(GlutMotionFunc);

	glutMainLoop();

	Shutdown();

	return 0;	
}

