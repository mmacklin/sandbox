#include <core/types.h>
#include <core/maths.h>
#include <core/platform.h>
#include <core/hashgrid.h>
#include <core/shader.h>
#include <core/tga.h>

#include "solve.h"

#include <iostream>

using namespace std;

const uint32_t kWidth = 800;
const uint32_t kHeight = 600;
const float kWorldSize = 2.0f;
const float kZoom = kWorldSize*2.5f;

int kNumParticles = 0;
const int kNumIterations = 5;

const float kDt = 1.0f/60.0f;
const float kRadius = 0.05f;

GrainSystem* g_grains;
GrainParams g_params;

vector<Vec2> g_positions;
vector<Vec2> g_velocities;
vector<float> g_radii;

vector<uint32_t> g_springIndices;
vector<float> g_springLengths;

bool g_pause = true;
bool g_step = false;

uint32_t g_scene = 1;

// mouse
static int lastx;
static int lasty;

Vec2 ScreenToScene(int x, int y)
{
	float aspect = float(kWidth)/kHeight;

	float left = -kZoom*aspect;
	float right =  kZoom*aspect;
	float bottom = -0.5;
	float top = 2*kZoom-0.5f;
	
	float tx = x / float(kWidth);
	float ty = 1.0f - (y / float(kHeight));

	return Vec2(left + tx*(right-left), bottom + ty*(top-bottom));
}

void Init(int scene)
{	
	g_positions.resize(0);
	g_velocities.resize(0);
	g_radii.resize(0);
	g_springIndices.resize(0);
	g_springLengths.resize(0);
		
	g_params.mGravity = Vec2(0.0f, -9.8f);
	g_params.mDamp = 0.0f;//powf(0.5f, float(kNumIterations));
	g_params.mBaumgarte = 0.2f;
	g_params.mFriction = 0.8f;
	g_params.mRestitution = 0.1f;
	g_params.mOverlap = kRadius*0.1f;
	g_params.mPlanes[2] = Vec3(1.0f, 0.0f, -5.0f);
	g_params.mPlanes[1] = Vec3(-1.0f, 0.0f, -5.0f);
	g_params.mPlanes[0] = Normalize(Vec3(0.0f, 1.0f, 0.0f));
	g_params.mNumPlanes = 3;

	switch (scene)
	{
	case 1:
	{
		for (int x=0; x < 64; ++x)
		{
			float s = -4.5f;

			const float sep = 1.0f*kRadius;

			for (int i=0; i < 32; ++i)
			{
				s += 2.0f*sep;// + Randf(0.0f, 0.05f)*kRadius;

				g_positions.push_back(Vec2(s, 1.0f + sep + 2.0f*x*sep));
				g_velocities.push_back(Vec2());
				g_radii.push_back(kRadius);// + kRadius*Randf(-0.1f, 0.0f));
			}
		}
		break;
	}	
	case 2:		
	case 3:
	case 4:
	case 5:
	{
		const char* file;
		if (scene == 2 || scene == 3)
			file = "bunny.tga";
		else if (scene == 4 || scene == 5)
			file = "armadillo.tga";

		TgaImage img;
		if (TgaLoad(file, img))
		{
			float xstart = -3.0f;

			float step = kRadius*1.0f;
			float x = xstart;
			float y = 1.5f;

			int dim = 80; 

			float dpx = float(img.m_width) / dim;
			float dpy = float(img.m_height) / dim;

			vector<uint32_t> lookup(dim*dim, -1);	

			for (int i=0; i < dim; ++i)
			{
				int py = i*dpy;

				for (int j=0; j < dim; ++j)
				{
					int px = j*dpx;

					uint32_t c = img.SampleClamp(px, py);
					
					if (c != 0)
					{
						uint32_t newIndex = g_positions.size(); 
						lookup[i*dim + j] = newIndex;

						float r = Randf(0.0f, 0.009f)*step;
						g_positions.push_back(Vec2(x + r , y));
						g_velocities.push_back(0.0f);//Vec2(10.0f, 0.0f));
						g_radii.push_back(kRadius);// + kRadius*Randf(-0.2f, 0.0f));

						// add springs
						if (scene & 1)
							{
							for (int ny=i-1; ny <= i+1; ++ny)
							{
								for (int nx=j-1; nx <= j+1; ++nx)
								{
									uint32_t r = lookup[ny*dim + nx];

									if (r != uint32_t(-1) && r != newIndex)
									{	
										g_springIndices.push_back(newIndex);
										g_springIndices.push_back(r);

										g_springLengths.push_back(Distance(g_positions[newIndex], g_positions[r]));
									}
								}
							}
						}
					}

					x += 2.0f*step;
				}

				x = xstart;

				y += 2.0f*step; 
			}	
		}
		break;
	}
	case 6:
	{
		g_positions.push_back(Vec2(0.0f, kRadius));
		g_velocities.push_back(Vec2(0.0f, 0.0f));
		g_radii.push_back(kRadius);// + kRadius*Randf(-0.1f, 0.0f));
		
			
		g_positions.push_back(Vec2(kRadius, kRadius + 2.0f*kRadius));
		g_velocities.push_back(Vec2(0.0f, 0.0f));
		g_radii.push_back(kRadius);// + kRadius*Randf(-0.1f, 0.0f));

		break;
	}
	case 7:
	{
		g_positions.push_back(Vec2(-0.2f, 1.0f));
		g_velocities.push_back(Vec2(1.0f, 0.0f));
		g_radii.push_back(kRadius);// + kRadius*Randf(-0.1f, 0.0f));
		
			
		g_positions.push_back(Vec2(0.2f, 1.0f));
		g_velocities.push_back(Vec2(-1.0f, 0.0f));
		g_radii.push_back(kRadius);// + kRadius*Randf(-0.1f, 0.0f));

		break;
	}	
	case 8:
	{
		g_params.mPlanes[0] = Normalize(Vec3(1.0f, 1.0f, 0.0f));

		g_positions.push_back(Vec2(0.0f, 1.0f));
		g_velocities.push_back(Vec2(0.0f, 0.0f));
		g_radii.push_back(kRadius);// + kRadius*Randf(-0.1f, 0.0f));

		break;
	}	
	case 9:
	{
		// pyramid
		const int kLevels = 10;

		for (int y=0; y < kLevels; ++y)
		{
			for (int x=0; x < kLevels-y; ++x)
			{
				g_positions.push_back(0.5f*Vec2(0.0f + 1.0f*y*kRadius + x*2.0f*kRadius, 0.5f + kRadius + 1.6f*y*kRadius));
				g_velocities.push_back(Vec2());
				g_radii.push_back(kRadius);// + kRadius*Randf(-0.1f, 0.0f));
			}
		}

		break;
	}
	default: 
		break;
	};

	kNumParticles = g_positions.size();

	g_grains = grainCreateSystem(kNumParticles);

	grainSetParams(g_grains, &g_params);
	grainSetPositions(g_grains, (float*)&g_positions[0], kNumParticles);
	grainSetVelocities(g_grains, (float*)&g_velocities[0], kNumParticles);
	grainSetRadii(g_grains, &g_radii[0]);

	if (!g_springIndices.empty())
		grainSetSprings(g_grains, &g_springIndices[0], &g_springLengths[0], g_springLengths.size());

}

void Shutdown()
{
	grainDestroySystem(g_grains);
}

void Reset()
{
	Shutdown();
	Init(g_scene);
}

void DrawCircle(const Vec2& p, float r, const Colour& c )
{
	glBegin(GL_TRIANGLE_FAN);
	glColor3f(c.r, c.g, c.b);
	glVertex2fv(p);
	
	const int kSegments = 40;
	for (int i=0; i < kSegments+1; ++i)
	{
		float theta = k2Pi*float(i)/kSegments;
		
		float y = p.y + r*Cos(theta);
		float x = p.x + r*Sin(theta);
		
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

void GLUTUpdate()
{
	//---------------------------

	glViewport(0, 0, kWidth, kHeight);
	glDisable(GL_CULL_FACE);
	glDisable(GL_DEPTH_TEST);
	
	glPointSize(5.0f);

	float aspect = float(kWidth)/kHeight;
	float viewWidth = kZoom*aspect;
	
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(OrthographicMatrix(-viewWidth, viewWidth, -0.5, 2*kZoom-0.5f, 0.0f, 1.0f));

	glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	double drawStart = GetSeconds();

	// Step
	GrainTimers timers;

	if (!g_pause || g_step)
	{
		glClear(GL_COLOR_BUFFER_BIT);

		grainSetParams(g_grains, &g_params);
		grainUpdateSystem(g_grains, kDt, kNumIterations, &timers);

		g_step = false;
	}

	double drawEnd = GetSeconds();

	for (int i=0; i < g_params.mNumPlanes; ++i)
	{	
		Vec2 p = g_params.mPlanes[i].z * Vec2(g_params.mPlanes[i].x, g_params.mPlanes[i].y);
		Vec2 d = Vec2(-g_params.mPlanes[i].y, g_params.mPlanes[i].x);
		
		glBegin(GL_LINES);
		glColor3f(1.0f, 1.0f, 1.0f);
		glVertex2fv(p - d*1000.0);
		glVertex2fv(p + d*1000.0);
		glEnd();		
	}
		
	// read-back data	
	grainGetPositions(g_grains, (float*)&g_positions[0]);
	grainGetVelocities(g_grains, (float*)&g_velocities[0]);
	
	glColor3f(0.7f, 0.7f, 0.8f);

	Colour colors[] = { Colour(0.5f, 0.5f, 1.0f),
						Colour(1.0f, 0.5f, 0.5f),
						Colour(0.5f, 1.0f, 0.5f) };

	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

	glPointSize(0.5f*kRadius*kWidth/viewWidth + 1.0f);

	glBegin(GL_POINTS);

	for (int i=0; i < kNumParticles; ++i)
	{
		glColor3fv(colors[i%3]);
		glVertex2fv(g_positions[i]);
	}

	glEnd();


	glDisable(GL_BLEND);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, kWidth, kHeight, 0);

	int x = 10;
	int y = 15;
	
	glColor3f(1.0f, 1.0f, 1.0f);
	DrawString(x, y, "Sim: %.2f", (drawEnd-drawStart)*1000.0f); y += 13;
	DrawString(x, y, "1-3: Select scene", ""); y += 13;
	DrawString(x, y, "r: Reset", ""); y += 13;
	DrawString(x, y, "p: Pause", ""); y += 13;
	DrawString(x, y, "o: Step", ""); y += 13;

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
	if (key > '0' && key <= '9')
	{
		g_scene = key-'0';
		Init(g_scene);
		return;
	}
	
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
			g_params.mNumPlanes--;
			break;
		}
		case 'o':
		{
			g_step = true;
			break;
		}
		case '=':
		{
			g_params.mNumPlanes++;
			break;
		}
		case '-':
		{
			g_params.mNumPlanes--;
			break;
		}
		case 'p':
		{
			g_pause = !g_pause;
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
	}
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

	
	g_positions[0] = ScreenToScene(x, y);
	g_velocities[0] = Vec2(dx*5.0f, 5.0f*(1.0f-dy));

	grainSetPositions(g_grains, (float*)&g_positions[0], 1);
	grainSetVelocities(g_grains, (float*)&g_velocities[0], 1);
	
}

int solveCuda(float* a, float* b, float* c, int n);
#include <xmmintrin.h>

int main(int argc, char* argv[])
{
	RandInit();
	Init(g_scene);
	
    // init gl
    glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
	
    glutInitWindowSize(kWidth, kHeight);
    glutCreateWindow("Granular2d");
    glutPositionWindow(200, 100);
		
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

