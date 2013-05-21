#include <core/types.h>
#include <core/maths.h>
#include <core/platform.h>
#include <core/hashgrid.h>
#include <core/maths.h>
#include <core/shader.h>
#include <core/mesh.h>

#include <iostream>

using namespace std;

#include "solve.h"
#include "shaders.h"

const uint32_t kWidth = 1024;
const uint32_t kHeight = 512;

int kNumParticles = 0;
const int kNumIterations = 5;
const float kRadius = 0.05f;

GrainSystem* g_grains;
GrainParams g_params;

vector<Vec3> g_positions;
vector<Vec3> g_velocities;
vector<float> g_radii;
vector<uint32_t> g_springIndices;
vector<float> g_springLengths;

Vec3 g_camPos(1.0f, 2.0f, 7.0f);
Vec3 g_camAngle(0.0f, 0.0f, 0.0f);

bool g_pause = false;
bool g_step = false;

// mouse
static int lastx;
static int lasty;

// render funcs
GLuint g_shadowTex;
GLuint g_shadowBuf; 

// size of initial grid of particles
const int kParticleHeight = 16;
const int kParticleWidth = 8;
const int kParticleLength = 16;

void Voxelize(const Mesh& mesh, uint32_t width, uint32_t height, uint32_t depth, uint32_t* volume);

int g_scene = 1;

void Init(int scene)
{	
	g_positions.resize(0);
	g_velocities.resize(0);
	g_radii.resize(0);
	g_springIndices.resize(0);
	g_springLengths.resize(0);

	if (g_grains)
		grainDestroySystem(g_grains);

	// sim params	
	g_params.mGravity = Vec3(0.0f, -9.8f, 0.0f);
	g_params.mDamp = 0.0f;//powf(0.3f, float(kNumIterations));
	g_params.mBaumgarte = 0.2f;
	g_params.mFriction = 0.5f;
	g_params.mRestitution = 0.1f;
	g_params.mOverlap = kRadius*0.05f;
	g_params.mPlanes[0] = Vec4(0.0f, 1.0f, 0.0f, 0.0f);
	g_params.mPlanes[1] = Vec4(1.0f, 0.0f, 0.0f, 2.0f);
	g_params.mNumPlanes = 2;
	g_params.mDissipation = 0.2f;
	
	float r = kRadius*0.95f;
	
	switch (scene)
	{
	case 1:
	{	
		g_params.mDissipation = 0.5f;

		float y = 0.0f + kRadius;

		for (int i=0; i < kParticleHeight; ++i)
		{	
			for (int z=0; z < kParticleLength; ++z)
			{
				for (int x=0; x < kParticleWidth; ++x)
				{
					g_positions.push_back(Vec3(x*2.0f*r + Randf(0.0f, 0.05f*r), y, z*2.0f*r + Randf(0.0f, 0.05f*r)));				
					g_velocities.push_back(Vec3());
					g_radii.push_back(kRadius);
				}
			}

			y += 2.f*r;
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
			file = "armadillo.ply";
		else if (scene == 4 || scene == 5)
			file = "bunny.ply";

		// voxelize mesh
		Mesh* m = ImportMeshFromPly(file);
	
		if (m)
		{
			// square cells
			Vec3 minExtents, maxExtents;
			m->GetBounds(minExtents, maxExtents);
			
			Vec3 extents(maxExtents-minExtents);

			// normalize scale
			float longestAxis = max(max(extents.x, extents.y), extents.z);
			float scale = 1.5f;			
			
			if (scene == 2 || scene == 3)
				m->Transform(RotationMatrix(kPi, Vec3(0.0f, 1.0f, 0.0f)));

			m->Transform(ScaleMatrix(scale/longestAxis));
			extents *= scale/longestAxis;

			Vec3 dim(extents/kRadius);

			uint32_t dx(dim.x);
			uint32_t dy(dim.y);
			uint32_t dz(dim.z);

			printf("%d %d %d\n", dx, dy, dz);
			printf("begin voxelize\n");

			uint32_t* volume = new uint32_t[dx*dy*dz];
			Voxelize(*m, dx, dy, dz, volume);

			printf("end voxelize\n");

			float yoff = 0.0f + r;

			if (scene & 1)
				yoff += 1.0f;

			// create a particle at each non-empty cell
			for (uint32_t x=0; x < dx; ++x)
			{
				for (uint32_t y=0; y < dy; ++y)
				{
					for (uint32_t z=0; z < dz; ++z)
					{
						uint32_t index = z*dx*dy + y*dx + x;

						if (volume[index])
						{
							g_positions.push_back(Vec3(x*r*2.01f, yoff + y*r*2.01f, z*r*2.01f));
							g_velocities.push_back(0.0f);//Vec3(-10.0f, 0.0f, 0.0f));
							g_radii.push_back(kRadius);						

							volume[index] = g_positions.size();

							if (scene & 1)
							{
								g_velocities.back() = Vec3(-10.0f, 0.0f, 0.0f);
								// add springs
								const int stride = 1;

								for (uint32_t i=max(0, int(x)-stride); i <= x; i+=stride)
								{
									for (uint32_t j=max(0, int(y)-stride); j <= y; j+=stride)
									{
										for (uint32_t k=max(0, int(z)-stride); k <= z; k+=stride)
										{
											uint32_t r = volume[k*dx*dy + j*dx + i];

											if (r && r != uint32_t(-1) && r != volume[index])
											{	
												g_springIndices.push_back(volume[index]-1);
												g_springIndices.push_back(r-1);

												g_springLengths.push_back(Distance(g_positions.back(), g_positions[r-1]));
											}
										}
									}
								}
							}
							else
							{
								g_params.mDissipation = 0.4f;
								g_positions.back() += Vec3(Randf(-0.001f, 0.001f), Randf(-0.001f, 0.001f), Randf(-0.001f, 0.001f));
							}
						}
					}
				}
			}

			delete m;
			delete[] volume;
		}
		else
			printf("Couldn't open %s for read\n", file); 

		break;
	}
	default:
		break;
	}

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


void GLUTUpdate()
{
	GrainTimers timers;

	if (!g_pause || g_step)
	{
		grainSetParams(g_grains, &g_params);
		grainUpdateSystem(g_grains, 1.0f/60.0f, kNumIterations, &timers);

		g_step = false;
	}

	//---------------------------

	glClearColor(0.492f, 0.684f, 0.999f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	
	glPointSize(5.0f);

	float aspect = float(kWidth)/kHeight;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0f, aspect, 0.1f, 1000.0f);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(RadToDeg(-g_camAngle.x), 0.0f, 1.0f, 0.0f);
	glRotatef(RadToDeg(-g_camAngle.y), cosf(-g_camAngle.x), 0.0f, sinf(-g_camAngle.x));	
	glTranslatef(-g_camPos.x, -g_camPos.y, -g_camPos.z);


	// read-back data	
	grainGetPositions(g_grains, (float*)&g_positions[0]);
	grainGetRadii(g_grains, (float*)&g_radii[0]);

	Vec3 grainsLower(FLT_MAX, FLT_MAX, FLT_MAX);
	Vec3 grainsUpper(-grainsLower);
	
	for (int i=0; i < kNumParticles; ++i)
	{
		grainsLower = Min(g_positions[i], grainsLower);
		grainsUpper = Max(g_positions[i], grainsUpper);
	}

	Vec3 grainsExtents = grainsUpper - grainsLower;
	Vec3 grainsCenter = 0.5f*(grainsUpper + grainsLower);

	// shadowing pass 
	//
	Vec3 lightDir = Normalize(Vec3(5.0f, 20.0f, 20.0f));
	Vec3 lightPos = grainsCenter + lightDir*Length(grainsExtents)*5.0f;
	Vec3 lightTarget = grainsCenter;

	// calculate tight bounds for shadow frustum
	float lightFov = 2.0f*atanf(Length(grainsUpper-grainsCenter)/Length(lightPos-grainsCenter));
	Matrix44 lightPerspective = ProjectionMatrix(RadToDeg(lightFov), 1.0f, 1.0f, 1000.0f);
	Matrix44 lightView = LookAtMatrix(Point3(lightPos), Point3(lightTarget));
	Matrix44 lightTransform = lightPerspective*lightView;

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadMatrixf(lightPerspective);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixf(lightView);

	ShadowBegin(g_shadowTex, g_shadowBuf);
	DrawPoints(&g_positions[0].x, kNumParticles, kRadius, 1024, 1.0f, lightFov, lightPos, lightTarget, lightTransform, g_shadowTex);
	ShadowEnd();

	glViewport(0, 0, kWidth, kHeight);
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	// lighting pass
	glColor3f(0.7f, 0.7f, 0.8f);

	double drawStart = GetSeconds();

	DrawPoints(&g_positions[0].x, kNumParticles, kRadius, kWidth, float(kWidth)/kHeight, kPi/4.0f, lightPos, lightTarget, lightTransform, g_shadowTex);

	glDisable(GL_LIGHTING);
	glDisable(GL_BLEND);

	// planes
	DrawPlanes(g_params.mPlanes, g_params.mNumPlanes, lightPos, lightTarget, lightTransform, g_shadowTex);
	
	double drawEnd = GetSeconds();

	glUseProgram(0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, kWidth, kHeight, 0);

	int x = 10;
	int y = 15;
	
	glColor3f(1.0f, 1.0f, 1.0f);
	DrawString(x, y, "Draw time: %.2fms", (drawEnd-drawStart)*1000.0f); y += 13;
	DrawString(x, y, "1-5: Scene Select"); y += 13;
	DrawString(x, y, "r: Reset"); y += 13;
		
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
	if (key > '0' && key <= '5')
	{
		g_scene = key-'0';
		Init(g_scene);
		return;
	}

	const float kSpeed = 0.5f;

	// update camera
	const Vec3 forward(-sinf(g_camAngle.x)*cosf(g_camAngle.y), sinf(g_camAngle.y), -cosf(g_camAngle.x)*cosf(g_camAngle.y));
	const Vec3 right(Normalize(Cross(forward, Vec3(0.0f, 1.0f, 0.0f))));
	
	Vec3 delta;

 	switch (key)
	{
		case 'w':
		{
			g_camPos += kSpeed*forward;
			break;
		}
		case 's':
		{
			g_camPos -= kSpeed*forward;
			break;
		}
		case 'a':
		{
			g_camPos -= kSpeed*right;
			break;
		}
		case 'd':
		{
			g_camPos += kSpeed*right;
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
		case 'u':
		{
			g_params.mPlanes[2].w -= 0.1f;
			break;
		}
		case ' ':
		{
			g_pause = !g_pause;
			break;
		}
		case 'o':
		{
			g_step = true;
			break;
		}
		case 'q':
		case 27:
			exit(0);
			break;
	};

	g_camPos += delta;
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

	const float kSensitivity = DegToRad(0.1f);

	g_camAngle.x -= dx*kSensitivity;
	g_camAngle.y -= dy*kSensitivity;
}

int solveCuda(float* a, float* b, float* c, int n);


int main(int argc, char* argv[])
{	
	RandInit();
	Init(g_scene);
	
    // init gl
    glutInit(&argc, argv);

#ifdef WIN32
    glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
#else
	glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH);
#endif
	
    glutInitWindowSize(kWidth, kHeight);
    glutCreateWindow("Granular");
    glutPositionWindow(200, 100);

#ifdef WIN32
	glewInit();
#endif

	ShadowCreate(g_shadowTex, g_shadowBuf);

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

