#include <iostream>

#include "core/types.h"
#include "core/shader.h"
#include "core/platform.h"
#include "core/mesh.h"
#include "core/voxelize.h"
#include "core/sdf.h"
#include "core/pfm.h"

#define STRINGIFY(A) #A

using namespace std;

uint32_t g_screenWidth = 800;
uint32_t g_screenHeight = 600;

const uint32_t g_dim = 128;

Mesh* g_mesh;
Vec3 g_lower;
Vec3 g_upper;

Vec3 g_voxelLower;
Vec3 g_voxelUpper;

const float g_scale = 20.0f;

std::vector<uint32_t> g_voxels;
std::vector<float> g_sdf;

bool g_showMesh = false;

float g_threshold = 0.0f;

void DrawMesh(const Mesh* m, Vec3 color)
{ 
	if (m)
	{
		glColor3fv(color);

		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, &m->m_positions.front());
		glEnableClientState(GL_NORMAL_ARRAY);
		glNormalPointer(GL_FLOAT, 0, &m->m_normals.front());

		glDrawElements(GL_TRIANGLES, m->GetNumFaces()*3, GL_UNSIGNED_INT, &m->m_indices.front());

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
	}
}

void Init()
{
	//const char* file = "../../data/TheArena.ply";
	//const char* file = "../../data/sphere.ply";
	//const char* file = "teapot.ply";
	//const char* file = "../pbf3d/data/claw.ply";
	//const char* file = "../../data/TheArena.obj";
	//const char* file = "armadillo.ply";
	//const char* file = "../../data/sphere.ply";
	const char* file = "../pbf3d/data/claw.ply";


	//g_mesh = ImportMeshFromPly(file);
	//g_mesh = ImportMeshFromObj("../../data/teapot.obj");
	g_mesh = ImportMeshFromPly("jotero.ply");
	//g_mesh = ImportMeshFromObj(file);
	g_mesh->GetBounds(g_lower, g_upper);
	
	
	/*
	// copy
	Mesh* m = ImportMeshFromPly(file);
	m->Transform(TranslationMatrix(Point3(0.0f, 0.0f, (g_upper.z-g_lower.z)*1.5f)));
	g_mesh->AddMesh(*m);
	g_mesh->GetBounds(g_lower, g_upper);
	// end test
	*/
	
	

	Vec3 edges = g_upper-g_lower;
	Vec3 center = (g_lower+g_upper)*0.5f;

	float maxEdge = max(max(edges.x, edges.y), edges.z);

	/*
	Matrix44 xform = RotationMatrix(sqrtf(2.0f), Vec3(0.0f, 1.0f, 0.0f))*ScaleMatrix(Vec3(g_scale/maxEdge))*TranslationMatrix(-Point3(center));
	g_mesh->Transform(xform);
	*/
	g_mesh->Transform(TranslationMatrix(Point3(-Vec3(g_lower))));
	g_mesh->Transform(RotationMatrix(-kPi*0.5f, Vec3(1.0f, 0.0f, 0.0f)));
	g_mesh->Transform(ScaleMatrix(g_scale/maxEdge));
	//g_mesh->Transform(TranslationMatrix(Point3(4.0f, -0.1f, 0.5f)));
	

	g_mesh->GetBounds(g_lower, g_upper);
		
	edges = g_upper-g_lower;
	maxEdge = max(max(edges.x, edges.y), edges.z);

	g_voxels.resize(g_dim*g_dim*g_dim);
	g_sdf.resize(g_dim*g_dim*g_dim);

	g_voxelLower = g_lower - Vec3(maxEdge*0.1f);
	g_voxelUpper = g_lower + Vec3(maxEdge*1.1f);

	//Voxelize(*g_mesh, g_dim, g_dim, g_dim, &g_voxels[0], g_lower-edges*float(2.0f/g_dim), g_upper+edges*float(2.0f/g_dim));
	Voxelize(*g_mesh, g_dim, g_dim, g_dim, &g_voxels[0], g_voxelLower, g_voxelUpper);

	MakeSDF(&g_voxels[0], g_dim, g_dim, g_dim, &g_sdf[0]);

	/*
	PfmImage sdf;
	if (PfmLoad("../pbf3d/data/claw.pfm", sdf))
		g_sdf.assign(sdf.m_data, sdf.m_data+sdf.m_width*sdf.m_height*sdf.m_depth);

	printf("Max d: %f\n", *std::max_element(g_sdf.begin(), g_sdf.end()));
	*/

	g_upper = g_lower + Vec3(maxEdge);
}


void GLUTUpdate()
{	
	glViewport(0, 0, g_screenWidth, g_screenHeight);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	const float fov = 45.0f;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fov, g_screenWidth/float(g_screenHeight), 0.01f, 1000.0f);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	Vec3 target = (g_upper+g_lower)*0.5f;
	Vec3 camera = target + Vec3(sinf(GetSeconds()*0.5f)*g_scale*2.0f, -0.f*g_scale, cosf(GetSeconds()*0.5f)*g_scale*2.0f);	

	gluLookAt(camera.x, camera.y, camera.z, target.x, target.y, target.z, 0.0f, 1.0f, 0.0f);

	if (g_showMesh)
	{
		DrawMesh(g_mesh, Vec3(1.0f));
	}
	else
	{
		Vec3 delta = (g_voxelUpper-g_voxelLower)/Vec3(float(g_dim));

		glEnable(GL_TEXTURE_2D);
		glEnable(GL_MULTISAMPLE);
		glEnable(GL_POINT_SPRITE);
		glEnable(GL_POINT_SMOOTH);
		glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glPointSize(4.0f);
		glBegin(GL_POINTS);

		glColor4f(1.0f, 1.0f, 1.0f, 0.025f);

		for (int z=0; z < g_dim; ++z)
		{
			for (int y=0; y < g_dim; ++y)
			{
				for (int x=0; x < g_dim; ++x)
				{
					
					if (g_voxels[z*g_dim*g_dim + y*g_dim + x])
					{
						Vec3 pos = g_lower + Vec3(x, y, z)*delta;
						glVertex3fv(pos);
					}
					
					/*
					float d = g_sdf[z*g_dim*g_dim + y*g_dim + x];

					if (d <= g_threshold)
					{
						Vec3 pos = g_voxelLower + Vec3(x, y, z)*delta;
						Vec3 col = Lerp(Vec3(1.0f, 0.5f, 0.5f), Vec3(0.5f, 1.0f, 0.5f), d*0.5f + 0.5f);

						glColor3fv(col);
						glVertex3fv(pos);
					}
					*/
				}
			}
		}				

		glEnd();
	}
		
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
			g_showMesh = !g_showMesh;
			break;
		}
		case '=':
		{
			g_threshold += 0.01f;
			break;
		}
		case '-':
		{
			g_threshold -= 0.01f;
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
	
    //int dx = x-lastx;
    //int dy = y-lasty;

	lastx = x;
	lasty = y;
	
}

void GLUTPassiveMotionFunc(int x, int y)
{
    //int dx = x-lastx;
    //int dy = y-lasty;
	
	lastx = x;
	lasty = y;

}


int main(int argc, char* argv[])
{	
    // init gl
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH);
	
    glutInitWindowSize(g_screenWidth, g_screenHeight);
    glutCreateWindow("Voxelize");
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

