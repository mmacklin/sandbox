// Sandbox.cpp : Defines the entry point for the console application.
//

#include "core/maths.h"
#include "core/mesh.h"
#include "core/aabbtree.h"
#include "core/threadgroup.h"
#include "core/platform.h"

#include "external/glut/glut.h"

#include "scene.h"
#include "camera.h"
#include "integrators.h"

#include "memory.h"

float g_statLightSampleTime;
float g_statLightGenSampleTime;
float g_statSceneTraceTime;
float g_statAABBTraceTime;
float g_statAABBIntersectTime;
float g_statTriIntersectTime;
float g_statPresentTime;

using namespace std;

AABBTree* g_tree;
Colour* g_pixels;
Colour* g_filtered;
uint32_t g_iterations;

Point3 g_camPos;
Rotation g_camDir;

Scene* g_scene;

// render dimensions
uint32_t g_width = 640;
uint32_t g_height = 360;

// window dimensions
uint32_t g_windowWidth = g_width;
uint32_t g_windowHeight = g_height;

float g_exposure = 0.1f;
float g_zoom = 1.0f;
float g_sunTheta = 0.41389135;
float g_sunPhi = 0.95993042;
float g_sunTurbidity = 2.0f;

// work harder!
const uint32_t g_numWorkers = 8;

float g_flySpeed = 0.5f;

enum RenderMode
{
	eNormals,
	eComplexity,
	ePathTrace
};

RenderMode g_mode = eNormals;

struct RayJob
{
    // raster space sub rect to process
    Rect m_rect;

    // our view
    Camera m_camera;

    // scene description
    const Scene* m_scene;

    // output, executing process must supply filtered radiance values
    Colour* m_output;
    uint32_t  m_outputPitch;

    // quality settings
    uint32_t m_samplesPerPixel;

};

void ValidateC(const Colour& c, uint32_t x, uint32_t y)
{
	if (!_finite(c.r) ||
		!_finite(c.g) ||
		!_finite(c.b) ||
		!_finite(c.a) ||
		_isnan(c.r) ||
		_isnan(c.g) ||
		_isnan(c.b) ||
		_isnan(c.a))
	{
		cout << "Failed validation at " << x << "," << y << " iteration: " << g_iterations << endl;
	}
}

uint32_t RayTraceThreadFunc(void* data)
{
    const RayJob& job = *((RayJob*)(data));

    const Rect& rect = job.m_rect;
    float t;
    Vec3 n;

    for (uint32_t j=rect.m_top; j < rect.m_bottom; ++j)
    {
        for (uint32_t i=rect.m_left; i < rect.m_right; ++i)
        {
            Point3 origin;
            Vec3 dir;

            // generate a ray         
            switch (g_mode)
            {
                case ePathTrace:
                {
					job.m_camera.GenerateRay(i, j, origin, dir);

			        job.m_output[(g_height-1-j)*g_width+i] += PathTrace(*job.m_scene, origin, dir);
					//ValidateC(job.m_output[(g_height-1-j)*g_width+i], i, j);
                    break;
                }
                case eNormals:
                {
					job.m_camera.GenerateRayNoJitter(i, j, origin, dir);

                    const Primitive* p;
                    if (job.m_scene->Trace(origin, dir, t, n, &p))
                    {
                        n = n*0.5f+0.5f;
                        job.m_output[(g_height-1-j)*g_width+i] = Colour(n.x, n.y, n.z, 1.0f);
                    }
                    else
                    {
                        job.m_output[(g_height-1-j)*g_width+i] = Colour(0.0f);
                    }
                    break;
                }
#if _WIN32				
                case eComplexity:
                {
					job.m_camera.GenerateRayNoJitter(i, j, origin, dir);

                    const Primitive* p;
                    job.m_scene->Trace(origin, dir, t, n, &p);

                    // visualise traversal
                    job.m_output[(g_height-1-j)*g_width+i] = Colour(AABBTree::GetTraceDepth() / 100.0f);

                    break;
                }
#endif				
            }
        }
    }

    return 0;
}


void RayTrace()
{
    const uint32_t kSamples = g_width*g_height;

    ThreadPool threadPool;

	const uint32_t kTileSize = 16;
    const uint32_t kNumJobs = ceil(g_width / float(kTileSize)) * ceil(g_height / float(kTileSize));
	vector<RayJob> jobs(kNumJobs);

    // create camera
    const Camera camera(TransformMatrix(g_camDir, g_camPos), 45.0f, 0.1f, 10000.0f, g_width, g_height);
    uint32_t jobIndex = 0;

    double startTime = GetSeconds();

	// create thread jobs
    for (uint32_t y=0; y < g_height; y += kTileSize)
    {
		uint32_t top = y;
		uint32_t bottom = min(top+kTileSize, g_height);

		for (uint32_t x=0; x < g_width; x += kTileSize)
		{
			RayJob& job = jobs[jobIndex++];

            job.m_rect = Rect(x, min(x+kTileSize, g_width), top, bottom);
            job.m_camera = camera;
            job.m_scene = g_scene;
            job.m_samplesPerPixel = 1;
            job.m_output = &g_pixels[0];
            job.m_outputPitch = g_width;

			threadPool.AddTask(RayTraceThreadFunc, &job);
		}
    }

	threadPool.Run(g_numWorkers);
    threadPool.Wait();

    // print out trace time every 10 frames
    double endTime = GetSeconds();

	/*
    cout << "Nodes checked: " << g_nodesChecked << endl;
    cout << "Tris checked: " << g_trisChecked << endl;

    g_trisChecked = 0;
    g_nodesChecked = 0;
	*/

    ++g_iterations;

    Colour* presentMem = g_pixels;

    if (g_mode == ePathTrace)
    {
        float s = g_exposure / g_iterations;

        for (uint32_t i=0; i < g_width*g_height; ++i)
        {
            g_filtered[i] = LinearToSrgb(g_pixels[i] * s);
        }

        presentMem = g_filtered;
    }

	static uint32_t s_counter=0;
    if (s_counter % 10)
    {
        cout << "Trace took: " << (endTime-startTime)*1000.0f << "ms" << " rays/s: " << g_width*g_height/(endTime-startTime) << endl;
    }
    ++s_counter;

    glDisable(GL_BLEND);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    
    glPixelZoom(float(g_windowWidth)/g_width, float(g_windowHeight)/g_height);
    glDrawPixels(g_width,g_height,GL_RGBA,GL_FLOAT, presentMem);
}

void InitFrameBuffer()
{
    delete[] g_pixels;
    delete[] g_filtered;

    g_width = g_windowWidth*g_zoom;
    g_height = g_windowHeight*g_zoom;

    g_pixels = new Colour[g_width*g_height];
    g_filtered = new Colour[g_width*g_height];

    g_iterations = 0;
}

DiscPrimitive* g_sun;

void Init()
{
    g_scene = new Scene();

	float yoffset = 0.0f;

    // default materials
    Material* lambert = new MatteMaterial(new ConstantTexture(Colour(0.6f, 0.55f, 0.53f)));	
    Material* blinnPhong = new PlasticMaterial(new ConstantTexture(Colour(0.5f, 0.3f, 0.3f)), 120.0f);
	
	Material* gold = new PlasticMaterial(new ConstantTexture(Colour(1.0f, 0.71f, 0.29f, 120.0f)));
	Material* silver = new PlasticMaterial(new ConstantTexture(Colour(0.95f, 0.93f, 0.88f, 120.0f)));
	Material* copper = new PlasticMaterial(new ConstantTexture(Colour(0.95f, 0.64f, 0.54f, 120.0f)));
	Material* iron = new PlasticMaterial(new ConstantTexture(Colour(0.56f, 0.57f, 0.58f, 120.0f)));
	Material* aluminum = new PlasticMaterial(new ConstantTexture(Colour(0.91f, 0.92f, 0.92f, 120.0f)));
	Material* plaster = new MatteMaterial(new ConstantTexture(Colour(0.94f, 0.94f, 0.94f)));

	//Material* checker = new MatteMaterial(new CheckerboardTexture(0.1f*edgeLength.x, Colour(0.2f, 0.2f, 0.2f), Colour(0.5f, 0.5f, 0.5f)));


    //Mesh* mesh = ImportMeshFromObj("../../data/happy.obj");
	Mesh* mesh = ImportMeshFromPly("../../data/thearena.ply");
    //Mesh* mesh = ImportMeshFromPly("models/bunny/reconstruction/bun_zipper_res4.ply");
    //Mesh* mesh = ImportMeshFromPly("models/happy_recon/happy_vrip_res3.ply");
    //Mesh* mesh = ImportMeshFromPly("models/dragon/xyzrgb_dragon.ply"); yoffset = 22.1f;
	//Mesh* mesh = ImportMeshFromPly("models/xyzrgb_statuette.ply"); 
    //Mesh* mesh = ImportMeshFromObj("models/elephant/elefant_from_jotero_com.obj"); mesh->Transform(RotationMatrix(DegToRad(-150.0f), Vector3(1.0f, 0.0f, 0.0f)));
    //Mesh* mesh = ImportMeshFromObj("models/ajax/Ajax_Jotero_com.obj");
    //Mesh* mesh = ImportMeshFromObj("models/aphrodite/Aphrodite_from_jotero_com.obj");
    //Mesh* mesh = ImportMeshFromObj("models/figure/figure_from_jotero_com.obj");
    //Mesh* mesh = ImportMeshFromObj("models/sanktkilian/Sankt_Kilian_Sockel_jotero_com.obj"); mesh->Transform(RotationMatrix(DegToRad(-90.0f), Vector3(1.0f, 0.0f, 0.0f)));
    //Mesh* mesh = ImportMeshFromPly("models/lucy.ply"); mesh->Transform(RotationMatrix(DegToRad(-90.0f), Vector3(1.0f, 0.0f, 0.0f)));

	//Mesh* lightMesh = ImportMeshFromPly("models/bunny/reconstruction/bun_zipper_res4.ply");
	//lightMesh->Transform(TranslationMatrix(Point3(-0.0f, 0.5f, 0.0f)));

	
	// import classroom
	/*
	vector<string> models;
	FileScan("models/classroom/*.ply", models);

	Mesh* mesh = new Mesh();

	for (uint32_t i=0; i < models.size(); ++i)
	{
		Mesh* submesh = ImportMeshFromPly(("models/classroom/" + models[i]).c_str());
		submesh->Transform(RotationMatrix(DegToRad(-90.0f), Vector3(1.0f, 0.0f, 0.0f)));

		if (models[i] == "sun.ply")
		{
			MeshPrimitive* m = new MeshPrimitive(lambert, submesh); m->m_emission = Colour(9.0f, 8.5f, 8.0f)*100.0f;
			Light* l = new Light(Light::eArea, 1, m);

			g_scene->AddPrimitive(m);
			g_scene->AddLight(l);
		}
		else
		{			
			mesh->AddMesh(*submesh);
		}
	}
	*/
	

    // move at 5% of diagonal
    Vector3 minExtent, maxExtent, edgeLength, center;
    mesh->GetBounds(minExtent, maxExtent);
    
	// offset by 
	maxExtent.y += yoffset;
	minExtent.y += yoffset;
	edgeLength = 0.5f*(maxExtent-minExtent);

	center = 0.5f*(minExtent+maxExtent);

    // set mouse speed based on mesh size
    g_flySpeed = 0.05f*Length(maxExtent-minExtent);   

    // add primitives
    Primitive* plane = new PlanePrimitive(lambert, Point3(0.0f, minExtent.y, 0.0f), Vector3(0.0f, 1.0f, 0.0f));
    Primitive* sphere = new SpherePrimitive(blinnPhong, Point3(0.0f, 1.0f, 0.0f), 1.0f);
    Primitive* disc = new DiscPrimitive(lambert, Point3(center.x, 3.0f*maxExtent.y, -2.0f*maxExtent.z), Vector3(0.0f, -1.0f, 0.0f), 1.0f*sqrt(edgeLength.x*edgeLength.x + edgeLength.y*edgeLength.y));    disc->m_emission = Colour(9.0f, 8.5f, 8.0f)*3000.0f;
 	Primitive* model = new MeshPrimitive(lambert, mesh);

    Light* light = new Light(Light::eArea, 1, disc);
    
    //g_scene->AddPrimitive(plane);
    g_scene->AddPrimitive(model);
	//g_scene->AddPrimitive(lightBun);
    g_scene->AddPrimitive(disc);
    //g_scene->AddPrimitive(sphere);
	g_scene->AddLight(light);
    
	g_sun = (DiscPrimitive*)disc;

    // set up camera
    g_camPos = Point3(center.x, center.y, center.z+0.5f);
    g_camDir = Rotation(0.0f, 0.0f, 0.0f);
    
    InitFrameBuffer();
}

void GLUTUpdate()
{
    RayTrace();

	// flip
	glutSwapBuffers();
}

void GLUTReshape(int width, int height)
{
    g_windowWidth = width;
    g_windowHeight = height;

    InitFrameBuffer();
}

void GLUTArrowKeys(int key, int x, int y)
{
}

void GLUTArrowKeysUp(int key, int x, int y)
{
}

void GLUTKeyboardDown(unsigned char key, int x, int y)
{
    Matrix44 v = TransformMatrix(g_camDir, g_camPos);

    bool resetFrame = false;

 	switch (key)
	{
    case 'w':
        g_camPos -= Vector3(v.columns[2])*g_flySpeed; resetFrame = true;
		break;
    case 's':
        g_camPos += Vector3(v.columns[2])*g_flySpeed; resetFrame = true; 
        break;
    case 'a':
        g_camPos -= Vector3(v.columns[0])*g_flySpeed; resetFrame = true;
        break;
    case 'd':
        g_camPos += Vector3(v.columns[0])*g_flySpeed; resetFrame = true;
        break;
	case '1':
		g_mode = eNormals;
		break;
	case '2':
		g_mode = eComplexity;
		break;
    case '3':
        g_mode = ePathTrace; resetFrame = true;
        break;
    case '+':
        g_zoom = min(1.0f, g_zoom+0.30f); resetFrame = true;
        break;
    case '-':
        g_zoom = max(0.1f, g_zoom-0.30f); resetFrame = true;
        break;
	case '[':
		g_exposure -= 0.01f;
		break;
	case ']':
		g_exposure += 0.01f;
		break;
	case '8': 
		g_sunTheta += DegToRad(1.0f); resetFrame = true;
		break;
	case '5':
		g_sunTheta -= DegToRad(1.0f); resetFrame = true;
		break;
	case '9': 
		g_sunPhi += DegToRad(1.0f); resetFrame = true;
		break;
	case '6':
		g_sunPhi -= DegToRad(1.0f); resetFrame = true;
		break;
	case '7':
		g_sunTurbidity += 0.01f; resetFrame = true;
		break;
	case '4':
		g_sunTurbidity -= 0.01f; resetFrame = true;
		break;
    case '0':
        g_zoom = 1.0f; resetFrame = true;
        break;
	case 27:
		exit(0);
		break;
	};

    // reset image if there are any camera changes
    if (resetFrame == true)
    {
		// update sun position
		Vector3 sunDir = SphericalToXYZ(g_sunTheta, g_sunPhi);
		Point3 sunPos = sunDir*100000.0f;

		g_sun->m_plane = Plane(sunPos, -sunDir);
		g_sun->m_position = sunPos;
		g_sun->m_objToWorld = TransformFromVector(sunDir, sunPos);

		g_scene->SetSkyParams(g_sunTheta, g_sunPhi, g_sunTurbidity);

        InitFrameBuffer();
    }
}

void GLUTKeyboardUp(unsigned char key, int x, int y)
{
// 	switch (key)
// 	{
// 	case 27:
// 		exit(0);
// 		break;
// 	};

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

    if (g_mode == ePathTrace)
    {
        InitFrameBuffer();
    }
}


/*
void Application::JoystickFunc(int x, int y, int z, unsigned long buttons)
{
g_app->JoystickFunc(x, y, z, buttons);
}
*/

int main(int argc, char* argv[])
{	
    RandInit();

	// init gl
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);

	glutInitWindowSize(g_width, g_height);
	glutCreateWindow("NanoRay");
	glutPositionWindow(200, 200);

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

	glutMainLoop();
}



 