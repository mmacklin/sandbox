#include <core/maths.h>
#include <core/pfm.h>
#include <core/mesh.h>
#include <core/shader.h>

#define STRINGIFY(A) #A

#include "sh.h"

#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <stdint.h>

// Work Items
// ----------
// Implement the spherical harmonic functions in C++
// Write basic unit tests for functions to check orthonormal properties, sanity check against mathematica
// Bring in light probe loading code from probe viewer app
// Write code to project probe into SH basis and convolve with Lambert function
// Write code to view projected map SH basis, check against references
// Write code to load and render teapot with projected probe

using namespace std;

const int kWidth = 800;
const int kHeight = 600;

uint32_t g_buffer[kWidth*kHeight];
Mesh* g_mesh;

Vec3 g_camPos(0.0f, 0.0f, 7.0f);
Vec3 g_camVel(0.0f);
Vec3 g_camAngle(0.0f, -kPi/12.0f, 0.0f);
Vec3 g_shDiffuse[9];
Vec3 g_shPhong[9];
float g_exposure = 0.175f;

GLuint g_shader;

// vertex shader
const char* vertexShader = STRINGIFY
(
	void main()
	{
		gl_Position = gl_ModelViewProjectionMatrix*vec4(gl_Vertex.xyz, 1.0);
		gl_TexCoord[0] = vec4(normalize(gl_Normal), 0.0);
		gl_TexCoord[1] = vec4(gl_Vertex.xyz, 1.0);
	}
);

// pixel shader
const char* fragmentShader = STRINGIFY
(
	uniform vec3 shDiffuse[9];
	uniform vec3 shPhong[9];
	uniform float exposure;

	vec3 shEval(vec3 dir, vec3 sh[9])
	{
		// evaluate irradiance
		vec3 e = sh[0];

		e += -dir.y*sh[1];
		e +=  dir.z*sh[2];
		e += -dir.x*sh[3];

		e +=  dir.x*dir.y*sh[4];
		e += -dir.y*dir.z*sh[5];
		e += -dir.x*dir.z*sh[7];

		e += (3.0*dir.z*dir.z-1.0)*sh[6];
		e += (dir.x*dir.x - dir.y*dir.y)*sh[8];

		return max(e, vec3(0.0));
	}

	vec3 fresnel(vec3 rf, float cosTheta)
	{
		return rf + pow(1.0-cosTheta, 5.0)*(vec3(1.0)-rf);
	}
	vec3 lerp(vec3 a, vec3 b, vec3 t)
	{
		return a + (b-a)*t;
	}

	void main()
	{
		vec3 rf = vec3(1.0, 0.71, 0.29);
		vec3 rd = rf*0.1;//vec3(0.2);

		vec3 n = gl_TexCoord[0].xyz;
		vec3 shadePos = gl_TexCoord[1].xyz;
		vec3 eyePos = gl_ModelViewMatrixInverse[3].xyz;
		vec3 eyeDir = normalize(eyePos-shadePos);

		// calculate reflected view direction
		vec3 r = normalize(2.0*dot(eyeDir, n)*n-eyeDir);

		// evaluate reflected radiance in reflected view dir
		vec3 s = shEval(r, shPhong);

		// calculate diffuse radiance in normal dir
		vec3 d = shEval(n, shDiffuse);
		
		// 
		vec3 f = fresnel(rf, clamp(dot(r, n), 0.0, 1.0));

		// lerp between two based on Fresnel function
		vec3 l = rd*d + f*s;
			
		gl_FragColor = vec4(pow(max(l*exposure, 0.0), vec3(1.0/2.0)), 1.0);
	}
);


double lambert(double theta, double phi)
{
	return max(cos(theta), 0.0);
}

double phong(double theta, double phi)
{
	double kExponent = 40.0f;

	return 0.5*(2.0 + kExponent)*pow(lambert(theta, phi), kExponent);
}

Colour light(double theta, double phi)
{
	Vec3 dir(float(sin(theta)*cos(phi)), float(sin(theta)*sin(phi)), float(cos(theta)));

	float f = Dot(dir, Normalize(Vec3(1.0f, 6.0f, 4.0f))) > cos(kPi/25.0);

	return Colour(f, f, f, 1.0f);
}

struct ProbeSampler
{
	ProbeSampler(const char* path)
	{
		if (!PfmLoad(path, mProbe))
		{
			cout << "Could not load probe " << path << endl;
		}
		else
		{
			cout << "Loaded probe: " << path << " (" << mProbe.m_width << ", " << mProbe.m_height << ")" << endl;
		}
	}

	Colour operator()(double theta, double phi) const
	{
		Vec3 dir(float(sin(theta)*cos(phi)), float(sin(theta)*sin(phi)), float(cos(theta)));

		// convert world space dir to probe space
		float c = (1.0f / kPi) * acosf(dir.z)/sqrt(dir.x*dir.x + dir.y*dir.y);
		
		uint32_t px = uint32_t((0.5f + 0.5f*(dir.x*c))*(mProbe.m_width-1));
		uint32_t py = uint32_t((0.5f + 0.5f*(dir.y*c))*(mProbe.m_height-1));
		
		float r = mProbe.m_data[py*mProbe.m_width*3 + (px*3)+0];
		float g = mProbe.m_data[py*mProbe.m_width*3 + (px*3)+1];
		float b = mProbe.m_data[py*mProbe.m_width*3 + (px*3)+2];

		return Colour(r, g, b, 1.0f);
	}

	PfmImage mProbe;
};

ostream& operator<<(ostream& s, const Colour& c) { s << c.r << ", " << c.g << ", " << c.b << ", " << c.a; return s; }

void shScale(const Colour in[9], Vec3 out[9])
{
	const float k0 = 0.5f * 1.0f/sqrtf(kPi); 
	const float k1 = 0.5f * sqrtf(3.0f / kPi);
	const float k2 = 0.5f * sqrtf(15.0f / kPi);
	const float k3 = 0.25f * sqrtf(5.0f / kPi);
	const float k4 = 0.25f * sqrtf(15.0f / kPi);

	// copy to global probe coefficients, pre-scaled to avoid work in shader
	out[0] = k0*Vec3(in[0]);

	out[1] = k1*Vec3(in[1]);
	out[2] = k1*Vec3(in[2]);
	out[3] = k1*Vec3(in[3]);

	out[4] = k2*Vec3(in[4]);
	out[5] = k2*Vec3(in[5]);
	out[6] = k3*Vec3(in[6]); 
	out[7] = k2*Vec3(in[7]);
	out[8] = k4*Vec3(in[8]);

	// dump to C-array format
	cout << "Vec3 out[] = {" << endl;
	for (int i=0; i < 9; ++i)
	{
		cout << "Vec3(" << out[i].x << ", " << out[i].y << ", " << out[i].z << ")," << endl;
	}
	cout << "}" << endl;
}

void shTests()
{
	const int lmax = 3;

	double lambertCoefficients[lmax*lmax] = { 0.0 };
	double phongCoefficients[lmax*lmax] = { 0.0 };

	shProject(lambert, lmax, lambertCoefficients); 
	shProject(phong, lmax, phongCoefficients);

	ProbeSampler sampler("../../data/beach.pfm");

	Colour probeCoefficients[lmax*lmax];
	shProject(sampler, lmax, probeCoefficients);

	Colour diffuseCoefficients[lmax*lmax];
	shConvolve(diffuseCoefficients, probeCoefficients, lambertCoefficients, lmax);
	shScale(diffuseCoefficients, g_shDiffuse);

	Colour specularCoefficients[lmax*lmax];
	shConvolve(specularCoefficients, probeCoefficients, phongCoefficients, lmax);
	shScale(specularCoefficients, g_shPhong);
//	shReduceRinging(g_shPhong, lmax, 0.01);

/*
	for (int i=0; i < (lmax*lmax); ++i)
		cout << probeCoefficients[i] << endl;

	// expand
	for (int y=0; y < kHeight; ++y)
	{
		for (int x=0; x < kWidth; ++x)
		{
			// calculate direction
			double u = 2.0*(double(x)/kWidth)-1.0;
			double v = 2.0*(double(y)/kHeight)-1.0;

			double theta = atan2(v, u);
			double phi = kPi*sqrt(u*u+v*v);

			// do not write pixels outside the sphere
			if (phi <= kPi)
			{
				Vec3 direction = SphericalToXYZ(theta, phi)*0.5 + Vec3(0.5);
				Colour c(direction.x, direction.y, direction.z);

				g_buffer[y*kWidth + x] = ColourToRGBA8(0.175*(shExpand(probeCoefficients, lmax, phi, theta)));//sampler(theta, phi));
			}
		}
	}
*/
}

void Update(void)
{
	const float dt = 1.0f/60.0f;

	// update camera
	const Vec3 forward(-sinf(g_camAngle.x)*cosf(g_camAngle.y), sinf(g_camAngle.y), -cosf(g_camAngle.x)*cosf(g_camAngle.y));
	const Vec3 right(Normalize(Cross(forward, Vec3(0.0f, 1.0f, 0.0f))));
	
	g_camPos += (forward*g_camVel.z + right*g_camVel.x)*dt;

	glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	/*
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0f, 1.0f, -2.0f, 2.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	glDrawPixels(kWidth, kHeight, GL_RGBA, GL_UNSIGNED_BYTE, g_buffer);
	*/

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0f, float(kWidth)/kHeight, 1.0f, 10000.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(RadToDeg(-g_camAngle.x), 0.0f, 1.0f, 0.0f);
	glRotatef(RadToDeg(-g_camAngle.y), cosf(-g_camAngle.x), 0.0f, sinf(-g_camAngle.x));	
	glTranslatef(-g_camPos.x, -g_camPos.y, -g_camPos.z);
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, &g_mesh->m_positions.front());
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, 0, &g_mesh->m_normals.front());

	glUseProgram(g_shader);
	GLint uDiffuse = glGetUniformLocation(g_shader, "shDiffuse");
	glUniform3fv(uDiffuse, 9, reinterpret_cast<float*>(&g_shDiffuse[0].x));
	
	GLint uPhong = glGetUniformLocation(g_shader, "shPhong");
	glUniform3fv(uPhong, 9, reinterpret_cast<float*>(&g_shPhong[0].x));

	GLint uExposure = glGetUniformLocation(g_shader, "exposure");
	glUniform1f(uExposure, g_exposure);

	glDrawElements(GL_TRIANGLES, g_mesh->GetNumFaces()*3, GL_UNSIGNED_INT, &g_mesh->m_indices.front());

	glUseProgram(0);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);

	glutSwapBuffers();

	g_mesh->Transform(RotationMatrix(kPi/180, Vec3(0.0f, 1.0f, 0.0f)));
}

int lastx = 0;
int lasty = 0;

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
	};
}

void GLUTKeyboardDown(unsigned char key, int x, int y)
{
	const float kSpeed = 50.0f;

	switch (key)
	{
		case 'w':
		{
			g_camVel.z = kSpeed;
			break;
		}
		case 's':
		{
			g_camVel.z = -kSpeed;
			break;
		}
		case 'a':
		{
			g_camVel.x = -kSpeed;
			break;
		}
		case 'd':
		{
			g_camVel.x = kSpeed;
			break;
		}
		case 'u':
		{
			g_exposure += 0.05f;
			break;
		}
		case 'j':
		{
			g_exposure -= 0.05f;
			break;
		}
		case 27:
		case 'q':
		{
			exit(0);
			break;
		}
	};
}

void GLUTKeyboardUp(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'w':
		{
			g_camVel.z = 0.0f;
			break;
		}
		case 's':
		{
			g_camVel.z = 0.0f;
			break;
		}
		case 'a':
		{
			g_camVel.x = 0.0f;
			break;
		}
		case 'd':
		{
			g_camVel.x = 0.0f;
			break;
		}
	};
}

void GLUTMotionFunc(int x, int y)
{	
	int dx = x-lastx;
	int dy = y-lasty;
	
	lastx = x;
	lasty = y;

	const float kSensitivity = DegToRad(0.2f);

	g_camAngle.x -= dx*kSensitivity;
	g_camAngle.y += dy*kSensitivity;
}


int main(int argc, char* argv[])
{
	// init gl
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
	
	glutInitWindowSize(kWidth, kHeight);
	glutCreateWindow("SHTest");
	glutPositionWindow(200, 100);

#if _WIN32
	glewInit();
#endif

	g_shader = CompileProgram(vertexShader, fragmentShader);
	if (g_shader == 0)
		return 0;

	g_mesh = ImportMeshFromObj("../../data/happy.obj");

	Vec3 minExtents, maxExtents;
	g_mesh->GetBounds(minExtents, maxExtents);
	
	Vec3 center = 0.5f*(maxExtents+minExtents);
	float scale = 10.0f / (maxExtents.y-minExtents.y);

	g_mesh->Transform(ScaleMatrix(Vector3(scale, scale, scale))*TranslationMatrix(Point3(-center)));	

	// center camera
	g_camPos = Vec3(0.0f, 5.0f, 20.0f); 

	shTests();

	glutIdleFunc(Update);	
	glutDisplayFunc(Update);
	glutMouseFunc(GLUTMouseFunc);
	glutMotionFunc(GLUTMotionFunc);
	glutKeyboardFunc(GLUTKeyboardDown);
	glutKeyboardUpFunc(GLUTKeyboardUp);
/*	

	glutReshapeFunc(GLUTReshape);
	glutDisplayFunc(GLUTUpdate);
   
	glutSpecialFunc(GLUTArrowKeys);
	glutSpecialUpFunc(GLUTArrowKeysUp);

*/
#if __APPLE__
	int swap_interval = 1;
	CGLContextObj cgl_context = CGLGetCurrentContext();
	CGLSetParameter(cgl_context, kCGLCPSwapInterval, &swap_interval);
#endif

	glutMainLoop();
	return 0;
}





