#include <iostream>

#include "core/types.h"
#include "core/shader.h"
#include "core/platform.h"

#define STRINGIFY(A) #A

using namespace std;

uint32_t g_screenWidth = 800;
uint32_t g_screenHeight = 600;

GLuint g_pointShader;
GLuint g_solidShader;

//--------------------------------------------------------
// Solid shaders
//
const char *vertexShader = STRINGIFY(

void main()
{
    // calculate window-space point size
	gl_Position = gl_ModelViewProjectionMatrix*vec4(gl_Vertex.xyz, 1.0);
	gl_TexCoord[0] = gl_ModelViewMatrixInverseTranspose*vec4(gl_Normal.xyz, 0.0);   
}
);

const char *fragmentShader = STRINGIFY(

void main()
{
	gl_FragColor = gl_TexCoord[0];
}
);

//--------------------------------------------------------
// Point shaders
//
const char *vertexPointShader = "#version 120\n"STRINGIFY(

// rotation matrix in xyz, scale in w
attribute vec4 q1;
attribute vec4 q2;
attribute vec4 q3;

// returns 1.0 for x==0.0 (unlike glsl)
float Sign(float x) { return x < 0.0 ? -1.0: 1.0; }

bool solveQuadratic(float a, float b, float c, out float minT, out float maxT)
{
	if (a == 0.0 && b == 0.0)
	{
		minT = maxT = 0.0;
		return false;
	}

	float discriminant = b*b - 4.0*a*c;

	if (discriminant < 0.0)
	{
		return false;
	}

	float t = -0.5*(b + Sign(b)*sqrt(discriminant));
	minT = t / a;
	maxT = c / t;

	if (minT > maxT)
	{
		float tmp = minT;
		minT = maxT;
		maxT = tmp;
	}

	return true;
}

float DotInvW(vec4 a, vec4 b) {	return a.x*b.x + a.y*b.y + a.z*b.z - a.w*b.w; }

void main()
{
	// construct quadric matrix
	mat4 q;
	q[0] = vec4(q1.xyz*q1.w, 0.0);
	q[1] = vec4(q2.xyz*q2.w, 0.0);
	q[2] = vec4(q3.xyz*q3.w, 0.0);
	q[3] = vec4(gl_Vertex.xyz, 1.0);

	// transforms a normal to parameter space (inverse transpose of (q*modelview)^-T)
	mat4 invClip = transpose(gl_ModelViewProjectionMatrix*q);

	// solve for the right hand bounds in homogenous clip space
	float a1 = DotInvW(invClip[3], invClip[3]);
	float b1 = -2.0f*DotInvW(invClip[0], invClip[3]);
	float c1 = DotInvW(invClip[0], invClip[0]);

	float xmin;
	float xmax;
 	solveQuadratic(a1, b1, c1, xmin, xmax);	

	// solve for the right hand bounds in homogenous clip space
	float a2 = DotInvW(invClip[3], invClip[3]);
	float b2 = -2.0f*DotInvW(invClip[1], invClip[3]);
	float c2 = DotInvW(invClip[1], invClip[1]); 

	float ymin;
	float ymax;
 	solveQuadratic(a2, b2, c2, ymin, ymax);

	gl_Position = vec4(gl_Vertex.xyz, 1.0);
	gl_TexCoord[0] = vec4(xmin, xmax, ymin, ymax);

	// construct inverse quadric matrix (used for ray-casting in parameter space)
	mat4 invq;
	invq[0].xyz = q1.xyz/q1.w;
	invq[1].xyz = q2.xyz/q2.w;
	invq[2].xyz = q3.xyz/q3.w;
	invq[3].w = 1.0;

	invq = transpose(invq);
	invq[3] = -(invq*gl_Position);

	// transform a point from view space to parameter space
	invq = invq*gl_ModelViewMatrixInverse;

	// pass down
	gl_TexCoord[1] = invq[0];
	gl_TexCoord[2] = invq[1];
	gl_TexCoord[3] = invq[2];
	gl_TexCoord[4] = invq[3];
}
);

const char* geometryPointShader = 
"#version 120\n"
"#extension GL_EXT_geometry_shader4 : enable\n"
STRINGIFY(
void main()
{
	vec3 pos = gl_PositionIn[0].xyz;
	vec4 bounds = gl_TexCoordIn[0][0];

	float xmin = bounds.x;
	float xmax = bounds.y;
	float ymin = bounds.z;
	float ymax = bounds.w;

	// inv quadric transform
	gl_TexCoord[0] = gl_TexCoordIn[0][1];
	gl_TexCoord[1] = gl_TexCoordIn[0][2];
	gl_TexCoord[2] = gl_TexCoordIn[0][3];
	gl_TexCoord[3] = gl_TexCoordIn[0][4];

	gl_Position = vec4(xmin, ymax, 0.0, 1.0);
	EmitVertex();

	gl_Position = vec4(xmin, ymin, 0.0, 1.0);
	EmitVertex();

	gl_Position = vec4(xmax, ymax, 0.0, 1.0);
	EmitVertex();

	gl_Position = vec4(xmax, ymin, 0.0, 1.0);
	EmitVertex();
}
);

// pixel shader for rendering points as shaded spheres
const char *fragmentPointShader = STRINGIFY(

uniform vec3 invViewport;
uniform vec3 invProjection;

float Sign(float x) { return x < 0.0 ? -1.0: 1.0; }

bool solveQuadratic(float a, float b, float c, out float minT, out float maxT)
{
	if (a == 0.0 && b == 0.0)
	{
		minT = maxT = 0.0;
		return true;
	}

	float discriminant = b*b - 4.0*a*c;

	if (discriminant < 0.0)
	{
		return false;
	}

	float t = -0.5*(b + Sign(b)*sqrt(discriminant));
	minT = t / a;
	maxT = c / t;

	if (minT > maxT)
	{
		float tmp = minT;
		minT = maxT;
		maxT = tmp;
	}

	return true;
}

float sqr(float x) { return x*x; }

void main()
{
	// transform from view space to parameter space
	mat4 invQuadric;
	invQuadric[0] = gl_TexCoord[0];
	invQuadric[1] = gl_TexCoord[1];
	invQuadric[2] = gl_TexCoord[2];
	invQuadric[3] = gl_TexCoord[3];

	vec4 ndcPos = vec4(gl_FragCoord.xy*invViewport.xy*vec2(2.0, 2.0) - vec2(1.0, 1.0), -1.0, 1.0);
	vec4 viewDir = gl_ProjectionMatrixInverse*ndcPos; 

	// ray to parameter space
	vec4 dir = invQuadric*vec4(viewDir.xyz, 0.0);
	vec4 origin = invQuadric[3];

	// set up quadratric equation
	float a = sqr(dir.x) + sqr(dir.y) + sqr(dir.z);// - sqr(dir.w);
	float b = dir.x*origin.x + dir.y*origin.y + dir.z*origin.z - dir.w*origin.w;
	float c = sqr(origin.x) + sqr(origin.y) + sqr(origin.z) - sqr(origin.w);

	float minT;
	float maxT;

	if (solveQuadratic(a, 2.0*b, c, minT, maxT))
	{
		vec3 hitPos = viewDir.xyz*minT;
		vec3 dx = dFdx(hitPos);
		vec3 dy = dFdy(hitPos);

		gl_FragColor.xyz = normalize(cross(dx, dy));
		gl_FragColor.w = 1.0;

		return;
	}
	//else
//		discard;	

	gl_FragColor = vec4(0.5, 0.0, 0.0, 1.0);
}


);

// dot product with negative w
float DotInvW(const Vec4& a, const Vec4& b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z - a.w*b.w;
}

void Init()
{
	int maxVerts;
	glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES_EXT, &maxVerts);
	printf("%d\n", maxVerts);
		
	g_pointShader = CompileProgram(vertexPointShader, fragmentPointShader, geometryPointShader);
	g_solidShader = CompileProgram(vertexShader, fragmentShader);
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
	gluLookAt(sinf(GetSeconds())*5.0f, cosf(GetSeconds())*5.0f, 10.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);

	float radius = 1.0f;
	float aspect = float(g_screenWidth)/g_screenHeight;

	Point3 quadricPos = Point3(1.f*radius, 0.0f, 0.0f);
	Vec3 quadricScale = Vec3(0.5f, 1.0f, 1.0f);

	Matrix44 R = RotationMatrix(0.0f, Vec3(0.0f, 0.0f, 1.0f));
	Matrix44 T = TranslationMatrix(quadricPos)*R; 
	Matrix44 S = ScaleMatrix(quadricScale);

	Matrix44 TInv = AffineInverse(T);
	Matrix44 SInv = ScaleMatrix(Vec3(1.0f/quadricScale.x, 1.0f/quadricScale.y, 1.0f/quadricScale.z));

	// world space to parameter space
	Matrix44 Q = T*S;
	Matrix44 QInv = SInv*TInv;

	glUseProgram(g_solidShader);
	glPushMatrix();
	glTranslatef(-radius, 0.0f, 0.0f);
	glScalef(1.0, 1.0f, 1.0);
	glutSolidSphere(radius, 20, 20);
	glPopMatrix();

	const float viewHeight = tanf(DegToRad(fov)/2.0f);

	glUseProgram(g_pointShader);
	glUniform1f( glGetUniformLocation(g_pointShader, "pointScale"), g_screenHeight/viewHeight);
	glUniform1f( glGetUniformLocation(g_pointShader, "pointRadius"), radius);
	glUniform3fv( glGetUniformLocation(g_pointShader, "invViewport"), 1, Vec3(1.0f/g_screenWidth, 1.0f/g_screenHeight, 1.0f));
	glUniform3fv( glGetUniformLocation(g_pointShader, "invProjection"), 1, Vec3(aspect*viewHeight, viewHeight, 1.0f));
	
	glDepthMask(GL_TRUE);
	glEnable(GL_DEPTH_TEST);		
	glDisable(GL_CULL_FACE);
	
	// pass in rotation and scale packed into 3x4 matrix
	Vec4 q1 = R.columns[0]; q1.w = quadricScale.x;
	Vec4 q2 = R.columns[1]; q2.w = quadricScale.y;
	Vec4 q3 = R.columns[2]; q3.w = quadricScale.z;	

	// ellipsoid eigenvectors
	int s1 = glGetAttribLocation(g_pointShader, "q1");
	assert(s1 != -1);
	glEnableVertexAttribArray(s1);
	glVertexAttribPointer(s1, 4, GL_FLOAT, GL_FALSE, 0, q1);

	int s2 = glGetAttribLocation(g_pointShader, "q2");
	assert(s2 != -1);
	glEnableVertexAttribArray(s2);
	glVertexAttribPointer(s2, 4, GL_FLOAT, GL_FALSE, 0, q2);

	int s3 = glGetAttribLocation(g_pointShader, "q3");
	assert(s3 != -1);
	glEnableVertexAttribArray(s3);
	glVertexAttribPointer(s3, 4, GL_FLOAT, GL_FALSE, 0, q3);
	
	glEnableClientState(GL_VERTEX_ARRAY);			
	glVertexPointer(3, GL_FLOAT, sizeof(float)*3, quadricPos);

	glDrawArrays(GL_POINTS, 0, 1);
	

	glDisableVertexAttribArray(s1);
	glDisableVertexAttribArray(s2);
	glDisableVertexAttribArray(s3);
	
	/*

	// build billboard
	Matrix44 view, projection;
	glGetFloatv(GL_MODELVIEW_MATRIX, (float*)&view);
	glGetFloatv(GL_PROJECTION_MATRIX, (float*)&projection);

	// transform a normal to parameter space
	Matrix44 invClip = Transpose(projection*view*Q);

	// solve for the right hand bounds in homogenous clip space
	float a1 = DotInvW(invClip.columns[3], invClip.columns[3]);
	float b1 = -2.0f*DotInvW(invClip.columns[0], invClip.columns[3]);
	float c1 = DotInvW(invClip.columns[0], invClip.columns[0]); 

	float xmin, xmax;
 	SolveQuadratic(a1, b1, c1, xmin, xmax);	

	// solve for the right hand bounds in homogenous clip space
	float a2 = DotInvW(invClip.columns[3], invClip.columns[3]);
	float b2 = -2.0f*DotInvW(invClip.columns[1], invClip.columns[3]);
	float c2 = DotInvW(invClip.columns[1], invClip.columns[1]); 

	float ymin, ymax;
 	SolveQuadratic(a2, b2, c2, ymin, ymax);	

	Vec4 quadVertices[4] = 
	{
	   	Vec4(xmin, ymax, 0.0f, 1.0f),
	   	Vec4(xmin, ymin, 0.0f, 1.0f),
	   	Vec4(xmax, ymin, 0.0f, 1.0f),
	   	Vec4(xmax, ymax, 0.0f, 1.0f)
	};

	glBegin(GL_QUADS);
	glVertex3fv(quadVertices[0]);
	glVertex3fv(quadVertices[1]);
	glVertex3fv(quadVertices[2]);
	glVertex3fv(quadVertices[3]);
	glEnd();	
		*/


	glUseProgram(0);
	
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

