#include <iostream>

#include "core/types.h"
#include "core/platform.h"
#include "core/aabbtree.h"
#include "core/mesh.h"
#include "core/meshuv.h"
#include "core/tga.h"

#include "graphics/rendergl/glutil.h"

ChartSet g_charts;
FaceSet g_faces;
EdgeSet g_edges;

using namespace std;

const uint32_t g_screenWidth = 512;
const uint32_t g_screenHeight = 512;

Matrix44 g_worldToView;
Matrix44 g_viewToScreen;
Matrix44 g_screenToRaster;
Matrix44 g_worldToRaster;

Mesh* g_mesh;

Point3 g_camPos = Point3(0.0f, -1.0f, -1000.5f);
Rotation g_camDir = Rotation(180.0f, 0.0f, 0.0f);

int g_highlightChart = 0;
GLuint g_checkerboard;

typedef float float2[2];
typedef float float3[3];

template <typename OutputFunc>
void RasterizeTri(const float2 p, const float2 q, const float2 r, OutputFunc f)
{
	float2 minP, maxP;
	minP[0] = min(min(p[0], q[0]), r[0]);
	minP[1] = min(min(p[1], q[1]), r[1]);
	
	maxP[0] = max(max(p[0], q[0]), r[0]);
	maxP[1] = max(max(p[1], q[1]), r[1]);

	int xmin = int(floorf(minP[0]));
	int ymin = int(floorf(minP[1]));

	int xmax = int(ceilf(maxP[0]));
	int ymax = int(ceilf(maxP[1]));

	// compute inverse of tbe 2x2
	float a = q[0]-p[0];
	float b = r[0]-p[0];
	float c = q[1]-p[1];
	float d = r[1]-p[1];
	
	float rdet = 1.0f / (a*d - b*c);

	float m11 = d * rdet;
	float m21 = -c * rdet;
	float m12 = -b * rdet;
	float m22 = a * rdet;
	
	bool zeroSamples = true;
	
	// switch to conservative rasterization if not covering any pixels	
	for (int y=ymin; y <= ymax; ++y)
	{
		for (int x=xmin; x <= xmax; ++x)
		{
			// transform 
			float px = x + 0.5f - p[0];
			float py = y + 0.5f - p[1];
			
			// compute barycentric coordinates
			float alpha = px*m11 + py*m12;
			float beta = px*m21 + py*m22;
			float gamma = 1.0f - alpha - beta;
			
			if (alpha >= 0.0f && beta >= 0.0f && gamma >= 0.0f)
			{
				f(x, y, alpha, beta);	
				zeroSamples = false;
			}
		}
	}
	
	if (zeroSamples)
	{
		RasterizeTriConservative(p, q, r, f);
	}
}

template <typename OutputFunc>
void RasterizeTriConservative(const float2 p, const float2 q, const float2 r, OutputFunc f)
{
	// check winding and flip to make ccw
	float flip = p[0]*q[1] + q[0]*r[1] + r[0]*p[1] - p[0]*r[1] - q[0]*p[1] - r[0]*q[1];
	
	float2 minP, maxP;
	minP[0] = min(min(p[0], q[0]), r[0]);
	minP[1] = min(min(p[1], q[1]), r[1]);
	
	maxP[0] = max(max(p[0], q[0]), r[0]);
	maxP[1] = max(max(p[1], q[1]), r[1]);
	
	int xmin = int(floorf(minP[0]));
	int ymin = int(floorf(minP[1]));
	
	int xmax = int(ceilf(maxP[0]));
	int ymax = int(ceilf(maxP[1]));
	
	// construct edge functions	
	float2 n1 = { -flip*(q[1]-p[1]), flip*(q[0]-p[0]) };
	float2 n2 = { -flip*(r[1]-q[1]), flip*(r[0]-q[0]) };
	float2 n3 = { -flip*(p[1]-r[1]), flip*(p[0]-r[0]) };
	
	float2 t1 = { n1[0]<0.0f?0.0f:1.0f, n1[1]<0.0f?0.0f:1.0f };
	float2 t2 = { n2[0]<0.0f?0.0f:1.0f, n2[1]<0.0f?0.0f:1.0f };
	float2 t3 = { n3[0]<0.0f?0.0f:1.0f, n3[1]<0.0f?0.0f:1.0f };
					
	float d1 = -(n1[0]*p[0] + n1[1]*p[1]) + n1[0]*t1[0] + n1[1]*t1[1];
	float d2 = -(n2[0]*q[0] + n2[1]*q[1]) + n2[0]*t2[0] + n2[1]*t2[1];
	float d3 = -(n3[0]*r[0] + n3[1]*r[1]) + n3[0]*t3[0] + n3[1]*t3[1];

	// compute inverse of tbe 2x2
	float a = q[0]-p[0];
	float b = r[0]-p[0];
	float c = q[1]-p[1];
	float d = r[1]-p[1];
	
	float rdet = 1.0f / (a*d - b*c);

	float m11 = d * rdet;
	float m21 = -c * rdet;
	float m12 = -b * rdet;
	float m22 = a * rdet;
	
	for (int y=ymin; y <= ymax; ++y)
	{
		// pixel center
		float2 sp = { xmin + 0.5f, y + 0.5f };

		// set up scan line edge functions and barycentric coordinates
		float e1 = sp[0]*n1[0] + sp[1]*n1[1] + d1;
		float e2 = sp[0]*n2[0] + sp[1]*n2[1] + d2;
		float e3 = sp[0]*n3[0] + sp[1]*n3[1] + d3;

		float alpha = (sp[0]-p[0])*m11 + (sp[1]-p[1])*m12;
		float beta = (sp[0]-p[0])*m21 + (sp[1]-p[1])*m22;

		for (int x=xmin; x <= xmax; ++x)
		{
			// test if point is in the positive half space of each edge
			if (e1 >= 0.0f && e2 >= 0.0f && e3 >= 0.0f)
			{
				f(x, y, alpha, beta);	
			}

			e1 += n1[0];
			e2 += n2[0];
			e3 += n3[0];

			alpha += m11;
			beta += m21;
		}
	}
}


const uint32_t g_bakeSize = 512;

uint32_t g_buffer[g_bakeSize*g_bakeSize];
Face* g_currentFace;

void WritePixel(uint32_t x, uint32_t y, float a, float b)
{
	if (x < g_bakeSize && y < g_bakeSize)
	{
		Vector3 n1 = g_mesh->m_normals[g_currentFace->v[0]];
		Vector3 n2 = g_mesh->m_normals[g_currentFace->v[1]];
		Vector3 n3 = g_mesh->m_normals[g_currentFace->v[2]];
		
		a = Max(a, 0.0f);
		b = Max(b, 0.0f);
		float c = Max(c, 1.0f-a-b);

		Vector3 n = (a*n1 + b*n2 + c*n3) * 0.5f + Vector3(0.5f);

		//g_buffer[y*g_bakeSize + x] = ColourToRGBA8(Colour(n.x, n.y, n.z, 0.0f));
		g_buffer[y*g_bakeSize + x] += ColourToRGBA8(Colour(0.0f, 0.0f, 0.0f, 0.1f));
		
		/*
		Colour q(0.2f, 0.2f, 0.2f);
		if ((x+y)&1)
			q = Colour(1.0f, 1.0f, 1.0f);
				
		g_buffer[y*g_bakeSize + x] = ColourToRGBA8(q);//g_bakeColor;//c; //0xffffffff;
		*/
	}
}

void Dilate(uint32_t width, uint32_t height, uint32_t* tex)
{
	uint32_t* out = new uint32_t[width*height];
	
	int offsets[8][2] = { {  0,  1 },
						  {  0, -1 },
						  {  1,  0 },
						  { -1,  0 },
						  { -1, -1 },
						  { -1,  1 },
					      {  1, -1 },
						  {  1,  1 } };
	
	for (int y=0; y < height; ++y)
	{
		for (int x=0; x < width; ++x)
		{	
			uint32_t& src = tex[y*width + x];
			uint32_t& dst = out[y*width + x];
			
			dst = src;
			
			// if sample is empty
			if (dst&0xff000000)
				continue;
			
			for (int i=0; i < 8; ++i)
			{
				int px = x + offsets[i][0];
				int py = y + offsets[i][1];
								
				if (px > 0 && px < width && py > 0 && py < height)
				{
					uint32_t tap = tex[py*width + px];
					
					if (tap&0xff000000)
					{
						dst = tap;
						break;
					}
				}
			}
		}
	}
						
	memcpy(tex, out, sizeof(uint32_t)*width*height);
	
	delete out;
}

void BakeNormals()
{
	Vector2 scale(g_bakeSize, g_bakeSize);

	//for (ChartSet::const_iterator cIt=g_charts.begin(), cEnd=g_charts.end(); cIt != cEnd; ++cIt)
	for (std::vector<Chart*>::const_iterator cIt=g_charts.begin(), cEnd=g_charts.end(); cIt != cEnd; ++cIt)
	{
		const Chart& chart = *(*cIt);
	
		// draw each polygon in the chart
		for (std::set<Face*>::const_iterator fIt=chart.faces.begin(), fEnd=chart.faces.end(); fIt != fEnd; ++fIt)
		{
			g_currentFace = *fIt;
									
			RasterizeTri(g_mesh->m_texcoords[1][g_currentFace->v[0]]*scale,
				 					 g_mesh->m_texcoords[1][g_currentFace->v[1]]*scale,
									 g_mesh->m_texcoords[1][g_currentFace->v[2]]*scale, WritePixel);
		}
	}
	
	Dilate(g_bakeSize, g_bakeSize, g_buffer);
	Dilate(g_bakeSize, g_bakeSize, g_buffer);
	
	TgaImage img;
	img.m_width = g_bakeSize;
	img.m_height = g_bakeSize;
	img.m_data = g_buffer;
	TgaSave("test.tga", img);
}

void Init()
{
	RandInit();
	
	//gMesh = ImportMeshFromPly("../../bunny/reconstruction/bun_zipper_res4.ply");
	g_mesh = ImportMeshFromObj("../../sponza/sponza.obj");
	//g_mesh = CreateDiscMesh(500.0f, 100);

	double start = GetSeconds();
	
	UvAtlas(g_mesh, g_edges, g_faces, g_charts);
	
	cout << "UvAtlas took: " << (GetSeconds()-start)*1000.0f << "ms" << endl;
	start = GetSeconds();
	
	//PackChartsStrip(g_charts.begin(), g_charts.end());
	PackChartsBrute(g_charts.begin(), g_charts.end(), g_bakeSize, 0.5);
	
	cout << "Pack took: " << (GetSeconds()-start)*1000.0f << "ms" << endl;
	start = GetSeconds();
	
	BakeCharts(g_charts.begin(), g_charts.end(), g_mesh->m_texcoords[0], g_mesh->m_texcoords[1]);
	
	cout << "Bake took: " << (GetSeconds()-start)*1000.0f << "ms" << endl;
	
	BakeNormals();
	
	g_checkerboard = GlCreateTexture(g_bakeSize, g_bakeSize, g_buffer);
	
	//GlVerify(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST));
	//GlVerify(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST));
	GlVerify(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT));
	GlVerify(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT));
	
	cout << "Num charts: " << g_charts.size() << endl;
	cout << "Num faces: " << g_faces.size() << endl;
}

#define RASTER 0

void GLUTUpdate()
{	
	g_screenToRaster = ScaleMatrix(Vector3(g_screenWidth*0.5f, g_screenHeight*0.5f, 1.0f))*TranslationMatrix(Point3(1.0f, 1.0f, 0.0f));
	g_worldToView = AffineInverse(TransformMatrix(g_camDir, g_camPos));
	g_viewToScreen = ProjectionMatrix(45.0f, 1.0f, 10.0f, 10000.0f);
	g_worldToRaster = g_screenToRaster*g_viewToScreen*g_worldToView;	
	
	GlVerify(glDisable(GL_CULL_FACE));
	GlVerify(glEnable(GL_DEPTH_TEST));
	GlVerify(glDisable(GL_LIGHTING));
	GlVerify(glDisable(GL_BLEND));
	
	GlVerify(glViewport(0, 0, g_screenWidth, g_screenHeight));
	GlVerify(glClearColor(0.5f, 0.5f, 0.5f, 1.0f));
	GlVerify(glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT));
	//GlVerify(glPolygonMode(GL_FRONT_AND_BACK, GL_LINE));
	
#if !RASTER
	
	GlVerify(glMatrixMode(GL_PROJECTION));
	GlVerify(glLoadMatrixf((float*)&g_viewToScreen));
	
	GlVerify(glMatrixMode(GL_MODELVIEW));
	GlVerify(glLoadMatrixf((float*)&g_worldToView));
	/*
		
	GlVerify(glEnableClientState(GL_VERTEX_ARRAY));
	GlVerify(glEnableClientState(GL_COLOR_ARRAY));

	GlVerify(glVertexPointer(3, GL_FLOAT, sizeof(Vector3), &g_mesh->m_positions.front()))
	GlVerify(glColorPointer(3, GL_FLOAT, sizeof(Vector3), &g_mesh->m_normals.front()));

	GlVerify(glDrawElements(GL_TRIANGLES, g_mesh->m_indices.size(), GL_UNSIGNED_INT, &g_mesh->m_indices[0]));

	GlVerify(glDisableClientState(GL_VERTEX_ARRAY));
	GlVerify(glDisableClientState(GL_COLOR_ARRAY));
	 */
		
	/*
	glBegin(GL_LINES);
	glColor3f(0.0f, 1.0f, 0.0f);

	for (EdgeSet::iterator eIt=g_edges.begin(), eEnd=g_edges.end(); eIt != eEnd; ++eIt)
	{	
		glVertex3fv(g_mesh->m_positions[eIt->v[0]]);
		glVertex3fv(g_mesh->m_positions[eIt->v[1]]);
	}
	glEnd();
	*/
	
	
	const uint32_t kPaletteSize = 16;
	Colour palette[kPaletteSize] = { Colour(0x05B0A900u),
									 Colour(0xC9033700u),
									 Colour(0xFF4A3500u),
									 Colour(0x38F20300u),
									 Colour(0xC6CF5200u),
									 Colour(0xFCF5BF00u),
									 Colour(0x9FE5C200u),
								 	 Colour(0x5EB29900u),
	Colour(0x3B29c300u),
	Colour(0x5EB29900u),
	Colour(0x6F723900u),
	Colour(0x01B39600u),
	Colour(0x58623900u),
	Colour(0x2EB56800u),
	Colour(0x64B57200u),
	Colour(0x52029900u)};
	
	GlVerify(glEnable(GL_TEXTURE_2D));
	GlVerify(glBindTexture(GL_TEXTURE_2D, g_checkerboard));
	
	glBegin(GL_TRIANGLES);

	int s=0;

	for (std::vector<Chart*>::const_iterator cIt=g_charts.begin(), cEnd=g_charts.end(); cIt != cEnd; ++cIt)
	{
		const Chart& chart = *(*cIt);
		
		/*
		glBegin(GL_LINES);
		glColor3fv(palette[(i++)%3]);
		 
		// draw each edge of the boundary
		for (std::set<Edge*>::const_iterator eIt=chart.boundary.begin(), eEnd=chart.boundary.end(); eIt != eEnd; ++eIt)
		{
			Edge* e = *eIt;
			
			glVertex3fv(g_mesh->m_positions[(*eIt)->v[0]]);
			glVertex3fv(g_mesh->m_positions[(*eIt)->v[1]]);	
		}
		*/

		//glColor3fv(palette[(i++)%kPaletteSize]);
		glColor3f(1.0f, 1.0f, 1.0f);
		if (g_highlightChart == s++)
			glColor3f(1.0f, 1.0f, 1.0f);
		
		if (1)
		{
			// draw each polygon in the chart
			for (std::set<Face*>::const_iterator fIt=chart.faces.begin(), fEnd=chart.faces.end(); fIt != fEnd; ++fIt)
			{
				Face* f = *fIt;
				

				glTexCoord2fv(g_mesh->m_texcoords[1][f->v[0]]);
				glVertex3fv(g_mesh->m_positions[f->v[0]]);

				glTexCoord2fv(g_mesh->m_texcoords[1][f->v[1]]);
				glVertex3fv(g_mesh->m_positions[f->v[1]]);
				
				glTexCoord2fv(g_mesh->m_texcoords[1][f->v[2]]);
				glVertex3fv(g_mesh->m_positions[f->v[2]]);
			}
		}

		float scale = 100.0f;

		for (std::set<Face*>::const_iterator fIt=chart.faces.begin(), fEnd=chart.faces.end(); fIt != fEnd; ++fIt)
		{
			Face* f = *fIt;
			
			glVertex2fv(g_mesh->m_texcoords[1][f->v[0]]*scale);
			glVertex2fv(g_mesh->m_texcoords[1][f->v[1]]*scale);
			glVertex2fv(g_mesh->m_texcoords[1][f->v[2]]*scale);
		}
		 
		
	}
	
	glEnd();
	
#else

	memset(gBuffer, 0x80808080, sizeof(gBuffer));

	for (uint32_t i=0; i < g_mesh->GetNumFaces(); ++i)
	{
		uint32_t a = g_mesh->m_indices[i*3+0];
		uint32_t b = g_mesh->m_indices[i*3+1];
		uint32_t c = g_mesh->m_indices[i*3+2];
		
		Vector4 sa = g_worldToRaster*g_mesh->m_positions[a].operator Vec4();
		Vector4 sb = g_worldToRaster*g_mesh->m_positions[b].operator Vec4();
		Vector4 sc = g_worldToRaster*g_mesh->m_positions[c].operator Vec4();
		
		sa /= sa.w;
		sb /= sb.w;
		sc /= sc.w;
		
		RasterizeTri(sa, sb, sc, WritePixel);
	}
    
    glDrawPixels(g_screenWidth, g_screenHeight, GL_RGBA, GL_UNSIGNED_BYTE, gBuffer);

#endif
	
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
	Matrix44 v = TransformMatrix(g_camDir, g_camPos);
	
	const float kFlySpeeed = 10.0f;
	
 	switch (key)
	{
		case 'e':
		{
			break;
		}
		case 'w':
			g_camPos -= Vector3(v.columns[2])*kFlySpeeed;
			break;
		case 's':
			g_camPos += Vector3(v.columns[2])*kFlySpeeed;
			break;
		case 'a':
			g_camPos -= Vector3(v.columns[0])*kFlySpeeed;
			break;
		case 'd':
			g_camPos += Vector3(v.columns[0])*kFlySpeeed;
			break;
		case 'o':
			g_highlightChart++;
			break;
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
    glutCreateWindow("Radiosity");
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

