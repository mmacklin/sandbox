#include <core/maths.h>
#include <core/shader.h>
#include <core/platform.h>
#include <core/shader.h>
#include <core/tga.h>
#include <core/mesh.h>

#define STRINGIFY(A) #A

#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <vector>
#include <stdint.h>

#include "tag.h"
#include "shaders.h"

#if __APPLE__
#include <mach-o/dyld.h>
#endif

using namespace std;

int gScreenWidth = 1280;
int gScreenHeight = 720;

// default file
string gAnimFile;
string gControlFile;

Vec3 gCamPos(0.0f);//, 150.0f, -357.0f);
Vec3 gCamVel(0.0f);
Vec3 gCamAngle(kPi, 0.0f, 0.0f);
float gTagWidth = 2.0f;
float gTagHeight = 4.0f;
Point3 gTagCenter;
Point3 gTagLower(-20.0f);
Point3 gTagUpper(20.0f);
float gTagSmoothing = 0.8f;
float gTagVelocityScale = 0.2f;
bool gShowHelp = true;
bool gShowCan = true;
bool gFullscreen = false;

vector<Brush*> gTagBrushes;
uint32_t gTagBrushIndex;

Mesh* gMeshCan;

GLuint gMainShader;
GLuint gDebugShader;
GLuint gShadowFrameBuffer;
GLuint gShadowTexture;

enum EventType
{
	eStopPaint,
	eStartPaint	
};

struct Control
{
	float time;
	EventType event;
};

vector<Control> gControlTrack; 
bool gControlRecord = false;
float gControlStartTime;

struct Frame
{
	Vec3 pos;
	Vec3 rot;
};

bool LoadBvh(const char* filename, std::vector<Frame>& frames, Point3& center, Point3& lower, Point3& upper)
{
	FILE* f = fopen(filename, "r");
	
	lower = Point3(FLT_MAX);
	upper = Point3(-FLT_MAX);

	if (f)
	{
		while (!feof(f))
		{
			Frame frame;
			int n = fscanf(f, "%f %f %f %f %f %f", 
					&frame.pos.x,	
					&frame.pos.y,	
					&frame.pos.z,	
					&frame.rot.z,	
					&frame.rot.x,	
					&frame.rot.y);	
	
			if (n == EOF)
				break;	

			if (n != 6)
			{
				char buf[1024];
				fgets(buf, 1024, f);
			}
			else
			{
				frames.push_back(frame);

				center += frame.pos;
				lower = Min(lower, Point3(frame.pos));
				upper = Max(upper, Point3(frame.pos));
			}
		}

		fclose(f);
	
		center /= frames.size();

		return true;
	}
	else
		return false;
}

void LoadControl(const char* path, std::vector<Control>& clicks)
{
	FILE* f = fopen(path, "r");

	if (f)
	{
		while (!feof(f))
		{
			Control c;
			int r = fscanf(f, "%f %d", &c.time, (int*)&c.event);

			if (r != 2)
				break;

			clicks.push_back(c);
		}

		fclose(f);
	}
}

vector<Frame> gFrames;
float gFrameRate = 0.01f;
float gFrameTime = 0.0f;
bool gPause = false;
bool gWireframe = false;
bool gShowNormals = true;
bool gFreeCam = true;

size_t CurrentFrame()
{
	return gFrameTime / gFrameRate;
}

const char* GetPath(const char* file)
{
#if __APPLE__

	static char path[PATH_MAX];
	uint32_t size = sizeof(path);
	_NSGetExecutablePath(path, &size);

	char* lastSlash = strrchr(path, '/');
	strcpy(lastSlash+1, file);
	return path;
#else
	return file;
#endif
}

const char* GetFile(const char* path)
{
	const char* p = strrchr(path, '\\');
	if (p)
		return p+1;
	else
		return path;
}

void Init()
{
	gFrames.clear();
	gFrameTime = 0.0f;
	
	if (!gAnimFile.empty())
	{
		const char* path = GetPath(gAnimFile.c_str());
		LoadBvh(path, gFrames, gTagCenter, gTagLower, gTagUpper);

		printf("Finished loading BVH: %s.\n", path); 
	}

	if (!gControlFile.empty())
	{
		const char* path = GetPath(gControlFile.c_str());
		LoadControl(path, gControlTrack);

		printf("Finished loading Control: %s.\n", path); 
	}
	
	gCamPos = Vec3(gTagCenter - Vec3(0.0f, 0.0f, gTagUpper.z-gTagLower.z));
	gCamAngle = Vec3(kPi, 0.0f, 0.0f);

	// create brushes
	
	gTagBrushes.push_back(new SquareBrush());
	gTagBrushes.push_back(new TriangleBrush());
	gTagBrushes.push_back(new CircleBrush());
	gTagBrushes.push_back(new SquareColorBrush(colors1));
	gTagBrushes.push_back(new SquareColorBrush(colors2));
	gTagBrushes.push_back(new SquareColorBrush(colors3));

	gMeshCan = ImportMeshFromObj(GetPath("can.obj"));	

	// translate the can so it's nozzle is at the origin
	Vec3 minExtents, maxExtents;
	gMeshCan->GetBounds(minExtents, maxExtents);
	gMeshCan->Transform(TranslationMatrix(Point3(0.0f, -(maxExtents.y-minExtents.y), 0.0f)));

}

void ToggleFullscreen()
{
	if (!gFullscreen)
	{
		glutFullScreen();
		gFullscreen = true;
	}
	else
	{
		glutReshapeWindow(1280, 720);
		gFullscreen = false;
	}
}

void DrawBasis(const Matrix44& m)
{
	glPushMatrix();
	glMultMatrixf(m);
	
	glBegin(GL_LINES);

	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3fv(Vec3(0.0f));
	glVertex3fv(Vec3(10.0f, 0.0f, 0.0f));

	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3fv(Vec3(0.0f));
	glVertex3fv(Vec3(0.0f, 10.0f, 0.0f));

	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3fv(Vec3(0.0f));
	glVertex3fv(Vec3(0.0f, 0.0f, 10.0f));

	glEnd();	
	glPopMatrix();
}

void DrawPlane(Vec3 n, float d)
{
 	const float size = 1000.0f;

	// calculate point on the plane
	Vec3 p = n*d;

	// two vectors that span the plane
	Vec3 a, b;
	BasisFromVector(n, &a, &b);
	
	a *= size;
	b *= size;

	glBegin(GL_QUADS);

	glNormal3fv(n);
	glVertex3fv(p + a + b);
	glVertex3fv(p - a + b);
	glVertex3fv(p - a - b);
	glVertex3fv(p + a - b);

	glEnd();
}

void ShadowCreate(GLuint& texture, GLuint& frameBuffer)
{
	glVerify(glGenFramebuffers(1, &frameBuffer));
	glVerify(glGenTextures(1, &texture));
	glVerify(glBindTexture(GL_TEXTURE_2D, texture));

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 
	 
	// This is to allow usage of shadow2DProj function in the shader 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_R_TO_TEXTURE); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL); 
	glTexParameteri(GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_INTENSITY); 

	glVerify(glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, 1024, 1024, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, NULL));

}

void ShadowBegin(GLuint texture, GLuint frameBuffer)
{
	glVerify(glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer));
	glVerify(glDrawBuffer(GL_NONE));
	glVerify(glReadBuffer(GL_NONE));
	glVerify(glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, texture, 0));

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, 1024, 1024);
}

void ShadowEnd()
{
	glVerify(glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0));
	glVerify(glBindFramebuffer(GL_FRAMEBUFFER, 0));
	glViewport(0, 0, gScreenWidth, gScreenHeight);
}

Tag* CreateTag(const vector<Frame>& frames, const vector<Control> controls, float smoothing, float width, float height, float velscale, uint32_t maxFrame)
{
	// re-create the tag each frame
	Tag* tag = new Tag(smoothing, width, height, velscale, gTagBrushes[gTagBrushIndex]);

	uint32_t c = 0;

	// build tag mesh
	for (size_t i=0; i < maxFrame; ++i)
	{
		// always draw if no control track 
		if (controls.empty() && i == 10)
		{
			tag->Start();
		}
		else
		{
			// check if we need to apply any control track events
			float t = i*gFrameRate;
			while (c < controls.size() && controls[c].time < t)
			{
				if (controls[c].event == eStartPaint)
				{
					tag->Start();
				}
				else
					tag->Stop();

				++c;
			}
		}

		Matrix44 m = TranslationMatrix(Point3(frames[i].pos));
		tag->PushSample(0.0f, m);
	}	
	
	// output end cap if in draw mode
	if (tag->draw)
		tag->Stop();

	return tag;
}

void DrawMesh(Mesh* m, Point3 p, GLuint shader)
{
	if (!gShowCan || gFrames.empty())
		return;

	glPushMatrix();	
		
	// use current BVH rotation frame for can
	Matrix44 x = RotationMatrix(DegToRad(gFrames[CurrentFrame()].rot.z), Vec3(0.0f, 0.0f, 1.0f))*
				 RotationMatrix(DegToRad(gFrames[CurrentFrame()].rot.x), Vec3(1.0f, 0.0f, 0.0f))*
				 RotationMatrix(DegToRad(gFrames[CurrentFrame()].rot.y), Vec3(0.0f, 1.0f, 0.0f));	

	// comment out to flip the can 180 degress back
	x *= RotationMatrix(DegToRad(180.0f), Vec3(1.0f, 0.0f, 0.0f));

	// set to anim head
	x.SetTranslation(p);
	
	glMultMatrixf(x);
	glScalef(2.0f, 2.0f, 2.0f);

	if (shader == gMainShader)
	{
		GLint uworldTransform = glGetUniformLocation(shader, "worldTransform");
		glUniformMatrix4fv(uworldTransform, 1, false, x);
	}

	glColor3f(1.5f, 1.5f, 1.5f);

	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, &m->m_positions.front());
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, 0, &m->m_normals.front());
	
	glDrawElements(GL_TRIANGLES, m->GetNumFaces()*3, GL_UNSIGNED_INT, &m->m_indices.front());

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);

	glPopMatrix();

	glColor3f(1.0f, 1.0f, 1.0f);

	// reset world transform
	if (shader == gMainShader)
	{
		GLint uworldTransform = glGetUniformLocation(shader, "worldTransform");
		glUniformMatrix4fv(uworldTransform, 1, false, Matrix44::kIdentity);
	}
}

void Advance(float dt)
{
	if (!gPause)
		gFrameTime += dt;

	size_t currentFrame = CurrentFrame();

	if (currentFrame >= gFrames.size())
	{
		gFrameTime = 0.0f;
		currentFrame = CurrentFrame();
	}

	Tag* tag = CreateTag(gFrames, gControlTrack, gTagSmoothing, gTagWidth, gTagHeight, gTagVelocityScale, currentFrame);

	glPolygonMode(GL_FRONT_AND_BACK, gWireframe?GL_LINE:GL_FILL);

	Vec3 tagExtents = gTagUpper-gTagLower;
	Point3 lightPos = gTagCenter + Vec3(0.0f, 2.0f*tagExtents.y, -1.2f*tagExtents.z);
	Point3 lightTarget = gTagCenter;

	Matrix44 lightPerspective = ProjectionMatrix(60.0f, 1.0f, 1.0f, 1000.0f);
   	Matrix44 lightView = LookAtMatrix(lightPos, lightTarget);
	Matrix44 lightTransform = lightPerspective*lightView;
	
	// shadowing pass 
	//
	ShadowBegin(gShadowTexture, gShadowFrameBuffer);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadMatrixf(lightPerspective);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixf(lightView);

	tag->Draw();

	DrawMesh(gMeshCan, tag->basis.GetTranslation(), -1);

	ShadowEnd();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	// lighting pass
	//
	glUseProgram(gMainShader);

	GLint uDiffuse = glGetUniformLocation(gMainShader, "shDiffuse");
	glUniform3fv(uDiffuse, 9, reinterpret_cast<float*>(&gShDiffuse[0].x));

	GLint uLightTransform = glGetUniformLocation(gMainShader, "lightTransform");
	glUniformMatrix4fv(uLightTransform, 1, false, lightTransform);

	GLint uLightPos = glGetUniformLocation(gMainShader, "lightPos");
	glUniform3fv(uLightPos, 1, lightPos);
	
	GLint uLightDir = glGetUniformLocation(gMainShader, "lightDir");
	glUniform3fv(uLightDir, 1, Normalize(lightTarget-lightPos));
	
	GLint uColor = glGetUniformLocation(gMainShader, "color");
	glUniform3fv(uColor, 1, Vec3(235.0f/255.0f, 244.0f/255.0f, 223.0f/255.0f));	

	GLint uworldTransform = glGetUniformLocation(gMainShader, "worldTransform");
	glUniformMatrix4fv(uworldTransform, 1, false, Matrix44::kIdentity);
	
	const Vec2 taps[] = 
	{ 
	   	Vec2(-0.326212,-0.40581),Vec2(-0.840144,-0.07358),
		Vec2(-0.695914,0.457137),Vec2(-0.203345,0.620716),
		Vec2(0.96234,-0.194983),Vec2(0.473434,-0.480026),
		Vec2(0.519456,0.767022),Vec2(0.185461,-0.893124),
		Vec2(0.507431,0.064425),Vec2(0.89642,0.412458),
		Vec2(-0.32194,-0.932615),Vec2(-0.791559,-0.59771) 
	};

	GLint uShadowTaps = glGetUniformLocation(gMainShader, "shadowTaps");
	glUniform2fv(uShadowTaps, 12, &taps[0].x);
	
	glEnable(GL_TEXTURE_2D);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, gShadowTexture);

	Point3 lower, upper;
	tag->GetBounds(lower, upper);

	lower = gTagLower;
	upper = gTagUpper;

	const Vec3 margin(0.5f, 0.1f, 0.5f);

	Vec3 edge = upper-lower;
	lower -= edge*margin;
	upper += edge*margin;

	// draw walls 
	glColor3f(1.0f, 1.0f, 1.0f);

	DrawPlane(Vec3(0.0f, 1.0f, 0.0f), lower.y);	
	DrawPlane(Vec3(0.0f, 0.0f, -1.0f), -upper.z);
	DrawPlane(Vec3(1.0f, 0.0f, 0.0f), lower.x);	
	DrawPlane(Vec3(-1.0f, 0.0f, 0.0f), -upper.x);
	//DrawPlane(Vec3(0.0f, -1.0f, 0.0f), -gTagUpper.y); 

	DrawMesh(gMeshCan, tag->basis.GetTranslation(), gMainShader);

	// draw the tag, in Mountain Dew green 
	glUniform3fv(uColor, 1, Vec3(1.0f));//Vec3(0.1f, 1.0f, 0.1f));	
	
	if (gShowNormals)
	{
		glUseProgram(gDebugShader);

		GLint uworldTransform = glGetUniformLocation(gDebugShader, "worldTransform");
		glUniformMatrix4fv(uworldTransform, 1, false, Matrix44::kIdentity);
	}

	tag->Draw();
		
	glDisable(GL_TEXTURE_2D);
	glUseProgram(0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	//DrawBasis(tag->basis);

	delete tag;
}

void DrawTimeline(int x, int y, int w, const vector<Control>& controls, float cursor, float startTime, float endTime)
{
	glDisable(GL_DEPTH_TEST);
	glBegin(GL_LINES);


	// time line
	glColor3f(1.0f, 1.0f, 1.0f);
	glVertex2i(x, y);
	glVertex2i(x+w, y);
	// terminators
	glVertex2i(x, y+5);
	glVertex2i(x, y-5);
	glVertex2i(x+w, y+5);
	glVertex2i(x+w, y-5);

	// cursor
	glColor3f(0.0f, 0.0f, 1.0f);
	int c = w*cursor;
	for (int i=-1; i<=1; ++i)
	{
		glVertex2i(x + c + i, y+5);
		glVertex2i(x + c + i, y-5);
	}

	// draw control marks
	glColor3f(1.0f, 0.0f, 0.0f);

	for (uint32_t i=0; i < controls.size(); ++i)
	{
		int p = (controls[i].time / (endTime-startTime))*w;
	
		glVertex2i(x + p, y);
	}

	glEnd();
}

void Update()
{
	const float dt = 1.0f/60.0f;
	static float t = 0.0f;
	//t += dt;

	// update camera
	const Vec3 forward(-sinf(gCamAngle.x)*cosf(gCamAngle.y), sinf(gCamAngle.y), -cosf(gCamAngle.x)*cosf(gCamAngle.y));
	const Vec3 right(Normalize(Cross(forward, Vec3(0.0f, 1.0f, 0.0f))));
	
	gCamPos += (forward*gCamVel.z + right*gCamVel.x)*dt;

	glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0f, float(gScreenWidth)/gScreenHeight, 10.0f, 10000.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	if (gFreeCam)
	{
		glRotatef(RadToDeg(-gCamAngle.x), 0.0f, 1.0f, 0.0f);
		glRotatef(RadToDeg(-gCamAngle.y), cosf(-gCamAngle.x), 0.0f, sinf(-gCamAngle.x));	
		glTranslatef(-gCamPos.x, -gCamPos.y, -gCamPos.z);
	}
	else
	{
		const float radius = 250.0f;
		const float speed = 0.15f;
		const float alpha = sinf(t*speed)*kPi*0.18f;		

		Vec3 eye = Vec3(gTagCenter) + Vec3(sinf(alpha)*radius, 50.0f, -cosf(alpha)*radius);
		gluLookAt(eye.x, eye.y, eye.z, gTagCenter.x, gTagCenter.y, gTagCenter.z, 0.0f, 1.0f, 0.0f);

		gCamPos = eye;
	}

	double startTime = GetSeconds();

	Advance(dt);
	
	double elapsedTime = GetSeconds()-startTime;

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, gScreenWidth, gScreenHeight, 0);

	if (gShowHelp)
	{
		int x = 10;
		int y = 15;
	
		char line[1024];

		glColor3f(1.0f, 1.0f, 1.0f);
		sprintf(line, "Draw Time: %.2fms", float(elapsedTime)*1000.0f);
		DrawString(x, y, line); y += 13;

		sprintf(line, "Anim Time: %.2fs", gFrameTime);
		DrawString(x, y, line); y += 13;
		sprintf(line, "Anim Frame: %d", int(CurrentFrame()));
		DrawString(x, y, line); y += 26;

		glColor3f(1.0f, 1.0f, 0.0f);
		DrawString(x, y, "Anim : %s", GetFile(gAnimFile.c_str())); y += 13;
		DrawString(x, y, "Mouse Track: %s", GetFile(gControlFile.c_str())); y += 26;
		glColor3f(1.0f, 1.0f, 1.0f);

		DrawString(x, y, "1 - Open Animation (.bvh)"); y += 13;
		DrawString(x, y, "2 - Open Mouse Track (.txt)"); y += 13;
		DrawString(x, y, "3 - Export Mesh (dew.obj)"); y += 13;
		
		if (gControlRecord)
		{
			DrawString(x, y, "4 - Stop Mouse Record  (control.txt)"); y += 26;			
		}
		else
		{
			DrawString(x, y, "4 - Start Mouse Record (control.txt)"); y += 26;
		}

		DrawString(x, y, "space - Pause/Play"); y += 13;
		DrawString(x, y, "r - Goto Start"); y += 13;
		DrawString(x, y, "f - Goto End"); y += 13;
		DrawString(x, y, "u,j - Speed (%.3f)", gFrameRate); y += 13;
		DrawString(x, y, "+/- Adjust Mouse Track Sync"); y += 26;

		DrawString(x, y, "i,k - Tag Width (%.2f)", gTagWidth); y += 13;
		DrawString(x, y, "o,l - Tag Height (%.2f)", gTagHeight); y += 13;
		DrawString(x, y, "p,; - Tag Velocity Scale (%.2f)", gTagVelocityScale); y += 13;
		DrawString(x, y, "/   - Next Brush"); y += 26;

		DrawString(x, y, "h - Hide Help"); y += 13;
		DrawString(x, y, "c - Hide Can"); y += 13;
		DrawString(x, y, "b - Show Wireframe", gTagWidth); y += 13;
		DrawString(x, y, "n - Show Normals", gTagHeight); y += 13;
		DrawString(x, y, "b - Toggle Fullscreen"); y += 13;
		

		if (gControlRecord)
		{
			glColor3f(1.0f, 0.0f, 0.0f);
			DrawString(x, y, "Recording"); 
			glColor3f(1.0f, 1.0f, 1.0f);
		}

		if (gPause)
		{
			glColor3f(1.0f, 0.0f, 0.0f);
			DrawString(x, y, "Paused"); 
			glColor3f(1.0f, 1.0f, 1.0f);
		}

		// time lines
		y = gScreenHeight-60;
		x = 60;
		//DrawTimeline(x, y, gScreenWidth-x*2, vector<Control>(), CurrentFrame()/float(gFrames.size()), 0.0f, gFrames.size()*gFrameRate); y += 20;
		DrawTimeline(x, y, gScreenWidth-x*2, gControlTrack, CurrentFrame()/float(gFrames.size()), 0.0f, gFrames.size()*gFrameRate); y += 20;

		if (gControlRecord)
		{
			glColor3f(1.0f, 0.0f, 0.0f);
			DrawString(x, y, "recording");
		}
	}


	glutSwapBuffers();

	/* enable to dump frames for video
	
	static int i=0;
	char buffer[255];
	sprintf(buffer, "dump/frame%d.tga", ++i);
		
	TgaImage img;
	img.m_width = gScreenWidth;
	img.m_height = gScreenHeight;
	img.m_data = new uint32_t[gScreenWidth*gScreenHeight];
		
	glReadPixels(0, 0, gScreenWidth, gScreenHeight, GL_RGBA, GL_UNSIGNED_BYTE, img.m_data);
		
	TgaSave(buffer, img);
		
	delete[] img.m_data;
	*/
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
		
			if (gControlRecord && b == GLUT_LEFT_BUTTON)
			{	
				Control c;
				c.time = GetSeconds()-gControlStartTime;
			    c.event = eStopPaint;

				gControlTrack.push_back(c);
			}
			break;
		}	
		case GLUT_DOWN:
		{
			lastx = x;
			lasty = y;

			if (gControlRecord && b == GLUT_LEFT_BUTTON)
			{	
				Control c;
				c.time = GetSeconds()-gControlStartTime;
			    c.event = eStartPaint;	

				gControlTrack.push_back(c);
			}
			break;
		}
	};
}

void GLUTKeyboardDown(unsigned char key, int x, int y)
{
	const float kSpeed = 120.0f;

	switch (key)
	{
		case 'w':
		{
			gCamVel.z = kSpeed;
			break;
		}
		case 's':
		{
			gCamVel.z = -kSpeed;
			break;
		}
		case 'a':
		{
			gCamVel.x = -kSpeed;
			break;
		}
		case 'd':
		{
			gCamVel.x = kSpeed;
			break;
		}
		case 'u':
		{
			gFrameRate += 0.001f;
			break;
		}
		case 'j':
		{
			gFrameRate -= 0.001f;
			gFrameRate = max(0.001f, gFrameRate);
			break;
		}
		case 'h':
		{
			gShowHelp = !gShowHelp;
			break;
		}
		case 'c':
		{
			gShowCan = !gShowCan;
			break;
		}
		case 'r':
		{
			gFrameTime = 0.0f;
			gPause = true;
			break;
		}
		case 'f':
		{
			gFrameTime = gFrames.size()*gFrameRate;
			gPause = true;
			break;
		}
		case ' ':
		{
			gPause = !gPause;
			break;
		}
		case 'b':
		{
			gWireframe = !gWireframe;
			break;
		}
		case 'n':
		{
			gShowNormals = !gShowNormals;
			break;
		}
		case 'i':
		{
			gTagWidth += 0.1f;
			break;
		}
		case 'k':
		{
			gTagWidth -= 0.1f;
			break;
		}
		case 'o':
		{
			gTagHeight += 0.1f;
			break;
		}
		case 'l':
		{
			gTagHeight -= 0.1f;
			break;
		}
		case 'p':
		{
			gTagVelocityScale += 0.01f;
			break;
		}
		case ';':
		{
			gTagVelocityScale -= 0.01f;
			break;
		}
		case '=':
		{
			for (uint32_t i=0; i < gControlTrack.size(); ++i)
				gControlTrack[i].time += 0.01f;
			break;
		}
		case '-':
		{
			for (uint32_t i=0; i < gControlTrack.size(); ++i)
				gControlTrack[i].time -= 0.01f;
			break;
		}
		case 'm':
		{
			ToggleFullscreen();
			break;
		}
		case '/':
		{
			gTagBrushIndex = (gTagBrushIndex+1)%gTagBrushes.size();
			break;
		}


		/*
		case 'g':
		{
			gFreeCam = !gFreeCam;
			break;
		}
		*/
		case '1':
		{
#ifndef _WIN32
			if (gFullscreen)
				break;	
#endif
			string path = FileOpenDialog();

			if (!path.empty())
			{
				gAnimFile = path;

				// re-init
				Init();
			}
//#endif
			break;
		}

		case '2':
		{
#ifndef _WIN32
			if (gFullscreen)
				break;	
#endif

			string path = FileOpenDialog();

			if (!path.empty())
			{
				gControlFile = path;

				// re-init
				Init();
			}
			break;
		}
		case '3':
		{
			// create tag copy
			Tag* tag = CreateTag(gFrames, gControlTrack, gTagSmoothing, gTagWidth, gTagHeight, gTagVelocityScale, CurrentFrame());
			
			// export
			tag->ExportToObj("dew.obj");
			delete tag;

			printf("Exported dew.obj mesh\n");

			break;
		}
		case '4':
		{			
			if (!gControlRecord)
			{
				gFrameTime = 0.0f;

				gControlTrack.clear();
				gControlRecord = true;
				gControlStartTime = GetSeconds();
			}
			else
			{
				const char* path = GetPath("control.txt");

				FILE* f = fopen(path, "w");
		
				if (f)
				{	
					for (uint32_t i=0; i < gControlTrack.size(); ++i)
					{
						fprintf(f, "%f %i\n", gControlTrack[i].time, gControlTrack[i].event);	
					}

					fclose(f);
				}

				gControlRecord = false;				
				gControlFile = "control.txt";
			}
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
			gCamVel.z = 0.0f;
			break;
		}
		case 's':
		{
			gCamVel.z = 0.0f;
			break;
		}
		case 'a':
		{
			gCamVel.x = 0.0f;
			break;
		}
		case 'd':
		{
			gCamVel.x = 0.0f;
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

	if (!gControlRecord)
	{
		const float kSensitivity = DegToRad(0.1f);

		gCamAngle.x -= dx*kSensitivity;
		gCamAngle.y += dy*kSensitivity;
	}
}

void GLUTReshape(int x, int y)
{
	gScreenWidth = x;
	gScreenHeight = y;
}

int main(int argc, char* argv[])
{
	if (argc < 2)
		printf("BVH file not specified, defaulting to %s.\n", gAnimFile.c_str());
	else
	{
		gAnimFile = argv[1];

		if (argc > 2)
			gControlFile = argv[2];
	}

	// init gl
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
	
	glutInitWindowSize(gScreenWidth, gScreenHeight);
	glutCreateWindow("Dew");
	glutPositionWindow(200, 100);

#if _WIN32
	glewInit();	
#endif

	gMainShader = CompileProgram(vertexShader, fragmentShaderMain);
	if (gMainShader == 0)
		return 0;

	gDebugShader = CompileProgram(vertexShader, fragmentShaderDebug);
	if (gDebugShader == 0)
		return 0;

	ShadowCreate(gShadowTexture, gShadowFrameBuffer);

	Init();

	glutIdleFunc(Update);	
	glutDisplayFunc(Update);
	glutMouseFunc(GLUTMouseFunc);
	glutMotionFunc(GLUTMotionFunc);
	glutKeyboardFunc(GLUTKeyboardDown);
	glutKeyboardUpFunc(GLUTKeyboardUp);
	glutReshapeFunc(GLUTReshape);

#if __APPLE__
	// enable v-sync
	int swap_interval = 1;
	CGLContextObj cgl_context = CGLGetCurrentContext();
	CGLSetParameter(cgl_context, kCGLCPSwapInterval, &swap_interval);
#endif

	glutMainLoop();
	return 0;
}





