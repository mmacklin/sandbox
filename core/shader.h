#pragma once

#if _WIN32

#include <external/glew/include/gl/glew.h>

#if USE_FREEGLUT
#include <external/freeglut/include/GL/freeglut.h>
#else
#include <external/glut/glut.h>
#endif

#elif __APPLE__
#define GL_DO_NOT_WARN_IF_MULTI_GL_VERSION_HEADERS_INCLUDED 
#include <opengl/gl3.h>
#include <glut/glut.h>
#elif PLATFORM_IOS

#if OGL1
#import <OpenGLES/EAGL.h>
#import <OpenGLES/ES1/gl.h>
#import <OpenGLES/ES1/glext.h>
#else
#import <OpenGLES/ES2/gl.h>
#import <OpenGLES/ES2/glext.h>
#endif

#endif

#include "core/maths.h"

#define glVerify(x) {x; glAssert(#x, __LINE__, __FILE__);}
void glAssert(const char* msg, long line, const char* file);

GLuint CompileProgramFromFile(const char *vertexPath, const char *fragmentPath);
GLuint CompileProgram(const char *vsource=NULL, const char *fsource=NULL, const char* gsource=NULL);
GLuint CompileProgram(const char *vsource, const char* csource, const char* esource, const char* fsource);

void DrawPlane(const Vec4& p, bool color=true);
void DrawString(int x, int y, const char* s, ...);
void DrawFrustum(const Matrix44& projToWorld);
