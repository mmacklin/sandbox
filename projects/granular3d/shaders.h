#pragma once

#define STRINGIFY(A) #A

#include <core/shader.h>
#include <core/maths.h>

void ShadowCreate(GLuint& texture, GLuint& frameBuffer);
void ShadowBegin(GLuint texture, GLuint frameBuffer);
void ShadowEnd();

void DrawPlanes(Vec4* planes, int n, Vec3 lightPos, Vec3 lightTarget, Matrix44 lightTransform, GLuint shadowTex);
void DrawPoints(float* positions, int n, float radius, float screenWidth, float screenAspect, float fov, Vec3 lightPos, Vec3 lightTarget, Matrix44 lightTransform, GLuint shadowTex);

