#pragma once

#include "core/maths.h"
#include "core/tga.h"

#include <vector>

void TriangulateDelaunay(const Vec2* points, uint32_t numPoints, std::vector<Vec2>& outPoints, std::vector<uint32_t>& outTris);
void TriangulateVariational(const Vec2* points, uint32_t numPoints, const Vec2* bpoints, uint32_t numBPoints, uint32_t iterations, std::vector<Vec2>& outPoints, std::vector<uint32_t>& outTris);

void CreateDonut(std::vector<Vec2>& points, std::vector<uint32_t>& indices, float inner, float outer, uint32_t segments);




