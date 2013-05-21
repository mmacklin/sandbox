#pragma once

#include "core/core.h"
#include "core/maths.h"

// 2d and 3d signed distance field computation using fast marching method (FMM), output array
// should be the same size as input, non-zero input pixels will have distance < 0.0f, resulting 
// distance is scaled by 1 / max(dimension)
void MakeSDF(const uint32_t* input, uint32_t width, uint32_t height, float* output);
void MakeSDF(const uint32_t* input, uint32_t width, uint32_t height, uint32_t depth, float* output);
