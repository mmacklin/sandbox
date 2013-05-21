#pragma once

#include "core/core.h"
#include "core/maths.h"

#include "scene.h"

Colour PathTrace(const Scene& s, const Point3& rayOrigin, const Vector3& rayDir);
Colour ForwardTraceImportance(const Scene& scene, const Point3& startOrigin, const Vector3& startDir);
Colour ForwardTraceUniform(const Scene& scene, const Point3& startOrigin, const Vector3& startDir);
Colour Whitted(const Scene& s, const Point3& rayOrigin, const Vector3& rayDir);
Colour Debug(const Scene& s, const Point3& rayOrigin, const Vector3& rayDir);

