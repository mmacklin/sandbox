#include "integrators.h"

float InScatter(const Point3& start, const Vector3& dir, const Point3& lightPos, float t)
{
	// calculate quadratic coefficients a,b,c
	Vector3 q = start - lightPos;

	float a = 1.0f;
	float b = 2.0f*Dot3(dir, q);
	float c = Dot3(q, q);

	// evaluate integral
	float s = 1.0f / sqrtf(c - b*b*0.25f);

	float l = s * (atanf( (t + b*0.5f) * s) - atanf( b*0.5f*s ));

	return l;	
}

extern _declspec(thread) MemoryArena* g_memArena;

Colour PathTrace(const Scene& scene, const Point3& startOrigin, const Vector3& startDir)
{	
	// path throughput
	Colour pathThroughput(1.0f, 1.0f, 1.0f, 1.0f);
	// accumulated radiance along the path
	Colour totalRadiance(0.0f);

	Point3 rayOrigin = startOrigin;
	Vector3 rayDir = startDir;

	float t = 0.0f;
	Vector3 n;
	const Primitive* hit;

	const uint32_t kMaxPathDepth = 4;

    g_memArena->Reset();

	for (uint32_t i=0; i < kMaxPathDepth; ++i)
	{
		// find closest hit
		if (scene.Trace(rayOrigin, rayDir, t, n, &hit))
		{			
			// update position and path direction
			Point3 p = rayOrigin + rayDir*t;

			// first trace is our only chance to add contribution from directly visible light sources
            if (i == 0)
			{
				totalRadiance += hit->m_emission;
			}

			const BRDF* brdf = hit->m_material->GetBRDF(p, n);

			// sample light sources
			totalRadiance += pathThroughput * scene.SampleLights(p, n, -rayDir, brdf);
			
			// generate new path direction by sampling BRDF
			float pdf = 1.f;
			Vector3 wi;
			brdf->Sample(-rayDir, wi, pdf);
		
			// evaluate brdf
			Colour f = brdf->F(-rayDir, wi);

			// update path throughput 
			pathThroughput *= f * Abs(Dot(n, wi)) / pdf;

			// update path start and direction
			rayOrigin = p;
			rayDir = wi;

			//delete brdf;
		}
            
		else
		{
			/*
			// sample volumes
			Colour volumeEmission;
			Colour volumeTransmission;
			scene.SampleVolumes(rayOrigin, rayDir, t, volumeTransmission, volumeEmission); 

			pathThroughput *= volumeTransmission;
			totalRadiance += pathThroughput * (volumeEmission);
			*/
			// hit nothing, evaluate background li and end loop
			totalRadiance += pathThroughput * scene.SampleSky(rayDir);						
			
			break;
		}
	}

	return totalRadiance;
}

// reference, no light sampling, uniform hemisphere sampling
Colour ForwardTraceUniform(const Scene& scene, const Point3& startOrigin, const Vector3& startDir)
{	
    // path throughput
    Colour pathThroughput(1.0f, 1.0f, 1.0f, 1.0f);
    // accumulated radiance
    Colour totalRadiance(0.0f);

    Point3 rayOrigin = startOrigin;
    Vector3 rayDir = startDir;

    float t = 0.0f;
    Vector3 n(rayDir);
    const Primitive* hit;

    const uint32_t kMaxPathDepth = 8;
    float pdf = 1.0f;

    g_memArena->Reset();

    for (uint32_t i=0; i < kMaxPathDepth; ++i)
    {
        // find closest hit
        if (scene.Trace(rayOrigin, rayDir, t, n, &hit))
        {			
            // calculate a basis for this hit point
            const Matrix44 localToWorld = TransformFromVector(n);

            const Point3 p = rayOrigin + rayDir*t;

            // update position and path direction
            const Vector3 outDir = localToWorld*UniformSampleHemisphere();

            // update total radiance
            totalRadiance += hit->m_emission * pathThroughput;

            // reflectance
            BRDF* brdf = hit->m_material->GetBRDF(p, n);
            Colour f = brdf->F(-rayDir, outDir);

            // update throughput with primitive reflectance
            pathThroughput *= f * Clamp(Dot(n, outDir), 0.0f, 1.0f) / kInv2Pi;

            rayDir = outDir;
            rayOrigin = p;
        }
        else
        {
            // hit nothing, terminate loop
            break;
        }
    }

    return totalRadiance;
}

// reference, no light sampling but does cosine weighted sampling
Colour ForwardTraceImportance(const Scene& scene, const Point3& startOrigin, const Vector3& startDir)
{	
	// path throughput
	Colour pathThroughput(1.0f, 1.0f, 1.0f, 1.0f);
	// accumulated radiance
	Colour totalRadiance(0.0f);

	Point3 rayOrigin = startOrigin;
	Vector3 rayDir = startDir;

	float t = 0.0f;
	Vector3 n(rayDir);
	const Primitive* hit;

	const uint32_t kMaxPathDepth = 8;

	for (uint32_t i=0; i < kMaxPathDepth; ++i)
	{
		// find closest hit
		if (scene.Trace(rayOrigin, rayDir, t, n, &hit))
		{			
				// update position and path direction
			Point3 p = rayOrigin + rayDir*t;

			totalRadiance += pathThroughput * hit->m_emission;

			const BRDF* brdf = hit->m_material->GetBRDF(p, n);

				// generate new path direction by sampling BRDF
			float pdf = 1.f;
			Vector3 wi;
			brdf->Sample(-rayDir, wi, pdf);
		
			// evaluate brdf
			Colour f = brdf->F(-rayDir, wi);

			// update path throughput 
			pathThroughput *= f * Abs(Dot(n, wi)) / pdf;

			// update path start and direction
			rayOrigin = p;
			rayDir = wi;

			delete brdf;
		}
		else
		{
			// hit nothing, terminate loop
			break;
		}
	}

	return totalRadiance;
}


/*

Colour Whitted(const Scene& s, const Point3& rayOrigin, const Vector3& rayDir)
{
	// TODO:

	return Colour();
}


*/

Colour Debug(const Scene& scene, const Point3& rayOrigin, const Vector3& rayDir)
{
	// find closest hit
	float t;
	Vector3 n;
	const Primitive* p;
	if (scene.Trace(rayOrigin, rayDir, t, n, &p))
	{
		return Colour(0.5f*n.x+0.5f, 0.5f*n.y+0.5f, 0.5f*n.z+0.5f, 1.0);
	}

	return Colour(0.0f);
}
