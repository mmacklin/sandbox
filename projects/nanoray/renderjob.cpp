#include "renderjob.h"

#include "integrators.h"

#include <iostream>

using namespace std;

extern _declspec(thread) MemoryArena* g_memArena;

uint32_t RenderTileThreadFunc(void* data)
{
	RenderJob& job = *(RenderJob*)(data);

	//TODO: this should just be a member
	SampleBuffer out(job.m_output, job.m_renderRect);

	const Point3 rayOrigin = job.m_camera.m_cameraToWorld.GetTranslation();
	const Rect rect = job.m_renderRect;

	// take image samples
	for (uint32_t i=rect.Top(); i < rect.Bottom(); ++i)
	{
		for (uint32_t j=rect.Left(); j < rect.Right(); ++j)
		{
			Colour radiance(0.0f);

			for (uint32_t s=0; s < job.m_samplesPerPixel; ++s)
			{
				Vector3 rayDir;
				job.m_camera.GenerateRay(j, i, rayDir);			

#if _DEBUG
				if (j == 22 && i == 16 && s == 48)
				{
					cout << "Blah" << endl;
				}
#endif
				Colour L = PathTrace(*job.m_scene, rayOrigin, rayDir);
				//Colour L = ForwardTraceImportance(*m_scene, rayOrigin, rayDir);
				//Colour L = Debug(*m_scene, rayOrigin, rayDir);

#if _DEBUG
				if (!_finite(L.r) || _isnan(L.r) || 
					!_finite(L.g) || _isnan(L.g) ||
					!_finite(L.b) || _isnan(L.b))
				{
					cout << "Error" << endl;
				}
#endif
				radiance += L;

				// reset the memory arena
				g_memArena->Reset();
			}

			out.AddSample(float(j), float(i), radiance * (1.0f / job.m_samplesPerPixel));
		}
	}
	
	return 0;
}