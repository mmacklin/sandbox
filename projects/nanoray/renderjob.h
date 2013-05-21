#pragma once

#include "scene.h"
#include "camera.h"


// holds samples
class SampleBuffer
{
public:

	SampleBuffer(Colour* memory, const Rect& r) : m_data(memory), m_rect(r)
	{

	}

	// accumulates samples specified in continuous raster space
	void AddSample(float x, float y, const Colour& c)
	{
		assert(m_rect.Contains(uint32_t(x), uint32_t(y)));

		// calculate offset into data
		const uint32_t row = uint32_t(y) - m_rect.Top();
		const uint32_t col = uint32_t(x) - m_rect.Left();
		
		const uint32_t i = (m_rect.Width() * row) + col;
		
		m_data[i] = c;	
	}
	
	Rect m_rect;
	Colour* m_data;

} ;


// perform rendering
struct RenderJob
{

	// raster space sub rect to process
	Rect m_renderRect;
	// raster space rect to filter over
	Rect m_filterRect;

	// our view
	Camera m_camera;

	// scene description
	const Scene* m_scene;

	// output, executing process must supply filtered radiance values
	Colour* m_output;
	uint32_t  m_outputPitch;

	// quality settings
	uint32_t m_samplesPerPixel;

};



