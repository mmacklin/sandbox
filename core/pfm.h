#include "types.h"

struct PfmImage
{
	// set m_depth to 1 for standard pfm compatability, > 1 will act as a volume texture (non-standard)
	uint32_t m_width;
	uint32_t m_height;
	uint32_t m_depth;

	// optional
	float m_maxDepth;
	
	float* m_data;
};

bool PfmLoad(const char* filename, PfmImage& image);
void PfmSave(const char* filename, const PfmImage& image);
