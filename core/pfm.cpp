#include "pfm.h"

#include <cassert>
#include <stdio.h>
#include <string.h>
#include <algorithm>

namespace
{
	// RAII wrapper to handle file pointer clean up
	struct FilePointer
	{
		FilePointer(FILE* ptr) : p(ptr) {}
		~FilePointer() { if (p) fclose(p); }

		operator FILE*() { return p; }

		FILE* p;
	};
}

bool PfmLoad(const char* filename, PfmImage& image)
{
	FilePointer f = fopen(filename, "rb");
	if (!f)
		return false;

	memset(&image, 0, sizeof(PfmImage));
	
	const uint32_t kBufSize = 1024;
	char buffer[kBufSize];
	
	if (!fgets(buffer, kBufSize, f))
		return false;
	
	if (strcmp(buffer, "PF\n") != 0)
		return false;
	
	if (!fgets(buffer, kBufSize, f))
		return false;

	image.m_depth = 1;
	sscanf(buffer, "%d %d %d", &image.m_width, &image.m_height, &image.m_depth);

	if (!fgets(buffer, kBufSize, f))
		return false;
	
	sscanf(buffer, "%f", &image.m_maxDepth);
	
	uint32_t dataStart = ftell(f);
	fseek(f, 0, SEEK_END);
	uint32_t dataEnd = ftell(f);
	fseek(f, dataStart, SEEK_SET);
	
	uint32_t dataSize = dataEnd-dataStart;

	// must be 4 byte aligned
	assert((dataSize&0x3) == 0);
	
	image.m_data = new float[dataSize/4];
	
	if (fread(image.m_data, dataSize, 1, f) != 1)
		return false;
	
	return true;
}

void PfmSave(const char* filename, const PfmImage& image)
{
	FILE* f = fopen(filename, "wb");
	if (!f)
		return;

	fprintf(f, "PF\n");
	if (image.m_depth > 1)
		fprintf(f, "%d %d %d\n", image.m_width, image.m_height, image.m_depth);
	else
		fprintf(f, "%d %d\n", image.m_width, image.m_height);

	fprintf(f, "%f\n", *std::max_element(image.m_data, image.m_data+(image.m_width*image.m_height*image.m_depth)));

	fwrite(image.m_data, image.m_width*image.m_height*image.m_depth*sizeof(float), 1, f);
}





