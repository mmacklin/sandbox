#include "core.h"

void radix_sort(uint32_t* begin, uint32_t* end, uint32_t* auxBegin, uint32_t* auxEnd)
{
	static uint32_t tables[2][1 << 16];
	memset(tables, 0, sizeof(tables));
		
	// build histograms
	for (Iter iter=begin; iter != end; ++iter)
	{
		const uint16_t low = (*iter) & 0xffff;
		const uint16_t high = (*iter) >> 16;
		
		++tables[0][low];
		++tables[1][high];
	}
	
	// convert histograms to offset tables in-place
	uint32_t offlow = 0;
	uint32_t offhigh = 0;
	
	for (uint32_t i=0; i < 65536; ++i)
	{
		const uint32_t newofflow = offlow + tables[0][i];
		const uint32_t newoffhigh = offhigh + tables[1][i];
		
		tables[0][i] = offlow;
		tables[1][i] = offhigh;
		
		offlow = newofflow;
		offhigh = newoffhigh;
	}
		
	// pass 1 - sort by low 16 bits
	for (uint32_t* iter=begin; iter != end; ++iter)
	{
		// lookup offset of input
		const uint32_t x = *iter;
		const uint32_t i = x & 0xffff;
		
		// find offset and increment
		const uint32_t offset = tables[0][i]++;
		
		auxBegin[offset] = x;
	}	
		
	// pass 2 - sort by high 16 bits
	for (uint32_t* iter=auxBegin; iter != auxEnd; ++iter)
	{
		// lookup offset of input
		const uint32_t x = *iter;
		const uint32_t i = x >> 16;
		
		const uint32_t offset = tables[1][i]++;
		
		begin[offset] = x;
	}	
}
