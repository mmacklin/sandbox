#pragma once

#include "types.h"			

#include <cassert>

inline void* AlignedMalloc(size_t size,int byteAlign)
{
    void *mallocPtr = malloc(size + byteAlign + sizeof(void*));
    size_t ptrInt = (size_t)mallocPtr;

    ptrInt = (ptrInt + byteAlign + sizeof(void*)) / byteAlign * byteAlign;
    *(((void**)ptrInt) - 1) = mallocPtr;

    return (void*)ptrInt;
}

inline void AlignedFree(void *ptr)
{
    free(*(((void**)ptr) - 1));
}

// temporary allocation storage, used to reduce load on crt malloc
class MemoryArena
{
public:

	MemoryArena(const uint32_t sizeInBytes)
	{
		m_mem = (uint8_t*)AlignedMalloc(sizeInBytes, 16);
		m_size = sizeInBytes;
		m_head = m_mem;
	}

	~MemoryArena()
	{
		AlignedFree(m_mem);
	}

	uint8_t* Allocate(uint32_t size)
	{
		if ((m_head+size)-m_mem > ptrdiff_t(m_size))
		{
			assert(!"Arena ran out of memory");
			return NULL;
		}

		uint8_t* p = m_head;
		m_head += size;

		return p;
	}
	
	void Reset()
	{
		m_head = m_mem;
	}

	uint8_t* m_mem;
	uint8_t* m_head;
	uint32_t m_size;
};

// version of placement new to support usage as: new (MemArena) Type(); 
inline void* operator new (size_t s, MemoryArena& a)
{
	return a.Allocate(s);
}

// not used but compiler will complain without it
inline void operator delete(void*, MemoryArena&)
{
}
