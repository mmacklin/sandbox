#include "threading.h"

#include <cassert>

#if USE_PTHREADS

CriticalSection::CriticalSection()
{
	int err = pthread_mutex_init(&m_semaphore, NULL);
	
	assert(err == 0);
}

CriticalSection::~CriticalSection()
{
	pthread_mutex_destroy(&m_semaphore);
}

void CriticalSection::Enter()
{
	pthread_mutex_lock(&m_semaphore);
}

void CriticalSection::Leave()
{
	pthread_mutex_unlock(&m_semaphore);
}

#else

CriticalSection::CriticalSection()
{
	InitializeCriticalSection(&m_semaphore);
}

CriticalSection::~CriticalSection()
{
	DeleteCriticalSection(&m_semaphore);
}

void CriticalSection::Enter()
{
	EnterCriticalSection(&m_semaphore);
}

void CriticalSection::Leave()
{
	LeaveCriticalSection(&m_semaphore);
}

#endif
