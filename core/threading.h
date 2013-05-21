#pragma once

#if _WIN32
#define NOMINMAX
#include <windows.h>
typedef CRITICAL_SECTION CriticalSectionHandle;
#else
#define USE_PTHREADS 1
#include <pthread.h>
typedef pthread_mutex_t CriticalSectionHandle;
#endif

class CriticalSection
{
public:

	CriticalSection();
	~CriticalSection();

	void Enter();
	void Leave();

private:

	CriticalSectionHandle m_semaphore;
};



// RAII style critical section lock
class CriticalSectionScopeLock
{
	CriticalSectionScopeLock& operator=(const CriticalSectionScopeLock&);

public:
	
	CriticalSectionScopeLock(CriticalSection& cs) : m_cs(cs)
	{
		m_cs.Enter();
	}

	~CriticalSectionScopeLock()
	{
		m_cs.Leave();
	}

	CriticalSection& m_cs;
};	

