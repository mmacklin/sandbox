#pragma once

#include <pthread.h>

class Thread
{
public:

	pthread_t m_handle;

	Thread(ThreadFunc func, void* param)
	{
		int retval = pthread_create(&m_handle, NULL, func, param);
		assert(retval);		
	}
};

void WaitForThreads(const Thread* threads, uint32_t n)
{
	// just join all the threads
	for (uint32_t i=0; i < n; ++i)
	{
		pthread_join(threads[i], NULL);
	}
}
