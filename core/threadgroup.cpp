#include "threadgroup.h"

#include "platform.h"

#include <iostream>
#include <cassert>

#include "memory.h"

#if defined(_WIN32)

using namespace std;

// these two things should not be here
_declspec(thread) MemoryArena* g_memArena;

DWORD WINAPI ThreadPool::WorkerMain(void* data)
{

	if (data == NULL)
		return 0;

	ThreadPool* pool = (ThreadPool*)(data);

	//double start = GetSeconds();
	//uint32_t tasksExecuted = 0;
	
    // initialize tls
    g_memArena = new MemoryArena(128*1024);

    // pull task from thread pool
    for(;;)
    {
        Task t;
        t.m_func = NULL;

		{
			CriticalSectionScopeLock cs(pool->m_mutex);
			
			if (!pool->m_tasks.empty())
			{			
				t = pool->m_tasks.back();
				pool->m_tasks.pop_back();
			}
			else
			{
				//double totalTime = GetSeconds() - start;
				//cout << "Total worker time: " << totalTime << " Tasks: " << tasksExecuted << endl;

				// no tasks so return
				delete g_memArena;

				return 0;
			}
		}

        // execute job
        if (t.m_func)
        {
            t.m_func(t.m_param);
			
			//++tasksExecuted;
        }
    }

	delete g_memArena;

	return 0;
}

ThreadPool::ThreadPool()
{
}

ThreadPool::~ThreadPool()
{
    for (uint32_t i=0; i < m_threads.size(); ++i)
    {
        CloseHandle(m_threads[i]);
    }
}

void ThreadPool::AddTask(ThreadFunc func, void* param)
{
    CriticalSectionScopeLock cs(m_mutex);

    Task t;
    t.m_func = func;
    t.m_param = param;

    m_tasks.push_back(t);
}

#include "maths.h"

void ThreadPool::Run(uint32_t workerCount)
{
	// create workers
    for(uint32_t i=0; i < workerCount; ++i)
    {
        HANDLE newThread = CreateThread(NULL, 0, WorkerMain, this, 0, NULL);
        if (newThread == NULL)
		{
			cout << "CreateThread failed with: " << GetLastError() << endl;
			return;
		}

        m_threads.push_back(newThread);
    }
}

void ThreadPool::Wait()
{
	WaitForMultipleObjects(m_threads.size(), &m_threads[0], TRUE, INFINITE);
}


#if USE_PTHREADS

void ThreadPool::AddWorker(ThreadFunc func, void* param)
{
    pthread_t newThread;

    int retval = pthread_create(&newThread, NULL, WorkerMain, this);

    if (retval)
    {
        cerr << "Could not create pthread, error code: " << retval << endl;
        return;
    }

    m_threads.push_back(newThread);
}	

void ThreadPool::Wait()
{
    for (uint32_t i=0; i < m_threads.size(); ++i)
    {
        pthread_join(m_threads[i], NULL);
    }		
}

#endif // USE_PTHREADS

#endif // _WIN32
