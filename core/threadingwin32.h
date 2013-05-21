#pragma once

class Thread
{
public:

	Thread(ThreadFunc func, void* param)
	{
		m_handle = CreateThread(NULL, 0, func, param, 0, NULL);
	}

	HANDLE m_handle;
}

void WaitForThreads(Thread* threads, uint32_t n)
{

}