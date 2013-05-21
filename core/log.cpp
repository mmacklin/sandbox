#include "log.h"

#include "threading.h"
#include "platform.h"

#include <algorithm>
#include <iostream>
#include <cassert>
#include <string>

using namespace std;

Log Log::Info;
Log Log::Warn;
Log Log::Error;

void Log::RegisterListener(LogListener* l)
{
	std::vector<LogListener*>::iterator iter = std::find(mListeners.begin(), mListeners.end(), l);

	if (iter == mListeners.end())
		mListeners.push_back(l);

}

void Log::UnRegisterListener(LogListener* l)
{
	std::vector<LogListener*>::iterator iter = std::find(mListeners.begin(), mListeners.end(), l);

	if (iter != mListeners.end())
		mListeners.erase(iter);
}

void Log::Flush(void)
{
	// get string from string stream
	string d = m_stream.str();
	
	// reset stream
	m_stream.str("");

	// send to all listeners
	for (size_t i=0; i < mListeners.size(); i++)
		mListeners[i]->Output(d);

	// default output
	//std::cout << d.c_str();
}



// log to file
class LogFileListener : public LogListener
{
public:

	LogFileListener(const char* filename)
	{
		m_file = fopen(filename, "w");
		//std::cout << "temp: " << GetTempDirectory() << std::endl;
		std::cout << filename << std::endl;
		assert(m_file);
	}

	virtual ~LogFileListener()
	{
		fclose(m_file);
	}

	virtual void Output(const string& msg)
	{	
		CriticalSectionScopeLock lock(m_mutex);

		fwrite(msg.c_str(), msg.length(), 1, m_file);
	}

	CriticalSection m_mutex;
	FILE* m_file;
};

#if WIN32

#include <windows.h>


// outputs to debugger console
class LogDebugListener : public LogListener
{
public:

	virtual void Output(const string& msg)
	{
		OutputDebugString(msg.c_str());
		std::cout << msg.c_str();
	}
};

#else
// outputs to debugger console
class LogDebugListener : public LogListener
{
public:
	
	virtual void Output(const string& msg)
	{
		std::cout << msg.c_str();
	}
};

#endif

namespace
{
	LogFileListener* gLogFileListener;
	LogDebugListener* gDebugListener;
}

// initialize default log listeners
void Log::Init(const char* logFile)
{
	gLogFileListener = new LogFileListener(logFile);
	gDebugListener = new LogDebugListener();
										   
	Log::Info.RegisterListener(gLogFileListener);
	Log::Warn.RegisterListener(gLogFileListener);
	Log::Error.RegisterListener(gLogFileListener);

	Log::Info.RegisterListener(gDebugListener);
	Log::Warn.RegisterListener(gDebugListener);
	Log::Error.RegisterListener(gDebugListener);
}

void Log::Shutdown()
{
	Log::Info.UnRegisterListener(gLogFileListener);
	Log::Warn.UnRegisterListener(gLogFileListener);
	Log::Error.UnRegisterListener(gLogFileListener);

	Log::Info.UnRegisterListener(gDebugListener);
	Log::Warn.UnRegisterListener(gDebugListener);
	Log::Error.UnRegisterListener(gDebugListener);

	delete gLogFileListener;
	delete gDebugListener;	
}
