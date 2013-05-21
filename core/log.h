#pragma once

#include <vector>
#include <sstream>
#include <iomanip>

#include "types.h"

class LogListener;

class Log
{
public:

	Log() {};
	virtual ~Log() {}

	//! register stream
	void RegisterListener(LogListener* l);
	void UnRegisterListener(LogListener* l);

	// handles std::endl which is used to mark termination of log messages
	//std::stringstream& operator<<(std::basic_ostream<char, std::char_traits<char> >& (*_Pfn)(std::basic_ostream<char, std::char_traits<char> >&));
	std::stringstream& operator<<(std::basic_ostream<char, std::char_traits<char> >& (*_Pfn)(std::basic_ostream<char, std::char_traits<char> >&))
	{
		_Pfn(m_stream);	
		Flush();
		
		return m_stream;
	}
	
	// logs any old thing it can convert to a string
	template <typename T>
	Log& operator << (const T& v)
	{
		m_stream << v;
		return *this;
	}

	// create default loggers
	static void Init(const char* logFile="log.txt");
	static void Shutdown();

	// log sinks
	static Log Info;
	static Log Warn;
	static Log Error;

private:

	void Flush();

	// stringstream used for temp conversion
	std::stringstream m_stream;

	typedef std::vector<LogListener*> ListenerArray;
	ListenerArray mListeners;

};

// interface for listener classes
class LogListener
{
public:
	
	virtual ~LogListener() {}
	virtual void Output(const std::string& msg)=0;

private:


};
