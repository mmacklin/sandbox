#pragma once

#include "types.h"
#include "log.h"

#include <vector>
#include <string>

// system functions
double GetSeconds();
void Sleep(double seconds);

// helper function to get exe path
std::string GetExePath();
std::string GetWorkingDirectory();

// shows a file open dialog
std::string FileOpenDialog(char *filter = "All Files (*.*)\0*.*\0");

// pulls out an option in the form option=value, must have no spaces
template <typename T>
bool GetCmdLineArg(const char* arg, T& out, int argc, char* argv[])
{
	// iterate over all the arguments
	for (int i=0; i < argc; ++i)
	{
		const char* s1 = arg;
		const char* s2 = argv[i];

		while (*s1 && *s2 && *s1 == *s2)
		{
			++s1;
			++s2;
		}
		
		// we've found the flag we're looking for
		if (*s1 == 0 && *s2 == '=')
		{
			++s2;

			// build a string stream and output			
			std::istringstream is(s2);
			if (is >> out)
			{
				return true;
			}
			else
				return false;
		}
	}

	return false;
}
// return the full path to a file
std::string ExpandPath(const char* path);
// takes a full file path and returns just the folder (with trailing slash)
std::string StripFilename(const char* path);
// strips the path from a file name
std::string StripPath(const char* path);
// strips the extension from a path
std::string StripExtension(const char* path);
// returns the file extension (excluding period)
std::string GetExtension(const char* path);
// normalize path
std::string NormalizePath(const char* path);

// loads a file to a text string
std::string LoadFileToString(const char* filename);
// loads a file to a binary buffer (free using delete[])
uint8_t* LoadFileToBuffer(const char* filename, uint32_t* sizeRead=NULL);
// save whole string to a file
bool SaveStringToFile(const char* filename, const char* s);

bool FileMove(const char* src, const char* dest);
bool FileScan(const char* pattern, std::vector<std::string>& files);

// file system stuff
const uint32_t kMaxPathLength = 2048;

#ifdef WIN32

// defined these explicitly because don't want to include windowsx.h
#define GET_WPARAM(wp, lp)                      (wp)
#define GET_LPARAM(wp, lp)                      (lp)

#define GET_X_LPARAM(lp)                        ((int)(short)LOWORD(lp))
#define GET_Y_LPARAM(lp)                        ((int)(short)HIWORD(lp))

#define vsnprintf _vsnprintf
#define snprintf _snprintf
#define vsnwprintf _vsnwprintf

#if _MSC_VER >= 1400 //vc8.0 use new secure
#define snwprintf _snwprintf_s
#else
#define snwprintf _snwprintf
#endif // _MSC_VER

#endif // WIN32

#if PLATFORM_IOS
inline std::string ExpandPath(const char* p)
{
	NSString *imagePath = [NSString stringWithUTF8String:p];
	NSString *fullPath = [[NSBundle mainBundle] pathForResource:[imagePath lastPathComponent] ofType:nil inDirectory:[imagePath stringByDeletingLastPathComponent]];
	
	if (fullPath)
	{
		std::string s = [fullPath cStringUsingEncoding:1];		
		return s;
	}
	else 
	{
		Log::Info << "Failed to map path for : " << p << std::endl;
		return std::string("");
	}
}
inline std::string GetTempDirectory()
{
	NSString* tmp = NSTemporaryDirectory();
	std::string s = [tmp cStringUsingEncoding:1];
	
	return s;
}

inline std::string DataPath(const char* p)
{
	return ExpandPath((std::string("DataCooked/") + p).c_str());
}

#else

inline std::string ExpandPath(const char* p)
{
	return p;
}

inline std::string DataPath(const char* p)
{
	return ExpandPath((std::string("Data/") + p).c_str());

}

#endif
