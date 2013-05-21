#import <Cocoa/Cocoa.h>

#include <string>
#include "../../core/platform.h"

using namespace std;

string FileOpenDialog(char* filter)
{
	NSOpenPanel *dialogOpen = [NSOpenPanel openPanel];
	[dialogOpen setCanChooseFiles:YES];
	[dialogOpen setCanChooseDirectories:NO];
	[dialogOpen setAllowsMultipleSelection:NO];

	NSMutableArray *typesArray = nil;

	NSString *directory = nil;
	int resultCode = [dialogOpen runModalForDirectory:directory file:nil types:typesArray];	

	if( resultCode == NSOKButton ) {
		NSString *result = [[dialogOpen filenames] objectAtIndex:0];
		return StripPath(string( [result UTF8String] ).c_str());
	}
	else
		return string();
}

