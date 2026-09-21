#ifndef _WASIMH_
#define _WASIMH_



//artificial change to set the new version label in svn

#if defined __TIME__ 
#define _wasimuzr_version_english_ "Richards 11.01.01, last change: July, 26, 2026 09:32 compiled at " __DATE__ ", " __TIME__ 
#define _wasim_version_english_    "Topmodel 11.01.01, last change: July, 26, 2026 09:32 compiled at " __DATE__ ", " __TIME__ 
#else														 		   
#define _wasim_version_english_    "Topmodel 11.01.01, last change: July, 26, 2026 09:32"
#define _wasimuzr_version_english_ "Richards 11.01.01, last change: July, 26, 2026 09:32"
#endif


void SetSignalAbfangfunktion();
void SetSignalDefaults();

#endif

