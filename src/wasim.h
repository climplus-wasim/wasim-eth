#ifndef _WASIMH_
#define _WASIMH_



//artificial change to set the new version label in svn

#if defined __TIME__ 
#define _wasimuzr_version_english_ "Richards 11.02.04, last change: November, 04, 2026 12:15 compiled at " __DATE__ ", " __TIME__ 
#define _wasim_version_english_    "Topmodel 11.02.04, last change: November, 04, 2026 12:15 compiled at " __DATE__ ", " __TIME__ 
#else														 		   
#define _wasim_version_english_    "Topmodel 11.02.04, last change: November, 04, 2026 12:15"
#define _wasimuzr_version_english_ "Richards 11.02.04, last change: November, 04, 2026 12:15"
#endif


void SetSignalAbfangfunktion();
void SetSignalDefaults();

#endif

