/*   $Revision: 4 $:  $ 
	$Date: 16.05.10 0:09 $  */

/*
(c) 1995-2007 by Joerg Schulla
more version history: see source safe comments
*/

#include <cstring>
#include <cstdlib>
#include "modconst.h"
#include "wasim.h"

using namespace std;


int index(char *string,const char *substr)
/* liefert den Index von string, an welchem zuerst substr auftritt */
{
	if (strstr(string,substr) != NULL)
	return (int)((long)strstr(string,substr)-(long)string); else return(0);
}


char *substr(char *string,int start,int stop)
/* liefert den substring von string zurueck, welcher an der Stelle
start beginnt und an der Stelle stop endet. String wird mit diesem 
Ergebnis ueberschrieben */
{
	char *temp;

	temp = (char *)calloc(_stringsize_,sizeof(char));
	if ((start>=0)&&(stop>=0))
	{
		string[stop+start+1] = '\0';
		strcpy(temp,string+start);
		strcpy(string,temp);
	}
	//  else string[0]='\0';
	free(temp);
	return(string); /* if start or stop are less than zero, the string is not altered */
}

