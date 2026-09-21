/*   $Revision: 21 $ 
	$Date: 7.03.11 0:36 $  */

/* Modellierung des Wasserkreislaufes. Voll objektorientiert, von modelobject
abgeleitete Typen sorgen fuer die Interpolation der Eingangsdaten sowie fuer
die Modellierung der Teilprozesse 

(c) 1995-2007 by Joerg Schulla
more version history: see source safe comments

*/

#ifndef _MODAPUZRH_
#define _MODAPUZRH_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <csignal>
#include <ctime>
#include "stringp.h"   /* String-funktionen */
#include "modconst.h"  /* Konstanten fuer das Modell, z.B. Anzahl der Stationen, der Groessen usw. */
#include "modef.h"     /* Definition des Ausgangsobjectes fuer Interpolation und Teilmodelle */
#include "teilmod.h"   /* Definition der Teilmodelle, Code in teilmod.C */
#include "siltingup.h" // definitions for siltingup-module
#include "surfrout.h"  // definitions for surface routing module
#include "dataassim.h"  // definitions for surface routing module

class ModellApplication:
public memory_handling
{
public:

	long lastCouplingHour;        // when was the last coupling hour --> when CouplingMode='H'
	long lastCouplingDay;         // when was the last coupling day --> when CouplingMode='D'
	long lastCouplingMonth;       // when was the last coupling month --> when CouplingMode='M'
	long lastCouplingYear;        // when was the last coupling year --> when CouplingMode='Y'
	long lastCouplingInterval;    // not used 
	long CouplingIntervalCounter; // not used (should be the interval counter, counting n intervals until the next coupling)
	char CouplingMode;            // 'I' for Interval, 'H' for hour, 'D' for day, 'M' for Month and 'Y' for year

	bool	   routingmode;
	int        controll_c_flag;
	int        i;
	char       taste;
	char       **meteonamen;
	int        progaus,warning_level,runoff_measure;
	double      minutes_to_interval_end, MinutesToInvEndRouting;
	int		    timestep;
	int			routingTimestep;
	int			routingTimestepFactor;
	int        numsubb;
	int        *listsubb;
	int        anzahlgroessen;
	double      breite,laenge,bezugsmeridian;
	unsigned short int    programmende;
	int		 iLogLevel;
	int  interceptionMethod;
	bool bInterceptionBeforeSnow;

	char       *steuerdateinamearg, /* der als Argument uebergebene Name der Steuerdatei */
	*steuerdateiname;    /* der Name der nach dem Ersetzen der Symbole in der Steuerdatei
									genutzten Datei mit den ersetzten eintraegen */
	FILE       *steuerdateiarg,
	*steuerdatei;

	char       *satz,*tsatz;

	/* Variablen, die an initial::initial uebergeben werden, um von dieser Routine
	mit daten belegt zu werden. Anschliessend werden diese Daten an alle neu 
	einzurichtenden Instanzen der Teilmodelle uebergeben. Somit haben alle Teilmodelle
	Zugriff auf die wichtigsten globalen Daten */

	long   minute,tag,monat,stunde,jahr,
	mina,tga,sta,mta,jra,
	mine,tge,ste,mte,jre;

	TLastResults  *LastResults;
	int        maxLastResults;

	LUSEtabletype *LUtable;
	SOILtabletype *SOILtable;
	MultiLayerLUSEType *MultiLayerLUTable;
	irrigationtabletype *irrigationtable;
	statrecord    *zone;
	int           *zarray,zonmax,zonmaxwert,ncols,nrows,mpistartrow,mpiendrow;
	int           *rarray,regiomax,regiomaxwert;
	center_of_mass_type *center_of_mass;
	double         xll,yll,cellsize,nodata;
	firstlast     *FirstLastCol;
	gridinfotype  *gridinfo ; /* array mit den Namen und Adressen der Memory-Grids */
	int           gridanzahl;
	int           resultnr;
	double         timekorr; /* Verschiebung der Zeitzone gegenueber der fuer die Notation der Meteodaten genutzten Zeit
							+1 bedeutet, dass die Daten in UTC angegeben sind (MEZ = UTC + 1) */
	int           tracer_count;
	tracer_info_type *tracer_info; /* array with tracer infos for each tracer */
	bool WriteASCIIGrids;
	bool WriteDoubleGrids;
	bool WriteNetcdfGrids;
	char *DefaultInputDirectory;
	char *DefaultOutputDirectory;
	/* Deklarierung der noetigen Object-Variablen (Instanzen) */
	//  void controll_c(int handlerflag); /* control_c handler */

	long *v_mpiregioncnt; // a vector of counts of cells for each MPI region used by Gather and Scatter methods later in the sub models
	std::vector<mpihalo> *sendHalo,*recvHalo; // used for fast eccess of the halo data. mpihalo is of type struct {int row,col,cnt;}

	int doSpinUp;
	int doTemperaturePreSpinUp;
	int repeatCnt;
	std::vector<HRUMapElement*> *vHRU;

	std::vector< std::vector<RiverLinkCell> > *RivLnkCellsToSend;
	std::vector< std::vector<RiverLinkCell> > *RivLnkCellsToReceive;
	std::vector<RiverLink*> *RiverLinkVector;


	defaultini *initial; /* Klasse, die die Standardwerte der Modellierung einliest, Speicherarrays
			definiert (gridinfo) und die Zonenfolge analysiert. Da die Ergebnisse
			dieser Klasse in allen Interpolatoren und in allen Teilmodellen genutzt
			werden, werden die Ergebnisse oder deren Adressen als Initialisierungs-
			werte an die anderen Routinen uebergeben (Zeiger auf die Arrays, Anfangs-
			und Endzeiten usw.) */

	std::vector<SP_class*> reservoirs;
	long *lakeCellsPerMPIRank;

	meteointerpol **meteo; /* ein Array von Zeigern auf
					InterpolationsInstanzen, Code in meteoint.C,
					Definitionen in meteoint.h */
	precicorrclass *precicorr; /* eine Klasse zur Korrektur der Niederschlaege */

	data_assimilation           *data_assimilation_model;
	RegionalSuperposition       *superposition;
	schneemodellklasse          *schneemodell;
	interzeptionsmodellklasse   *interzeptionsmodell;
	siltingupclass		      *siltingupmodell;
	surfaceroutingclass		  *surfaceroutingmodell;
	unsatzonmodellklasse        *unsatzonmodell;
	irrigationclass             *irrigationmodell;
	abflussmodellklasse         *abflussmodell;
	verdunstungsmodellklasse    *verdunstungsmodell;
	expokorrklasse              *expokorrektur;
	gwflowclass                 *gwflow;
	SpecialOutputClass			*specialoutput;

	/* for substance transport */
	tracer_ini_class *tracer_ini;

	ModellApplication(wasim::mystream *Stream,int argc,char *argv[],int steuerdateinummer);
	virtual int run();
	virtual void ClearMemory();
	virtual void m_sleep(clock_t wait);
	virtual bool WaitForFile(char *filename,long waitinterval);
	void WriteSemaphore(const char *filename,char *content);
	std::vector<std::pair<int,int>>* getSubSteps(int timestep, int routingTimestep, int minutes_to_interval_end, int MinutesToInvEndRouting, int year, int month, int day, int hour, int minute);
	void RunPreSpinUp();

};

#endif
