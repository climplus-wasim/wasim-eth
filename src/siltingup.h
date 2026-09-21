/*   $Revision: 13 $:  $ 
	$Date: 7.03.11 0:36 $  */

/*

(c) 2007 by Joerg Schulla
header file, containing declarations for surface sealing or silting-up-module for WaSiM

*/

#ifndef _SILTINGUPH_
#define _SILTINGUPH_

#include <iostream>
#include <string>
#include <vector>
#include "parser.h"
#include "teilmod.h"

//using namespace std;

class siltingupclass:
public teilmodellklasse
{
public:

	int method;
	bool SurfaceRouting;
	bool SurfaceRoutingChecked;

	// Grids which are already in the gridinfo-array
	double   **DHM,  // digital elevation model --> not yet really used
	**P,    // precipitation (either output from interception, snow-melt or precipitation+irrigation)
	***LUSE,// land use
	***VCF, // vegetation covered fraction (from landuse) for each layer
	**SOIL, // soil types
	**SSSU; // surface storage for silting up module = water which was routed in surface routing module but never reinfiltrated or flew into a river

	// names of internal and or output grids
	double   **SDISP, // silting disposition = Verschlaemmungsneigung C -> should be output grid since it's time dependent for later re-starts with reading in old grids
	**RAINPAUSE,
	**EKIN,
	**TSBB,
	**QDSU; // QDSU=QD Silting Up -> this fraction must be handled in unsatzonmodel as additional input for the QDgrid like SnowRunoff
	// more grids may be required, so define them here....
	//double **GRID1,
	//      **GRID2,
	// a.s.o. ...	
	//      **GRIDN;
	double *RainPause,*ParB; // an example for a parameter which must be read in for each subbasin or zone
	int    *code;	// zone codes

	// writecodes for outputgrids
	int      SDISPwritegrid,
	RAINPAUSEwritegrid,
	EKINwritegrid,
	TSBBwritegrid,
	QDSUwritegrid;

	char *SDISPgridname,
	*RAINPAUSEgridname,
	*EKINgridname,
	*TSBBgridname,
	*QDSUgridname;

	int sdispindex,
	rainpauseindex,
	ekinindex,
	tsbbindex,
	qdsuindex;

	// variables required for statistic tables
	zonalstatistic  *SDISPstat,
	*QDSUstat;

	statrecord	  *zoneSDISP,
	*zoneQDSU;      

	int  SDIPSstatcode,
	QDSUstatcode;

	char *SDISPstatfilename,    
	*QDSUstatfilename;

	int iNumVegetationLayers;
	bool interceptionBeforeSnow;

	// define additional surface sealing or sliting-up-model parameters here
	// int    ...
	// double  ...
	// double ...
	double parA; // for method 2:  parameter A in I_0 = A (initial infiltration capacity, in method 0 defined as 65 mm/h)
	double parB; // for method 2:  parameter B in I_end = B*(dg^C)*(fd^D) (in method 0 defined as 12.2)
	double parC; // for method 2:  parameter C in I_end = B*(dg^C)*(fd^D) (in method 0 defined as 0.52)
	double parD; // for method 2:  parameter D in I_end = B*(dg^C)*(fd^D) (in method 0 defined as -0.64)
	double parE; // for method 2:  parameter E in Cv = E*(fd^F)*(dg^G)*(t_cult^H) (in method 0 defined as 0.013)
	double parF; // for method 2:  parameter F in Cv = E*(fd^F)*(dg^G)*(t_cult^H) (in method 0 defined as -1.03)
	double parG; // for method 2:  parameter G in Cv = E*(fd^F)*(dg^G)*(t_cult^H) (in method 0 defined as 0.7)
	double parH; // for method 2:  parameter H in Cv = E*(fd^F)*(dg^G)*(t_cult^H) (in method 0 defined as -0.19)

	std::vector<std::string> expressionlist;

	std::map<std::string,double> mapVars;

	// constructor of the class
	siltingupclass(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,int maxLastResults,
	double breitearg,double laengearg,double bezugsmeridianarg,char *bezeichnung,
	char *steuername,int LogLevel,int progaus,
	double p_min_to_int_end,int p_numsubb,int *p_listsubb,
	gridinfotype *p_gridinfo,
	int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
	int p_zonmax,  int p_zonmaxwert,
	int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
	double p_xll ,  double p_yll,
	double p_cellsize    , double      p_nodata,
	long mine,long ste , long tge , long mte , long jre, double p_timekorr,
	LUSEtabletype *p_LUtable,
	SOILtabletype *p_SOILtable,
	MultiLayerLUSEType *p_MultiLayerLUTable,
	bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
	long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
	int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU);

	virtual void   run(long me,long se,long tg,long mt,long jr);
	virtual void   ClearMemory();
	virtual void   flushbuffers();

	/* overwriting the inherited functions for substance transport from teilmodellklasse */
	virtual void substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults);
	virtual void substance_transport_delete(); /* free memeory */
	virtual void setInterceptionBeforeSnow(bool b_interceptionBeforeSnow);
};

#endif

