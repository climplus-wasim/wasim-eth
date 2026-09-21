/*   $Revision: 20 $:  $ 
	$Date: 8.06.10 0:05 $  */

/*

(c) 2008 by Joerg Schulla
header file, containing declarations for surface routing module for WaSiM

*/

#ifndef _SURFROUTH_
#define _SURFROUTH_

#include <iostream>
#include "teilmod.h"

//using namespace std;

typedef struct {
	double DT;
	double dhm[5];
	double roughx[5];
	double roughy[5];
	double rivdep[5];
	double rivwit[5];
	double head[5];
	double newhead;
	double h[5];
	double T[5];
	double f[5];
	double dx;
	double dy;
	double headold;
	double Tsum;
	double fsum;
	double atemp;
	double btemp;
	double ctemp;
} DumptypeSR;

class surfaceroutingclass:
public teilmodellklasse
{
public:

	int method;
	bool SiltingUp; // floag if silting up model is used. If not, then CONC_SSSU must be created in SurfRout itself and the decay has to be calculated.

	// Grids which are already in the gridinfo-array
	double   **DHM,  // digital elevation model, may be usefull for e.g. calculating special flow directions
	**T,    // temperature, will be usefull for e.g. Reynolds Number
	***LUSE,// land use, usefull for surface parameters for routing, like roughness
	***VCF,  // Vegetation covered Fraction
	**SOIL, // soil types
	**SLOPE, // slope, important for flow velocities on the surface
	**FLD, // flow direction
	**SSSU, // existing grid from variable grids. called also SSSU in silting up module and having the identifier _SURFSTORSILTINGUP_ 
	**QDgrid, // direct flow or surface runoff
	**QIgrid, // interflow grid
	***CONC_QD, // tracer concentrations for surface runoff
	***CONC_QI, // tracer concentration for interflow
	***CONC_SSSU, // tracer concentration from silting up model
	**TAUCRIT, // critical shear stress (in Pa) for starting channel creation
	**RIVDEP, // river depth
	**RIVWIT, // river witdh
	**RIVLNK, // river link grid for routing (used for making a cell a subcatchment outlet if no slope to a cell inside the catchment is there
	**FZT;    // flow time grid (accumulated mean flow travel times to the next subbasin outlet for surface and channel runoff in secounds)

	// the resulting water flows and tracer fluxes which will be given to the routing model as total runoff
	double	  **QB,  // not a grid but an array of arrays, containing zonal baseflow for water and tracer fluxes
	**QI,  // the same for interflow, which may be created within this module or taken from unsatzonmodel
	**QIi, // interflow with Single linear storage (ELS) effect applied
	**QD,  // and the same again for surface runoff
	**QDd, // direct runoff with Single linear storage (ELS) effect aplied
	**Qg;  // total runoff as sum of surface runoff, interflow and base flow

	// names of internal and or output grids which will be created within this module
	// more grids may be required, so define them here....
	double ***SRFCVEL,**SRFCVEL1; // a stack of grids, for each direction one grid with 1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW,7=W,8=NW
	double ***FLUX,**FLUX1; // a stack of grids, for each direction one grid with 1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW,7=W,8=NW
	double ***WEIGHT,**WEIGHT1; // a stack of grids, for each direction one grid with 1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW,7=W,8=NW
	double **HEAD,**HEADold,**PLANCURV; // total input into a cell for a subintervall (will be used for calculating slopes)
	double **SRQIN,**SRQOUT; // inflows and outflows within a sub time step
	double **ROUGHX,**ROUGHY,**MANNING; // used for new Gauss-Seidl approach 
	double **CellSizeX, **CellSizeY;
	//      **GRID2,
	// a.s.o. ...	
	//      **GRIDN;

	double global_mintime;
	int old_year,old_month, old_day, old_hour, old_minute;
	int moduleSubstepCounter;

	// writecodes for outputgrids
	// SRFCVEL = Surface Velocity
	int SRFCVELwritegrid,
	FLUXwritegrid;

	char *SRFCVELgridname,
	*FLUXgridname;

	int  srfcvelindex,
	fluxindex,
	headindex, 
	headoldindex, 
	roughxindex, 
	roughyindex, 
	manningindex,
	plancurvindex,
	weightindex,
	srqinindex,
	srqoutindex,
	taucritindex;

	int  concsssuindex[10];

	double *kelsqi; // parameter is taken from unsatzonmodell
	double *kelsqd; // parameter is taken from unsatzonmodell


	// variables required for statistic tables

	statrecord *zoneHEADStart, *zoneHEADEnd, *zoneQDFromFlux;
	double *QGScalings, *numCellsQD;
	double tempsum;

	statrecord ***zoneFZT_QD,
	***zoneFZT_QI;  /* fuer die Translation des Oberflaechenabflusses wird
							pro Teilgebiet pro Zeitflaechenzone ein Statistikregistersatz
							benoetigt, ausserdem ein Array mit der Information ueber die
							maximale Anzahl an Registern in jedem Teilgebiet (maxzeit) */

	zonalstatistic  *QDSRstat,
	*QISRstat,
	*QBSRstat,
	*QGSRstat,
	*SFSTSRstat;

	statrecord	  *zoneQDSR,
	*zoneQISR,
	*zoneQBSR,
	*zoneQGSR,
	*zoneSFSTSR;

	int QDSRstatcode,
	QISRstatcode,
	QBSRstatcode,
	QGSRstatcode,
	SFSTSRstatcode;

	char *QDSRstatfilename,
	*QISRstatfilename,
	*QBSRstatfilename,
	*QGSRstatfilename,
	*SFSTSRstatfilename;

	int iNumVegetationLayers;

	// define additional surface routing parameters here
	int    iMaxIterANL,
	iMaxIterVelocity;

	// double  ...

	double dMaxDeltaANL,
	dMaxDeltaVelocity,
	dMinSubTimeStep,
	dMaxSubTimeStep;
	double dConcFact,dMinFlowDepth;

	double       *maxzeit; // for flow time zones in FZT and ***zoneFZT_QD and ***zoneFZT_QI

	/* variables for the storage.ftz file */
	FILE *storagefile; /* a text file containing the water content of all flow time zones and
					   the actual content of the single linear reservoirs for all subbasins
					   for both, surface runoff and interflow */
	char *storagefilename;

	bool dumpmode;
	bool closeAfterDump;
	bool repairAfterDump;
	long maxDumpFiles;
	long actDumpFiles;

	// constructor of the class
	surfaceroutingclass(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,int maxLastResults,
						double breitearg,double laengearg,double bezugsmeridianarg,char *bezeichnung,
						char *steuername,int LogLevel,int progaus,
						double p_min_to_int_end,int p_numsubb,int *p_listsubb,
						gridinfotype *p_gridinfo,
						int *p_zarray,  statrecord *p_zone, firstlast *p_FirstLastCol ,
						int p_zonmax,  int p_zonmaxwert,
						int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
						double p_xll ,  double p_yll,
						double p_cellsize    , double      p_nodata,
						long mine, long ste , long tge , long mte , long jre, double p_timekorr,
						LUSEtabletype *p_LUtable,
						SOILtabletype *p_SOILtable,
						MultiLayerLUSEType *p_MultiLayerLUTable,
						bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
						long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
						int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU,int substepcnt=1);

	virtual void setSubStepCnt(int step); // inherited from teilmodellklasse and changed here

	virtual void   fzt_grid_auswerten();
	
	virtual void   pre_run(long me,long se,long tg,long mt,long jr, long actsubstep);
	virtual void   run(long me, long se,long tg,long mt,long jr, long actsubstep);
	virtual void   post_run(long me, long se,long tg,long mt,long jr, long actsubstep);

	virtual void   ClearMemory(bool appendmode);
	virtual void   flushbuffers();

	virtual double  **GetQgesAdresse();
	virtual void	 SetQBAddress(double **address);
	virtual void	 SetQIAddress(double **address);
	virtual void	 SetQDAddress(double **address);
	virtual void	 SetKelsQDAddress(double *address);
	virtual void	 SetKelsQIAddress(double *address);

	virtual void   getnachbar(int i,int j,double **Grid);
	virtual void   GetPlanCurvature();
	virtual void   GetSlopeAspectAndCelltype(int i,int j,double &slope,double &aspect,int &celltype,int step);
	virtual void   GetMultipleFlowPath(int i,int j,double slope,double aspect,int &celltype,int step);
	virtual void   GetFlowVelocities(int i,int j,int step,int NextMajorStep,int celltype);
	virtual double  MinimumFlowTime(int step);
	virtual void   CalculateWaterFlows(int i,int j,double mintime);
	//  virtual double  GetRoughness(int i,int j,double v);
	virtual double GetRoughness(int i,int j,double v,int luse,int SR_cultivation,double RH,double VCFlocal,double Ie);
	virtual void effective_roughness();
	virtual double gauss_seidel_iteration(double DT,int iter);
	virtual void GetFlowVelocities_AfterGaussSeidel();

	/* overwriting the inherited functions for substance transport from teilmodellklasse */
	virtual void   substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults);
	virtual void   substance_transport_delete(); /* free memeory */
	virtual void  ReadSlrContents(char *_storagefilename, int flag);
	virtual void  WriteSlrContents(bool appendmode);
	virtual void  ReadFloatFromStoragefile(double *FValue, const char *ErrorMessage);
	virtual void  StorageFileError(const char *message, double c1, double c2);

	virtual void  setDumpMode(bool p_dumpmode, bool p_closeAfterDump, bool p_repairAfterDump, long p_maxDumpFiles) { this->dumpmode = p_dumpmode; this->closeAfterDump = p_closeAfterDump; this->repairAfterDump = p_repairAfterDump; this->maxDumpFiles = p_maxDumpFiles; };

};

#endif