/*   $Revision: 24 $ 
	$Date: 16.05.10 0:09 $  */

/*
(c) 1995-2007 by Joerg Schulla
more version history: see source safe comments

Quelltexte fuer Interzeptionsmodell.
Alle Teilmodelle sind Ableitungen der
virtuellen Klasse 'modelobject', wobei vor der Bildung eigener Klassen fuer die
Teilmodelle noch die Definition einer gemeinsamen Teilmodellklasse steht,
welche allen Teilmodellen gemeinsame Routinen enthaelt.

Version 1: Juni 1995

Version 2: 19.9.97: 
				Die Verdunstung ETR wird nicht mehr in verdunstungsmodell::potentiell bereits bestimmt und hier durch die
				Interzeptionsverdunstung nur noch modifiziert, sondern nun wird die potentielle Verdunstung
				hier durch Interzeptionsverluste modifiziert

Version 3: 16./19. Maerz 1998: Stofftransport implementiert

date       : march, 9 1999
version    : 4, hour-entry can be given for start or for end of time step via min_to_int_end

date       : december, 11 1999
version    : stream output into all possible childs of baseclass 
				wasim::mystream now possible. A stream has to be given 
				to the Application object at initialization (constructor 
				parameter).
				Also new: Exception handling is used for exit from the program.
				All old calls of "exit(n);" are replaced by a "throw n;" in order
				to allow the modules to be used in a DLL under the WIN32-environment.

date       : January, 2006
version    : implementing multi layer interception in conjunction with the new multilayer-transpiration model

*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "stringp.h"   /* String-funktionen */
#include "modconst.h"
#include "modef.h"     /* Definition des Ausgangsobjectes der Interpolationen und Teilmodelle, Code in modef.C */
#include "teilmod.h"   /* Definition der Teilmodelle */
#include "wasim.h"

using namespace std;
using namespace snow_canopy;

/* ------------- Definition einer Klasse fuer das Interzeptionsmodell ----------------- */

GridStack::GridStack(int capacity, int mpistartrow, int mpiendrow, int rows, int cols) :modelobject()
{
	this->stack = new std::vector<double**>;
	this->mpistartrow = mpistartrow;
	this->mpiendrow = mpiendrow;
	this->nrows = rows;
	this->ncols = cols;
	this->n = 0; // counter of actually hold grids in the stack
	this->capacity = capacity;
	if (capacity < 1) capacity = 1;
}

bool GridStack::addToStack(double** grid_in)
{
	this->stack->push_back(grid_in);
	// keep the last 24 hours, older data will be deleted
	if (this->stack->size()>capacity)
	{
		free_matrix(this->stack->front(), mpistartrow, mpiendrow, 1, ncols, nrows);
		//delete this->stack->front();
		this->stack->erase(this->stack->begin());
	}
	this->n = this->stack->size();
	return true;
}

double GridStack::getAvg(int x, int y)
{
	double average = 0.;
	int n = this->stack->size();
	for (int i = 0; i<n; i++)
	{
		double** grid_i = this->stack->at(i);
		average += grid_i[x][y];
	}
	average /= n;
	return average;
}

GridStack::~GridStack()
{
	while (this->stack->size() > 0)
	{
		free_matrix(this->stack->back(), mpistartrow, mpiendrow, 1, ncols, nrows);
		//delete this->stack->back();
		this->stack->pop_back();
		this->n = stack->size();
	}
	delete this->stack;
}

interzeptionsmodellklasse::interzeptionsmodellklasse(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,
	int maxLastResults,double breitearg,double laengearg,double bezugsmeridianarg,
	char *bezeichnung,char *steuername,int LogLevel,
	int progaus,double p_min_to_int_end,int p_numsubb,int *p_listsubb,
	gridinfotype *p_gridinfo,int *p_zarray,statrecord *p_zone,firstlast *p_FirstLastCol ,
	int p_zonmax,int p_zonmaxwert,
	int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow, double p_xll, double p_yll, double p_cellsize, double p_nodata,
	long mine,long ste , long tge , long mte , long jre, double p_timekorr,
	LUSEtabletype *p_LUtable,SOILtabletype *p_SOILtable,MultiLayerLUSEType *p_MultiLayerLUTable,
	bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
	long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
	int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU):
teilmodellklasse(Stream,
	breitearg,laengearg,bezugsmeridianarg,bezeichnung,steuername,LogLevel,progaus,p_min_to_int_end,p_numsubb,p_listsubb,
	p_gridinfo,p_zarray,p_zone,p_FirstLastCol,p_zonmax,p_zonmaxwert,
	p_ncols, p_nrows, p_mpistartrow, p_mpiendrow, p_xll, p_yll, p_cellsize, p_nodata, mine,ste, tge, mte, jre, p_timekorr,
	p_LUtable, p_SOILtable, p_MultiLayerLUTable, WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory, v_mpiregioncnt, p_sendHalo, p_recvHalo,
	p_doSpinUp, p_routingmode, p_vHRU),
EIstat(0),SIstat(0),Qstat(0),EIstatfilename(0),SIstatfilename(0),Qstatfilename(0),
EIstatcode(0),SIstatcode(0),Qstatcode(0),EIgridname(0),SIgridname(0),Qgridname(0),
QInput(0),ETP(0),EIP(0),ETR(0),LUSE(0),EI(0),SI(0),Q(0),dhm(0),WURZ(0),LAI(0),VBG(0),ICE(0),
siindex(0),eiindex(0),qindex(0),EIwritegrid(0),SIwritegrid(0),Qwritegrid(0),
SImax(0),SIdicke(0),tagnr(0),lai(0),vbg(0),use(0),Z0(0),SLIQ(0),RSC(0),RSI(0),Albedo(0),SSTO(0),
zoneEI(NULL), zoneSI(NULL), zoneQ(NULL),
TInput(NULL), RInput(NULL), WSInput(NULL), RHInput(NULL), CLInput(NULL),
SCN_Albedo(NULL), SCN_LOAD(NULL), SCN_AGE(NULL), SCN_SUBL(NULL), SCN_Out_solid(NULL), SCN_Out_liquid(NULL),
SCN_Temp(NULL), SCN_RH(NULL), SCN_RAD(NULL), SCN_Wind(NULL), SCN_IR(NULL), CEB_CANENBAL(NULL), CEB_DUMMY(NULL), CEB_CANTEMP(NULL), CEB_CANTEMP_OLD(NULL),
grd_idx_scn_temp(0), grd_idx_scn_rh(0), grd_idx_scn_rad(0), grd_idx_scn_wind(0), grd_idx_scn_ir(0), grd_idx_scn_alb(0), grd_idx_scn_load(0),
grd_idx_scn_age(0), grd_idx_scn_subl(0), grd_idx_scn_out_s(0), grd_idx_scn_out_l(0), grd_idx_ceb_canenbal(0), grd_idx_ceb_dummy(0), grd_idx_ceb_cantemp(0), grd_idx_ceb_cantemp_old(0),
SCN_Temp_stats(NULL), SCN_RH_stats(NULL), SCN_RAD_stats(NULL), SCN_Wind_stats(NULL), SCN_Alb_stats(NULL), SCN_Age_stats(NULL), SCN_Load_stats(NULL), SCN_Sub_stats(NULL), SCN_Out_solid_stats(NULL), SCN_Out_liq_stats(NULL),
SCN_Temp_zone(NULL), SCN_RH_zone(NULL), SCN_RAD_zone(NULL), SCN_Wind_zone(NULL), SCN_Alb_zone(NULL), SCN_Age_zone(NULL), SCN_Load_zone(NULL), SCN_Sub_zone(NULL), SCN_Out_solid_zone(NULL), SCN_Out_liq_zone(NULL),
CEB_CANENBAL_stats(NULL), CEB_DUMMY_stats(NULL), CEB_CANTEMP_stats(NULL), CEB_CANTEMP_OLD_stats(NULL), CEB_CANENBAL_zone(NULL), CEB_DUMMY_zone(NULL), CEB_CANTEMP_zone(NULL), CEB_CANTEMP_OLD_zone(NULL),
bSCN_model(false), CEB_model(false),
SCN_Temp_filename(NULL), SCN_RH_filename(NULL), SCN_RAD_filename(NULL), SCN_Wind_filename(NULL),
SCN_ALB_filename(NULL), SCN_LOAD_filename(NULL), SCN_AGE_filename(NULL), SCN_SUBL_filename(NULL), SCN_Out_solid_filename(NULL),
SCN_Out_liq_filename(NULL), CEB_CANENBAL_filename(NULL), CEB_DUMMY_filename(NULL), CEB_CANTEMP_filename(NULL), CEB_CANTEMP_OLD_filename(NULL),
SCN_Temp_gridname(NULL), SCN_RH_gridname(NULL), SCN_RAD_gridname(NULL), SCN_Wind_gridname(NULL),
SCN_ALB_gridname(NULL), SCN_LOAD_gridname(NULL), SCN_AGE_gridname(NULL), SCN_SUBL_gridname(NULL), SCN_Out_solid_gridname(NULL), SCN_Out_liq_gridname(NULL),
CEB_CANENBAL_gridname(NULL), CEB_DUMMY_gridname(NULL), CEB_CANTEMP_gridname(NULL), CEB_CANTEMP_OLD_gridname(NULL),
SCN_Temp_code(0), SCN_RH_code(0), SCN_RAD_code(0), SCN_Wind_code(0),
SCN_ALB_code(0), SCN_LOAD_code(0), SCN_AGE_code(0), SCN_SUBL_code(0), SCN_Out_solid_code(0), SCN_Out_liq_code(0),
CEB_CANENBAL_code(0), CEB_DUMMY_code(0), CEB_CANTEMP_code(0), CEB_CANTEMP_OLD_code(0),
SCN_Temp_writegrid(0), SCN_RH_writegrid(0), SCN_RAD_writegrid(0), SCN_Wind_writegrid(0),
SCN_ALB_writegrid(0), SCN_LOAD_writegrid(0), SCN_AGE_writegrid(0), SCN_SUBL_writegrid(0), SCN_Out_solid_writegrid(0),
SCN_Out_liq_writegrid(0), CEB_CANENBAL_writegrid(0), CEB_DUMMY_writegrid(0), CEB_CANTEMP_writegrid(0), CEB_CANTEMP_OLD_writegrid(0),
SCN_factor_lai(0), SCN_min_lai(0), SCN_min_z0(0), transzone(0), T0R(0),
SURFTEMP(NULL), TempGrids(NULL), stunde_old(-1), tag_old(-1), monat_old(-1), jahr_old(-1),
do_ArcVeg(0), arcVegNumLayers(0), ARCVEGFRACTIONS(NULL)
{
	try
	{
		retvalue=0;
		int i,j,k;


		if (routingmode)
		{
			do_it = 0; // wether the submodels code for do_it is set to 0 or 1: when routingmode is active, this submodel should never run
		}

		if (do_it == 0) return;
#ifdef __TIME__
		(*WSstream) << "file " __FILE__ ", compiled at " __DATE__ ", " __TIME__ <<"[end-entry]";
#endif
		/* reservieren von Speicher fuer die Dateinamen */

		if (MPIREG == NULL)
		{
			if (ist_vorhanden(_MPIRegions_))
			{
				MPIREG = gridinfo[getindex(_MPIRegions_)].pointer;
			}
			else
			{
				MPIREG = NULL;
			}

			if (ist_vorhanden(_MPIHelper_))
			{
				MPIHLP = gridinfo[getindex(_MPIHelper_)].pointer;
			}
			else
			{
				MPIHLP = NULL;
			}
		}

		if (MultiLayerLUTable)
		{
			iNumVegetationLayers=MultiLayerLUTable[1].getNumLayers();
		}
		else
		{
			iNumVegetationLayers=1;
		}

		EIstat = new zonalstatistic*[iNumVegetationLayers+1];
		SIstat = new zonalstatistic*[iNumVegetationLayers+1];
		Qstat  = new zonalstatistic*[iNumVegetationLayers+1];

		EIstatfilename = new char*[iNumVegetationLayers+1];
		SIstatfilename = new char*[iNumVegetationLayers+1];
		Qstatfilename  = new char*[iNumVegetationLayers+1];

		EIgridname = new char*[iNumVegetationLayers+1];
		SIgridname = new char*[iNumVegetationLayers+1];
		Qgridname  = new char*[iNumVegetationLayers+1];

		LUSE  = new double**[iNumVegetationLayers+1];
		WURZ  = new double**[iNumVegetationLayers+1];
		LAI   = new double**[iNumVegetationLayers+1];
		VBG   = new double**[iNumVegetationLayers+1];
		ETP   = new double**[iNumVegetationLayers+1];
		EIP   = new double**[iNumVegetationLayers+1];
		ETR   = new double**[iNumVegetationLayers+1];
		EI    = new double**[iNumVegetationLayers+1];
		SI    = new double**[iNumVegetationLayers+1];
		Q     = new double**[iNumVegetationLayers+1];
		Z0    = new double**[iNumVegetationLayers+1];
		RSC   = new double**[iNumVegetationLayers+1];
		RSI   = new double**[iNumVegetationLayers+1];

		zoneEI	= new statrecord*[iNumVegetationLayers+1];
		zoneSI	= new statrecord*[iNumVegetationLayers+1];
		zoneQ	= new statrecord*[iNumVegetationLayers+1];

		siindex = new int[iNumVegetationLayers+1];
		eiindex = new int[iNumVegetationLayers+1];
		qindex  = new int[iNumVegetationLayers+1];

		int lay;
		for (lay=0;lay<=iNumVegetationLayers;lay++)
		{
			Qgridname[lay]      = newstring(_stringsize_+1,1);
			EIgridname[lay]     = newstring(_stringsize_+1,1);
			SIgridname[lay]     = newstring(_stringsize_+1,1);

			Qstatfilename[lay]  = newstring(_stringsize_+1,1);
			EIstatfilename[lay] = newstring(_stringsize_+1,1);
			SIstatfilename[lay] = newstring(_stringsize_+1,1);
		}
		
                
		// read parameters and filenames 
		fscanfln(steuerdatei,"%d",&methode," method for ",tstr);
		// methods for interception model, supported by version 10.04.00 and higher
		// 1 = traditional / legacy: using ETP(no separate EIP), first snow model, then interception
		// 2 = traditional / legacy: using separate EIP, first snow model, then interception
		// 11 = using ETP(no separate EIP), first interception model, then snow model
		// 12 = using separate EIP, first interception model, then snow model
		// 21 = using ETP(no separate EIP), first interception model, then snow model, +Snow canopy interception SCN
		// 22 = using separate EIP, first interception model, then snow model, +Snow canopy interception SCN
		// 31 = using ETP(no separate EIP), first interception model, then snow model, +Snow canopy interception SCN + Canopy energy balance CEB
		// 32 = using separate EIP, first interception model, then snow model, +Snow canopy interception SCN + Canopy energy balance CEB

		if ((methode != 1) && (methode != 2) && (methode != 11) && (methode != 12) && (methode != 21) && (methode != 22) && (methode != 31) && (methode != 32))
		{
			(*WSstream) << "FATAL: method "<<methode<<" is not supported in interception model (use one of 1, 2, 11, 12, 21, 22, 31, 32)" << "[end-entry]";
			throw(1161);
		}
		fscanfln(steuerdatei,"%s",&(Qgridname[0][0])         ,"name outflow grid for "    ,tstr);
		fscanfln(steuerdatei,"%d",&Qwritegrid      ,"write code outflow grid for "        ,tstr);
		fscanfln(steuerdatei,"%s%d",&(Qstatfilename[0][0]),&Qstatcode  ,"name outflow statistics file for "      ,tstr);
		fscanfln(steuerdatei,"%s",&(EIgridname[0][0])        ,"name potential interception evaporation grid for "         ,tstr);
		fscanfln(steuerdatei,"%d",&EIwritegrid    ,"write code potential interception evaporation grid for "              ,tstr);
		fscanfln(steuerdatei,"%s%d",&(EIstatfilename[0][0]),&EIstatcode,"name potential interception statistics file for "     ,tstr);
		fscanfln(steuerdatei,"%s",&(SIgridname[0][0])        ,"name interception storage grid for "         ,tstr);
		fscanfln(steuerdatei,"%d",&SIwritegrid    ,"write code interception storage grid for "              ,tstr);
		fscanfln(steuerdatei,"%s%d",&(SIstatfilename[0][0]),&SIstatcode,"name interception storage statistics file for "     ,tstr);

		fscanfln(steuerdatei,"%lf",&SIdicke        ,"parameter interception storage maximum water layer thickness for "                   ,tstr);
		fscanfln(steuerdatei,"%d",&gridslesen   ,"readgrid parameter for ",tstr);
		int run_scn = 0;
		int run_ceb = 0;

		// $kf 2016-07-18: initialize variables for input from control file
		if (methode >= 21)
		{
			// switch if model should be used
			run_scn = 1;

			this->SCN_LOAD_gridname = newstring(_stringsize_ + 1, 1);
			this->SCN_AGE_gridname = newstring(_stringsize_ + 1, 1);
			this->SCN_ALB_gridname = newstring(_stringsize_ + 1, 1);
			this->SCN_SUBL_gridname = newstring(_stringsize_ + 1, 1);
			this->SCN_Out_solid_gridname = newstring(_stringsize_ + 1, 1);
			this->SCN_Out_liq_gridname = newstring(_stringsize_ + 1, 1);
			this->SCN_Temp_gridname = newstring(_stringsize_ + 1, 1);
			this->SCN_RH_gridname = newstring(_stringsize_ + 1, 1);
			this->SCN_RAD_gridname = newstring(_stringsize_ + 1, 1);
			this->SCN_Wind_gridname = newstring(_stringsize_ + 1, 1);

			this->SCN_LOAD_filename = newstring(_stringsize_ + 1, 1);
			this->SCN_AGE_filename = newstring(_stringsize_ + 1, 1);
			this->SCN_ALB_filename = newstring(_stringsize_ + 1, 1);
			this->SCN_SUBL_filename = newstring(_stringsize_ + 1, 1);
			this->SCN_Out_solid_filename = newstring(_stringsize_ + 1, 1);
			this->SCN_Out_liq_filename = newstring(_stringsize_ + 1, 1);
			this->SCN_Temp_filename = newstring(_stringsize_ + 1, 1);
			this->SCN_RH_filename = newstring(_stringsize_ + 1, 1);
			this->SCN_RAD_filename = newstring(_stringsize_ + 1, 1);
			this->SCN_Wind_filename = newstring(_stringsize_ + 1, 1);

			// grids for modified sub-canopy meteorological fields
			fscanfln(steuerdatei, "%s", &(this->SCN_Temp_gridname[0]), "name sub-canopy temperature grid for ", tstr);
			fscanfln(steuerdatei, "%d", &this->SCN_Temp_writegrid, "write code sub-canopy temperature grid for ", tstr);
			fscanfln(steuerdatei, "%s%d", &(this->SCN_Temp_filename[0]), &this->SCN_Temp_code, "name sub-canopy stemperature statistics file for ", tstr);
			fscanfln(steuerdatei, "%s", &(this->SCN_RH_gridname[0]), "name sub-canopy relative humidity grid for ", tstr);
			fscanfln(steuerdatei, "%d", &this->SCN_RH_writegrid, "write code sub-canopy relative humidity grid for ", tstr);
			fscanfln(steuerdatei, "%s%d", &(this->SCN_RH_filename[0]), &this->SCN_RH_code, "name sub-canopy relative humidity statistics file for ", tstr);
			fscanfln(steuerdatei, "%s", &(this->SCN_RAD_gridname[0]), "name sub-canopy shortwave radiation grid for ", tstr);
			fscanfln(steuerdatei, "%d", &this->SCN_RAD_writegrid, "write code sub-canopy shortwave radiation grid for ", tstr);
			fscanfln(steuerdatei, "%s%d", &(this->SCN_RAD_filename[0]), &this->SCN_RAD_code, "name sub-canopy shortwave radiation statistics file for ", tstr);
			fscanfln(steuerdatei, "%s", &(this->SCN_Wind_gridname[0]), "name sub-canopy windspeed grid for ", tstr);
			fscanfln(steuerdatei, "%d", &this->SCN_Wind_writegrid, "write code sub-canopy windspeed grid for ", tstr);
			fscanfln(steuerdatei, "%s%d", &(this->SCN_Wind_filename[0]), &this->SCN_Wind_code, "name sub-canopy windspeed statistics file for ", tstr);
			// grids for state variables
			fscanfln(steuerdatei, "%s", &(this->SCN_LOAD_gridname[0]), "name snow canopy load grid for ", tstr);
			fscanfln(steuerdatei, "%d", &this->SCN_LOAD_writegrid, "write code snow canopy load grid for ", tstr);
			fscanfln(steuerdatei, "%s%d", &(this->SCN_LOAD_filename[0]), &this->SCN_LOAD_code, "name snow canopy load statistics file for ", tstr);
			fscanfln(steuerdatei, "%s", &(this->SCN_AGE_gridname[0]), "name snow canopy load grid for ", tstr);
			fscanfln(steuerdatei, "%d", &this->SCN_AGE_writegrid, "write code snow canopy load grid for ", tstr);
			fscanfln(steuerdatei, "%s%d", &(this->SCN_AGE_filename[0]), &this->SCN_AGE_code, "name snow canopy load statistics file for ", tstr);
			fscanfln(steuerdatei, "%s", &(this->SCN_ALB_gridname[0]), "name snow canopy load grid for ", tstr);
			fscanfln(steuerdatei, "%d", &this->SCN_ALB_writegrid, "write code snow canopy load grid for ", tstr);
			fscanfln(steuerdatei, "%s%d", &(this->SCN_ALB_filename[0]), &this->SCN_ALB_code, "name snow canopy load statistics file for ", tstr);
			fscanfln(steuerdatei, "%s", &(this->SCN_SUBL_gridname[0]), "name snow canopy load grid for ", tstr);
			fscanfln(steuerdatei, "%d", &this->SCN_SUBL_writegrid, "write code snow canopy load grid for ", tstr);
			fscanfln(steuerdatei, "%s%d", &(this->SCN_SUBL_filename[0]), &this->SCN_SUBL_code, "name snow canopy load statistics file for ", tstr);
			fscanfln(steuerdatei, "%s", &(this->SCN_Out_solid_gridname[0]), "name snow canopy unload (snow) grid for ", tstr);
			fscanfln(steuerdatei, "%d", &this->SCN_Out_solid_writegrid, "write code snow canopy unload (snow) grid for ", tstr);
			fscanfln(steuerdatei, "%s%d", &(this->SCN_Out_solid_filename[0]), &this->SCN_Out_solid_code, "name snow canopy unload (snow) statistics file for ", tstr);
			fscanfln(steuerdatei, "%s", &(this->SCN_Out_liq_gridname[0]), "name snow canopy liquid water output grid for ", tstr);
			fscanfln(steuerdatei, "%d", &this->SCN_Out_liq_writegrid, "write code snow canopy liquid water output grid for ", tstr);
			fscanfln(steuerdatei, "%s%d", &(this->SCN_Out_liq_filename[0]), &this->SCN_Out_liq_code, "name snow canopy liquid water output statistics file for ", tstr);
		}
		
		if (methode >= 21)
		{
			// model parameters for snow canopy model
			fscanfln(steuerdatei, "%lf", &this->SCN_min_lai, "minimum LAI value above which the snow canopy model will be started for", tstr);
			fscanfln(steuerdatei, "%lf", &this->SCN_min_z0, "minimum z0 value above which the snow canopy model will be started for", tstr);
			fscanfln(steuerdatei, "%lf", &this->SCN_factor_lai, "LAI scaling factor for the snow canopy model for", tstr);
		}

		if (methode >= 31) //ceb
		{
			run_ceb = 1;
			this->CEB_CANENBAL_gridname = newstring(_stringsize_ + 1, 1);
			this->CEB_CANTEMP_gridname = newstring(_stringsize_ + 1, 1);
			this->CEB_CANTEMP_OLD_gridname = newstring(_stringsize_ + 1, 1);
			this->CEB_DUMMY_gridname = newstring(_stringsize_ + 1, 1);

			this->CEB_CANENBAL_filename = newstring(_stringsize_ + 1, 1);
			this->CEB_CANTEMP_filename = newstring(_stringsize_ + 1, 1);
			this->CEB_CANTEMP_OLD_filename = newstring(_stringsize_ + 1, 1);
			this->CEB_DUMMY_filename = newstring(_stringsize_ + 1, 1);
			//grids for CEB
			fscanfln(steuerdatei, "%s", &(this->CEB_CANENBAL_gridname[0]), "name canopy energy balance grid for ", tstr);
			fscanfln(steuerdatei, "%d", &this->CEB_CANENBAL_writegrid, "write code canopy energy balance grid for ", tstr);
			fscanfln(steuerdatei, "%s%d", &(this->CEB_CANENBAL_filename[0]), &this->CEB_CANENBAL_code, "name canopy energy balance output statistics file for ", tstr);

			fscanfln(steuerdatei, "%s", &(this->CEB_CANTEMP_gridname[0]), "name canopy energy balance grid for ", tstr);
			fscanfln(steuerdatei, "%d", &this->CEB_CANTEMP_writegrid, "write code canopy energy balance grid for ", tstr);
			fscanfln(steuerdatei, "%s%d", &(this->CEB_CANTEMP_filename[0]), &this->CEB_CANTEMP_code, "name canopy energy balance output statistics file for ", tstr);

			fscanfln(steuerdatei, "%s", &(this->CEB_CANTEMP_OLD_gridname[0]), "name canopy energy balance grid for ", tstr);
			fscanfln(steuerdatei, "%d", &this->CEB_CANTEMP_OLD_writegrid, "write code canopy energy balance grid for ", tstr);
			fscanfln(steuerdatei, "%s%d", &(this->CEB_CANTEMP_OLD_filename[0]), &this->CEB_CANTEMP_OLD_code, "name canopy energy balance output statistics file for ", tstr);

			fscanfln(steuerdatei, "%s", &(this->CEB_DUMMY_gridname[0]), "name canopy energy balance grid for ", tstr);
			fscanfln(steuerdatei, "%d", &this->CEB_DUMMY_writegrid, "write code canopy energy balance grid for ", tstr);
			fscanfln(steuerdatei, "%s%d", &(this->CEB_DUMMY_filename[0]), &this->CEB_DUMMY_code, "name canopy energy balance output statistics file for ", tstr);
		}

		
//changes Kristian
		this->T0R = 0.;
		this->transzone = 0.;
		
		if(run_scn < 1){
			this->bSCN_model = false;
		}
		else{
			this->bSCN_model = true;
		}

		if(run_ceb < 1){
			this->CEB_model = false;
		}
		else{
			this->CEB_model = true;
		}

		for (lay=1;lay<=iNumVegetationLayers;lay++)
		{
			int p;
			std::string s;

			//		p=index(Qgridname[0],"."); 
			s=std::string(Qgridname[0]);
			p=s.find_last_of("."); 
			sprintf(Qgridname[lay],"%.*s_Layer_%d%s",p,Qgridname[0],lay,&(Qgridname[0][p]));

			//		p=index(EIgridname[0],"."); 
			s=std::string(EIgridname[0]);
			p=s.find_last_of("."); 
			sprintf(EIgridname[lay],"%.*s_Layer_%d%s",p,EIgridname[0],lay,&(EIgridname[0][p]));

			//		p=index(SIgridname[0],"."); 
			s=std::string(SIgridname[0]);
			p=s.find_last_of("."); 
			sprintf(SIgridname[lay],"%.*s_Layer_%d%s",p,SIgridname[0],lay,&(SIgridname[0][p]));

			//		p=index(Qstatfilename[0],"."); 
			s=std::string(Qstatfilename[0]);
			p=s.find_last_of("."); 
			sprintf(Qstatfilename[lay],"%.*s_Layer_%d%s",p,Qstatfilename[0],lay,&(Qstatfilename[0][p]));

			//		p=index(EIstatfilename[0],"."); 
			s=std::string(EIstatfilename[0]);
			p=s.find_last_of("."); 
			sprintf(EIstatfilename[lay],"%.*s_Layer_%d%s",p,EIstatfilename[0],lay,&(EIstatfilename[0][p]));

			//		p=index(SIstatfilename[0],"."); 
			s=std::string(SIstatfilename[0]);
			p=s.find_last_of("."); 
			sprintf(SIstatfilename[lay],"%.*s_Layer_%d%s",p,SIstatfilename[0],lay,&(SIstatfilename[0][p]));

		}

		if (fclose(steuerdatei)==EOF)
		{
			(*WSstream) << "FATAL: error closing "<<steuerdateiname<< "[end-entry]"; 
			throw(11);
		}


		// now reading the ArcVeg state (on or off)
		if ((steuerdatei = fopen(steuerdateiname, "rt")) == NULL)
		{
			(*WSstream) << "FATAL: error opening " << steuerdateiname << "[end-entry]";
			throw(28);
		}

		if (findeabschnitt(_ArcVegIntegration_, false)) // looks for a section  [ArcVegIntegration]
		{
			fscanfln(steuerdatei, "%d", &do_ArcVeg, "", tstr);

			if (do_ArcVeg > 0)
			{
				arcVegNumLayers = MultiLayerLUTable[1].getNumFractions();
				if (arcVegNumLayers != iNumVegetationLayers)
				{ // second try: assume the grids may be there but not the initial fractions in the multi layer landuse table
					// check for grids based on num layers
					arcVegNumLayers = iNumVegetationLayers;
				}

				bool bAllArcVegFractionsFound = true;
				if (arcVegNumLayers > 0)
				{
					char* layername;
					layername = new char[200];
					arcVegNumLayers = iNumVegetationLayers;
					ARCVEGFRACTIONS = new double** [arcVegNumLayers + 1];

					for (int layer = 1; layer <= arcVegNumLayers; layer++)
					{
						sprintf(layername, "%s%d", "ArcVeg_Fraction_", layer);
						if (ist_vorhanden(layername))
						{
							ARCVEGFRACTIONS[layer] = gridinfo[getindex(layername)].pointer;
						}
						else
						{
							bAllArcVegFractionsFound = false;
						}
					}

					delete[] layername;
				}
				else
				{
					bAllArcVegFractionsFound = false;
				}
				if (!bAllArcVegFractionsFound)
				{
					(*WSstream) << "FATAL: ERROR in initializing interception: [ArcVegIntegration] is active but there are no or not all fraction grids found as standard grids (ArcVeg_Fraction_1..n according to the number of land use layers)" << "[end-entry]";
					throw (1400);
				}
			}
		}

		if (fclose(steuerdatei) == EOF)
		{
			(*WSstream) << "FATAL: error closing " << steuerdateiname << "[end-entry]";
			throw(11);
		}

		/* zuordnen der Zeiger zu den stellvertretenden Variablen */
		char name[500];
		for (lay=1;lay<=iNumVegetationLayers;lay++)
		{
			sprintf(name,"%s%d",_Landnutzung_,lay);
			LUSE[lay]   = gridinfo[getindex(name)].pointer;

			sprintf(name,"%s%d",_Blattflaechenindex_,lay);
			LAI[lay]    = gridinfo[getindex(name)].pointer;

			sprintf(name,"%s%d",_Vegetationsbedeckungsgrad_,lay);
			VBG[lay]    = gridinfo[getindex(name)].pointer;

//js			// $kf 2017-06-09: get z0
			sprintf(name,"%s%d",_RoughnessLength_,lay);
			Z0[lay]     = gridinfo[getindex(name)].pointer;
			//$mk: get RSC + RSI
			sprintf(name,"%s%d",_SurfaceCanopyResistance_,lay);
			RSC[lay]    = gridinfo[getindex(name)].pointer;
			sprintf(name,"%s%d",_SurfaceInterzepResistance_,lay);
			RSI[lay]    = gridinfo[getindex(name)].pointer;
		}

		for (lay=0;lay<=iNumVegetationLayers;lay++)
		{
			sprintf(name,"%s%d",_Potentielle_Verdunstung_,lay);
			ETP[lay]    = gridinfo[getindex(name)].pointer;

			sprintf(name,"%s%d",_Potential_InterceptionEvap_,lay);
			EIP[lay]    = gridinfo[getindex(name)].pointer;

			sprintf(name,"%s%d",_Reale_Verdunstung_,lay);
			ETR[lay]    = gridinfo[getindex(name)].pointer;

		}

		if ((methode == 1) || (methode == 2)) // old traditional sequence snow -> interception
		{
			QInput = gridinfo[getindex_first_of_3(_Schneeausfluss_, _Niederschlag_, _ungueltiger_Name_)].pointer;
		}
		else // methode == 11 || methode == 12 || methode == 21 || methode == 22 || methode == 31 || methoe == 32  --> new sequence interception -> snow 
		{
			QInput = gridinfo[getindex( _Niederschlag_)].pointer; // here, we only look for precipitation from interpolation. There is no other sub model inbetween involved any more when using the new sequence
		}

		// $kf: 2016-06-15: Define additional grids that are not required for standard interception calculation
		// update 2016-08-22: read values also in case the snow model isnot activated since the grids must contain
		// values for the snow model
		if (bSCN_model)
		{
			this->TInput = gridinfo[this->getindex_first_of_3(_Korrigierte_Temperatur_, _Temperatur_, _Temperatur_14_)].pointer;
			this->RInput = gridinfo[this->getindex(_AspectCorrectedGlobalRadiation_)].pointer;
			this->WSInput = gridinfo[getindex(_Windgeschwindigkeit_)].pointer;
			this->RHInput = gridinfo[getindex(_Luftfeuchte_)].pointer;
			this->CLInput = gridinfo[getindex(_Cloudiness_)].pointer;
			this->TempGrids = new GridStack(int(1440 / InV), mpistartrow, mpiendrow, nrows, ncols);
		}
		if (CEB_model)
		{
			if (ist_vorhanden(_Albedo_))
			{
				this->Albedo = gridinfo[getindex(_Albedo_)].pointer;
			}
			else
			{
				(*WSstream) << "FATAL: Albedo not available in advanced snow interception model with canopy energy balance (should be available when an extended snow model is used)" << "[end-entry]";
				throw(1160);
			}

			if (ist_vorhanden(_Schneespeicher_Total_))
			{
				this->SSTO = gridinfo[getindex(_Schneespeicher_Total_)].pointer;
			}
			else
			{
				(*WSstream) << "FATAL: total snow storage not available in advanced snow interception model with canopy energy balance (use any snow model to have this grid been created for you)" << "[end-entry]";
				throw(1160);
			}

			if (ist_vorhanden(_surftemp_))
			{
				SURFTEMP = gridinfo[getindex(_surftemp_)].pointer;
			}
			else
			{
				(*WSstream) << "FATAL: snow surface temperatur not available in advanced snow interception model with canopy energy balance (in snow model, an energy balance approach should be used)" << "[end-entry]";
				throw(1160);
			}
		}

		/* Alias fuer das Hoehenmodell als Zeiger referenzieren */
		dhm = gridinfo[getindex(_Hoehenmodell_)].pointer;

		/* wo ist Platz fuer SI  im array gridinfo? suchen des ersten freien Zeigers */
		i=1;
		while (gridinfo[i].identification!=NULL)
		{
			i++;
			if (i==gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no free grid left"<<gridmaxanzahl << "[end-entry]";
				throw(12);
			}
		}

		int ii;
		char layername[2000];
		for (lay=0;lay<=iNumVegetationLayers;lay++)
		{
			sprintf(layername,"%s%d",_Interzeptionsspeicher_,lay);
			if (gridslesen!=1)
			/* Grid wird nicht von Platte eingelesen */
			{
				SI[lay]  = matrix(mpistartrow,mpiendrow,1,ncols,nrows);
				initialize_grid(i, SIgridname[lay], SI[lay], NULL, 0, layername, SIwritegrid, mpistartrow, mpiendrow);
				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int j = 1; j <= ncols; j++)
					{
						SI[lay][k][j] = zn[k][j] == nodata ? nodata : 0;
					}
				}
				siindex[lay] = i;
			}
			else
			{
				/* alte Dateien einlesen */
				readgridfromfile(SIgridname[lay], &(SI[lay]), 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
				initialize_grid(i, SIgridname[lay], SI[lay], NULL, 0, layername, SIwritegrid, mpistartrow, mpiendrow);
				siindex[lay] = i;
			}

			i++;

			// when running lake model, EI of layer 0 is required as full grid in rank 0
			int srei0 = (mpi_rank == 0) ? 1 : mpistartrow;
			int erei0 = (mpi_rank == 0) ? nrows : mpiendrow;

			EI[lay] = matrix(srei0, erei0, 1, ncols, nrows);
			sprintf(layername,"%s%d",_Interzeptionsverdunstung_,lay);
			initialize_grid(i, EIgridname[lay], EI[lay], NULL, 0, layername, EIwritegrid, srei0, erei0);
			for (int k = srei0; k <= erei0; k++)
			{
				for (int j = 1; j <= ncols; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						EI[lay][k][j] = 0;
					}
					else
					{
						EI[lay][k][j] = nodata;
					}
				}
			}
			eiindex[lay] = i;

			i++;
			sprintf(layername,"%s%d",_Bestandesniederschlag_,lay);
			Q[lay] = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(i, Qgridname[lay], Q[lay], NULL, 0, layername, Qwritegrid, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = 1; j <= ncols; j++)
				{
					Q[lay][k][j] = zn[k][j] == nodata ? nodata : 0;
				}
			}
			qindex[lay] = i;


			zoneEI[lay] = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zoneEI[lay]) speichererror(1);
			for (ii=0;ii<=zonmax;ii++)
			{
				zoneEI[lay][ii].n          = zoneEI[lay][ii].sx         = 0;
				zoneEI[lay][ii].zonenwert  = zone[ii].zonenwert;
				zoneEI[lay][ii].rel        = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				zoneEI[lay][ii].RelCellSum = zone[ii].RelCellSum;
			}
			
			sprintf(layername,"%s%d_","interception_evaporation_in_mm_layer_",lay);
			EIstat[lay] = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
			&(zoneEI[lay]),(int)zonmax,(int)EIstatcode,(int)5,(char*)EIstatfilename[lay],DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
			if (retvalue!=0) throw retvalue;
			EIstat[lay]->schreibzonenwerte(layername,instanzname,LastResults,maxLastResults);

			zoneSI[lay] = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zoneSI[lay]) speichererror(1);
			for (ii=0;ii<=zonmax;ii++)
			{
				zoneSI[lay][ii].n          = zoneSI[lay][ii].sx         = 0;
				zoneSI[lay][ii].zonenwert  = zone[ii].zonenwert;
				zoneSI[lay][ii].rel        = zone[ii].rel;  /* fuer die relativen AntSIle der Zonen am (gueltigen) Gesamtgrid */
				zoneSI[lay][ii].RelCellSum = zone[ii].RelCellSum;
			}
			
			sprintf(layername,"%s%d_","interception_storage_in_mm_layer_",lay);
			SIstat[lay] = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
			&(zoneSI[lay]),(int)zonmax,(int)SIstatcode,(int)5,(char*)SIstatfilename[lay],DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
			if (retvalue!=0) throw retvalue;
			SIstat[lay]->schreibzonenwerte(layername,instanzname,LastResults,maxLastResults);

			zoneQ[lay] = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
			if (!zoneQ[lay]) speichererror(1);
			for (ii=0;ii<=zonmax;ii++)
			{
				zoneQ[lay][ii].n          = zoneQ[lay][ii].sx         = 0;
				zoneQ[lay][ii].zonenwert  = zone[ii].zonenwert;
				zoneQ[lay][ii].rel        = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				zoneQ[lay][ii].RelCellSum = zone[ii].RelCellSum;
			}
			sprintf(layername,"%s%d_","throughfall_in_mm_layer_",lay);
			Qstat[lay]  = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
			&(zoneQ[lay]),(int)zonmax,(int)Qstatcode,(int)5,(char*)Qstatfilename[lay],DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
			if (retvalue!=0) throw retvalue;
			Qstat[lay]->schreibzonenwerte(layername,instanzname,LastResults,maxLastResults);
			i++;
		}
		
		if (this->bSCN_model)
		{
			// $kf 2016-06-20: Setup grids for the state variables of the snow canopy model
			/* wo ist Platz fuer Snow Canopy Grids im array gridinfo? suchen des ersten freien Zeigers */
			i = 1;
			while (gridinfo[i].identification != NULL)
			{
				i++;
				if (i + 14 >= gridmaxanzahl)
				{
					(*WSstream) << "FATAL: no free grid left" << gridmaxanzahl << "[end-entry]";
					throw(12);
				}
			}
			// define grids for output, @todo: implement reading as well
			// @todo: link definitions to control file
			// for testing aeach of the grids listed below will be written to a separate file
			this->grd_idx_scn_load = i + 0;
			this->grd_idx_scn_alb = i + 1;
			this->grd_idx_scn_age = i + 2;
			this->grd_idx_scn_subl = i + 3;
			this->grd_idx_scn_out_s = i + 4;
			this->grd_idx_scn_out_l = i + 5;
			this->grd_idx_scn_temp = i + 6;
			this->grd_idx_scn_rh = i + 7;
			this->grd_idx_scn_rad = i + 8;
			this->grd_idx_scn_wind = i + 9;
			this->grd_idx_scn_ir = i + 10;
			if (CEB_model)
			{
				this->grd_idx_ceb_canenbal = i + 11;
				this->grd_idx_ceb_cantemp = i + 12;
				this->grd_idx_ceb_cantemp_old = i + 13;
				this->grd_idx_ceb_dummy = i + 14;
			}


			// first attempt: initialize the grids using default values
			if (gridslesen != 1)
			{
				this->SCN_LOAD = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				this->SCN_AGE = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				this->SCN_Albedo = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			}
			else
			{
				readgridfromfile(SCN_LOAD_gridname, &(SCN_LOAD), 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
				readgridfromfile(SCN_ALB_gridname, &(SCN_Albedo), 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
				readgridfromfile(SCN_AGE_gridname, &(SCN_AGE), 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
			}
			this->SCN_SUBL = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			this->SCN_Out_solid = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			this->SCN_Out_liquid = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			this->SCN_Temp = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			this->SCN_RH = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			this->SCN_RAD = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			this->SCN_Wind = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			this->SCN_IR = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			//ceb
			if (CEB_model)
			{
				this->CEB_CANENBAL = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				this->CEB_CANTEMP = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				if (gridslesen != 1)
				{
					this->CEB_CANTEMP_OLD = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
				}
				else
				{
					readgridfromfile(CEB_CANTEMP_OLD_gridname, &(CEB_CANTEMP_OLD), 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
				}
				this->CEB_DUMMY = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			}
			// $kf 2016-07-08: @todo replace SIwritegrid by another variable dedicated to each grid
			initialize_grid(this->grd_idx_scn_load, this->SCN_LOAD_gridname, this->SCN_LOAD, NULL, 0, _Snow_Canopy_Load_, this->SCN_LOAD_writegrid, mpistartrow, mpiendrow);
			initialize_grid(this->grd_idx_scn_alb, this->SCN_ALB_gridname, this->SCN_Albedo, NULL, 0, _Snow_Canopy_Albedo_, this->SCN_ALB_writegrid, mpistartrow, mpiendrow);
			initialize_grid(this->grd_idx_scn_age, this->SCN_AGE_gridname, this->SCN_AGE, NULL, 0, _Snow_Canopy_Age_, this->SCN_AGE_writegrid, mpistartrow, mpiendrow);
			initialize_grid(this->grd_idx_scn_subl, this->SCN_SUBL_gridname, this->SCN_SUBL, NULL, 0, _Snow_Canopy_Sublim_, this->SCN_SUBL_writegrid, mpistartrow, mpiendrow);
			initialize_grid(this->grd_idx_scn_out_s, this->SCN_Out_solid_gridname, this->SCN_Out_solid, NULL, 0, _Snow_Canopy_Unload_solid_, this->SCN_Out_solid_writegrid, mpistartrow, mpiendrow);
			initialize_grid(this->grd_idx_scn_out_l, this->SCN_Out_liq_gridname, this->SCN_Out_liquid, NULL, 0, _Snow_Canopy_Unload_liquid_, this->SCN_Out_liq_writegrid, mpistartrow, mpiendrow);
			initialize_grid(this->grd_idx_scn_temp, this->SCN_Temp_gridname, this->SCN_Temp, NULL, 0, _Snow_Canopy_Temperature_, this->SCN_Temp_writegrid, mpistartrow, mpiendrow);
			initialize_grid(this->grd_idx_scn_rh, this->SCN_RH_gridname, this->SCN_RH, NULL, 0, _Snow_Canopy_RH_, this->SCN_RH_writegrid, mpistartrow, mpiendrow);
			initialize_grid(this->grd_idx_scn_rad, this->SCN_RAD_gridname, this->SCN_RAD, NULL, 0, _Snow_Canopy_RAD_, this->SCN_RAD_writegrid, mpistartrow, mpiendrow);
			initialize_grid(this->grd_idx_scn_wind, this->SCN_Wind_gridname, this->SCN_Wind, NULL, 0, _Snow_Canopy_Wind_, this->SCN_Wind_writegrid, mpistartrow, mpiendrow);
			initialize_grid(this->grd_idx_scn_ir, "SCN_IR", this->SCN_IR, NULL, 0, _Snow_Canopy_IR_, 0 /* do not write grid as it is an intermediate result */, mpistartrow, mpiendrow);
			//ceb
			if (CEB_model)
			{

				initialize_grid(this->grd_idx_ceb_canenbal, this->CEB_CANENBAL_gridname, this->CEB_CANENBAL, NULL, 0, _Canopy_Energy_Balance_, this->CEB_CANENBAL_writegrid, mpistartrow, mpiendrow);
				initialize_grid(this->grd_idx_ceb_cantemp, this->CEB_CANTEMP_gridname, this->CEB_CANTEMP, NULL, 0, _Canopy_Temperature_, this->CEB_CANTEMP_writegrid, mpistartrow, mpiendrow);
				initialize_grid(this->grd_idx_ceb_cantemp_old, this->CEB_CANTEMP_OLD_gridname, this->CEB_CANTEMP_OLD, NULL, 0, _Canopy_Temperature_Old_, this->CEB_CANTEMP_OLD_writegrid, mpistartrow, mpiendrow);
				initialize_grid(this->grd_idx_ceb_dummy, this->CEB_DUMMY_gridname, this->CEB_DUMMY, NULL, 0, _Canopy_Dummy_, this->CEB_DUMMY_writegrid, mpistartrow, mpiendrow);
			}
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = 1; j <= ncols; j++)
				{

					if (gridslesen != 1)
					{ // otherwise grids are already read in
						this->SCN_LOAD[k][j] = zn[k][j] == nodata ? nodata : 0;
						this->SCN_Albedo[k][j] = zn[k][j] == nodata ? nodata : 0.2;
						this->SCN_AGE[k][j] = zn[k][j] == nodata ? nodata : 0;
					}
					this->SCN_SUBL[k][j] = zn[k][j] == nodata ? nodata : 0;
					this->SCN_Out_solid[k][j] = zn[k][j] == nodata ? nodata : 0;
					this->SCN_Out_liquid[k][j] = zn[k][j] == nodata ? nodata : 0;
					this->SCN_Temp[k][j] = zn[k][j] == nodata ? nodata : 0;
					this->SCN_RH[k][j] = zn[k][j] == nodata ? nodata : 0;
					this->SCN_RAD[k][j] = zn[k][j] == nodata ? nodata : 0;
					this->SCN_Wind[k][j] = zn[k][j] == nodata ? nodata : 0;
					this->SCN_IR[k][j] = zn[k][j] == nodata ? nodata : 0;
					//ceb
					if (CEB_model)
					{
						this->CEB_CANENBAL[k][j] = zn[k][j] == nodata ? nodata : 0;
						this->CEB_CANTEMP[k][j] = zn[k][j] == nodata ? nodata : 0;
						if (gridslesen != 1)
						{
							this->CEB_CANTEMP_OLD[k][j] = zn[k][j] == nodata ? nodata : 0;
						}
						this->CEB_DUMMY[k][j] = zn[k][j] == nodata ? nodata : 0;
					}
				}
			}

			// $kf 2016-07-06: zonal statistics for snow canopy interactions
			this->SCN_Load_zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
			this->SCN_Age_zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
			this->SCN_Alb_zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
			this->SCN_Sub_zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
			this->SCN_Out_solid_zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
			this->SCN_Out_liq_zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
			this->SCN_Temp_zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
			this->SCN_RH_zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
			this->SCN_RAD_zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
			this->SCN_Wind_zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
			//ceb
			if (CEB_model)
			{
				this->CEB_CANENBAL_zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
				this->CEB_CANTEMP_zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
				this->CEB_CANTEMP_OLD_zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
				this->CEB_DUMMY_zone = (statrecord *)calloc(zonmax + 1, sizeof(statrecord));
			}
			if (!this->SCN_Load_zone) speichererror(1);
			if (!this->SCN_Age_zone) speichererror(1);
			if (!this->SCN_Alb_zone) speichererror(1);
			if (!this->SCN_Sub_zone) speichererror(1);
			if (!this->SCN_Out_solid_zone) speichererror(1);
			if (!this->SCN_Out_liq_zone) speichererror(1);
			if (!this->SCN_Temp_zone) speichererror(1);
			if (!this->SCN_RH_zone) speichererror(1);
			if (!this->SCN_RAD_zone) speichererror(1);
			if (!this->SCN_Wind_zone) speichererror(1);
			//ceb
			if (CEB_model)
			{

				if (!this->CEB_CANENBAL_zone) speichererror(1);
				if (!this->CEB_CANTEMP_zone) speichererror(1);
				if (!this->CEB_CANTEMP_OLD_zone) speichererror(1);
				if (!this->CEB_DUMMY_zone) speichererror(1);
			}
			for (ii = 0; ii <= zonmax; ii++)
			{
				this->SCN_Load_zone[ii].n = this->SCN_Load_zone[ii].sx = 0;
				this->SCN_Load_zone[ii].zonenwert = zone[ii].zonenwert;
				this->SCN_Load_zone[ii].rel = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				this->SCN_Load_zone[ii].RelCellSum = zone[ii].RelCellSum;

				this->SCN_Age_zone[ii].n = this->SCN_Age_zone[ii].sx = 0;
				this->SCN_Age_zone[ii].zonenwert = zone[ii].zonenwert;
				this->SCN_Age_zone[ii].rel = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				this->SCN_Age_zone[ii].RelCellSum = zone[ii].RelCellSum;

				this->SCN_Alb_zone[ii].n = this->SCN_Alb_zone[ii].sx = 0;
				this->SCN_Alb_zone[ii].zonenwert = zone[ii].zonenwert;
				this->SCN_Alb_zone[ii].rel = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				this->SCN_Alb_zone[ii].RelCellSum = zone[ii].RelCellSum;

				this->SCN_Sub_zone[ii].n = this->SCN_Sub_zone[ii].sx = 0;
				this->SCN_Sub_zone[ii].zonenwert = zone[ii].zonenwert;
				this->SCN_Sub_zone[ii].rel = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				this->SCN_Sub_zone[ii].RelCellSum = zone[ii].RelCellSum;

				this->SCN_Out_solid_zone[ii].n = this->SCN_Out_solid_zone[ii].sx = 0;
				this->SCN_Out_solid_zone[ii].zonenwert = zone[ii].zonenwert;
				this->SCN_Out_solid_zone[ii].rel = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				this->SCN_Out_solid_zone[ii].RelCellSum = zone[ii].RelCellSum;

				this->SCN_Out_liq_zone[ii].n = this->SCN_Out_liq_zone[ii].sx = 0;
				this->SCN_Out_liq_zone[ii].zonenwert = zone[ii].zonenwert;
				this->SCN_Out_liq_zone[ii].rel = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				this->SCN_Out_liq_zone[ii].RelCellSum = zone[ii].RelCellSum;

				this->SCN_Temp_zone[ii].n = this->SCN_Temp_zone[ii].sx = 0;
				this->SCN_Temp_zone[ii].zonenwert = zone[ii].zonenwert;
				this->SCN_Temp_zone[ii].rel = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				this->SCN_Temp_zone[ii].RelCellSum = zone[ii].RelCellSum;

				this->SCN_RH_zone[ii].n = this->SCN_RH_zone[ii].sx = 0;
				this->SCN_RH_zone[ii].zonenwert = zone[ii].zonenwert;
				this->SCN_RH_zone[ii].rel = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				this->SCN_RH_zone[ii].RelCellSum = zone[ii].RelCellSum;

				this->SCN_RAD_zone[ii].n = this->SCN_RAD_zone[ii].sx = 0;
				this->SCN_RAD_zone[ii].zonenwert = zone[ii].zonenwert;
				this->SCN_RAD_zone[ii].rel = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				this->SCN_RAD_zone[ii].RelCellSum = zone[ii].RelCellSum;

				this->SCN_Wind_zone[ii].n = this->SCN_Wind_zone[ii].sx = 0;
				this->SCN_Wind_zone[ii].zonenwert = zone[ii].zonenwert;
				this->SCN_Wind_zone[ii].rel = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
				this->SCN_Wind_zone[ii].RelCellSum = zone[ii].RelCellSum;

				if (CEB_model)
				{
					this->CEB_CANENBAL_zone[ii].n = this->CEB_CANENBAL_zone[ii].sx = 0;
					this->CEB_CANENBAL_zone[ii].zonenwert = zone[ii].zonenwert;
					this->CEB_CANENBAL_zone[ii].rel = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
					this->CEB_CANENBAL_zone[ii].RelCellSum = zone[ii].RelCellSum;

					this->CEB_CANTEMP_zone[ii].n = this->CEB_CANTEMP_zone[ii].sx = 0;
					this->CEB_CANTEMP_zone[ii].zonenwert = zone[ii].zonenwert;
					this->CEB_CANTEMP_zone[ii].rel = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
					this->CEB_CANTEMP_zone[ii].RelCellSum = zone[ii].RelCellSum;

					this->CEB_CANTEMP_OLD_zone[ii].n = this->CEB_CANTEMP_OLD_zone[ii].sx = 0;
					this->CEB_CANTEMP_OLD_zone[ii].zonenwert = zone[ii].zonenwert;
					this->CEB_CANTEMP_OLD_zone[ii].rel = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
					this->CEB_CANTEMP_OLD_zone[ii].RelCellSum = zone[ii].RelCellSum;

					this->CEB_DUMMY_zone[ii].n = this->CEB_DUMMY_zone[ii].sx = 0;
					this->CEB_DUMMY_zone[ii].zonenwert = zone[ii].zonenwert;
					this->CEB_DUMMY_zone[ii].rel = zone[ii].rel;  /* fuer die relativen Anteile der Zonen am (gueltigen) Gesamtgrid */
					this->CEB_DUMMY_zone[ii].RelCellSum = zone[ii].RelCellSum;
				}
			}

			this->SCN_Load_stats = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols,
				&(this->SCN_Load_zone), (int)zonmax, (int) this->SCN_LOAD_code, (int)5, (char*) this->SCN_LOAD_filename, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
			if (retvalue != 0) throw retvalue;
			this->SCN_Load_stats->schreibzonenwerte("SCN_LOAD_in_mm", instanzname, LastResults, maxLastResults);

			this->SCN_Age_stats = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols,
				&(this->SCN_Age_zone), (int)zonmax, (int) this->SCN_AGE_code, (int)5, (char*) this->SCN_AGE_filename, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
			if (retvalue != 0) throw retvalue;
			this->SCN_Age_stats->schreibzonenwerte("SCN_AGE_Number_of_time_steps_per_day", instanzname, LastResults, maxLastResults);

			this->SCN_Alb_stats = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols,
				&(this->SCN_Alb_zone), (int)zonmax, (int) this->SCN_ALB_code, (int)5, (char*) this->SCN_ALB_filename, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
			if (retvalue != 0) throw retvalue;
			this->SCN_Alb_stats->schreibzonenwerte("SCN_ALBEDO_Dimensionless", instanzname, LastResults, maxLastResults);

			this->SCN_Sub_stats = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols,
				&(this->SCN_Sub_zone), (int)zonmax, (int) this->SCN_SUBL_code, (int)5, (char*) this->SCN_SUBL_filename, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
			if (retvalue != 0) throw retvalue;
			this->SCN_Sub_stats->schreibzonenwerte("SCN_SUBLIMATION_mm_per_time_step", instanzname, LastResults, maxLastResults);

			this->SCN_Out_solid_stats = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols,
				&(this->SCN_Out_solid_zone), (int)zonmax, (int) this->SCN_Out_solid_code, (int)5, (char*) this->SCN_Out_solid_filename, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
			if (retvalue != 0) throw retvalue;
			this->SCN_Out_solid_stats->schreibzonenwerte("SCN_SNOW_Throughfall_mm_per_time_step", instanzname, LastResults, maxLastResults);

			this->SCN_Out_liq_stats = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols,
				&(this->SCN_Out_liq_zone), (int)zonmax, (int) this->SCN_Out_liq_code, (int)5, (char*) this->SCN_Out_liq_filename, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
			if (retvalue != 0) throw retvalue;
			this->SCN_Out_liq_stats->schreibzonenwerte("SCN_Output_rain_mm_per_time_step", instanzname, LastResults, maxLastResults);

			this->SCN_Temp_stats = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols,
				&(this->SCN_Temp_zone), (int)zonmax, (int) this->SCN_Temp_code, (int)5, (char*) this->SCN_Temp_filename, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
			if (retvalue != 0) throw retvalue;
			this->SCN_Temp_stats->schreibzonenwerte("SCN_sub-canopy_temperature", instanzname, LastResults, maxLastResults);

			this->SCN_RH_stats = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols,
				&(this->SCN_RH_zone), (int)zonmax, (int) this->SCN_RH_code, (int)5, (char*) this->SCN_RH_filename, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
			if (retvalue != 0) throw retvalue;
			this->SCN_RH_stats->schreibzonenwerte("SCN_sub-canopy_relative_humidity", instanzname, LastResults, maxLastResults);

			this->SCN_RAD_stats = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols,
				&(this->SCN_RAD_zone), (int)zonmax, (int) this->SCN_RAD_code, (int)5, (char*) this->SCN_RAD_filename, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
			if (retvalue != 0) throw retvalue;
			this->SCN_RAD_stats->schreibzonenwerte("SCN_sub-canopy_shortwave_radiation", instanzname, LastResults, maxLastResults);

			this->SCN_Wind_stats = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols,
				&(this->SCN_Wind_zone), (int)zonmax, (int) this->SCN_Wind_code, (int)5, (char*) this->SCN_Wind_filename, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
			if (retvalue != 0) throw retvalue;
			this->SCN_Wind_stats->schreibzonenwerte("SCN_sub-canopy_windspeed", instanzname, LastResults, maxLastResults);

			if (CEB_model)
			{
				this->CEB_CANENBAL_stats = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols,
					&(this->CEB_CANENBAL_zone), (int)zonmax, (int) this->CEB_CANENBAL_code, (int)5, (char*) this->CEB_CANENBAL_filename, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
				if (retvalue != 0) throw retvalue;
				this->CEB_CANENBAL_stats->schreibzonenwerte("CEB_canopy_energy_balance", instanzname, LastResults, maxLastResults);

				this->CEB_CANTEMP_stats = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols,
					&(this->CEB_CANTEMP_zone), (int)zonmax, (int) this->CEB_CANTEMP_code, (int)5, (char*) this->CEB_CANTEMP_filename, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
				if (retvalue != 0) throw retvalue;
				this->CEB_CANTEMP_stats->schreibzonenwerte("CEB_canopy_temperature", instanzname, LastResults, maxLastResults);

				this->CEB_CANTEMP_OLD_stats = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols,
					&(this->CEB_CANTEMP_OLD_zone), (int)zonmax, (int) this->CEB_CANTEMP_OLD_code, (int)5, (char*) this->CEB_CANTEMP_OLD_filename, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
				if (retvalue != 0) throw retvalue;
				this->CEB_CANTEMP_OLD_stats->schreibzonenwerte("CEB_canopy_temperature_old", instanzname, LastResults, maxLastResults);

				this->CEB_DUMMY_stats = new zonalstatistic(retvalue, WSstream, LastResults, maxLastResults, RelCellArea, nrows, ncols,
					&(this->CEB_DUMMY_zone), (int)zonmax, (int) this->CEB_DUMMY_code, (int)5, (char*) this->CEB_DUMMY_filename, DefaultOutputDirectory, (int)numsubb, (int*)listsubb, (InV < 60));
				if (retvalue != 0) throw retvalue;
				this->CEB_DUMMY_stats->schreibzonenwerte("CEB_dummy", instanzname, LastResults, maxLastResults);
			}
		}
			
		/* looking for ICE grid */
		if (ist_vorhanden(_ice_firn_)||ist_vorhanden(_GlacierizedCells_))
		{
			if (ist_vorhanden(_GlacierizedCells_))
				ICE=gridinfo[getindex(_GlacierizedCells_)].pointer;
			else
				ICE=gridinfo[getindex(_ice_firn_)].pointer;
		}
		else
		{
			ICE=NULL;
		}

		if (doSpinUp!=0 && HRU==NULL)
		{
			HRU       = this->gridinfo[getindex(_HydrologicResponseUnits_)].pointer;
			HRUMaster = this->gridinfo[getindex(_MasterHydrologicResponseUnits_)].pointer;
		}


	}
	catch (int ex)
	{
		retvalue=ex;
	}
	catch (...)
	// catching all exceptions in order to allow the object to be constructed correct. This is required 
	// in order to allow a dextructor call in case of an error during construction
	{
		retvalue=103;
	}
}

void interzeptionsmodellklasse::run(long me,long se,long tg,long mt,long jr)
{
	int i,j,k,z;
	double icevalue;
	double trms,trmsmax; /* for non-evaporating tracers: tracer-mass, max. allowed tracer mass (dependent on max. concentration) */
	if (do_it == 0) return;
	minute = me;
	stunde = se;
	tag    = tg;
	monat  = mt;
	jahr   = jr;
	teilmodellklasse::run(me,se,tg,mt,jr,NULL,NULL);

	tagnr=0;
	for (k=1;k<monat;k++) tagnr += mon[k];
	tagnr += tag;

	int lay;
	for (lay=0;lay<=iNumVegetationLayers;lay++)
	{
		EIstat[lay]->init();
		SIstat[lay]->init();
		Qstat[lay]->init();
	}
	for (z=1;z<=tracer_count;z++) tracer_module[z].stat->init();

	// $kf 2016-06-20 Temperature stack to calculate moving averages:
	// An average temperature representing a daily mean temperature
	// is required to compute the damped diurnal temperature course
	// in the canopy. This task is performed using the GridStack
	// class which has a memory of the last 24 grids.
	if (methode > 20) // (bSCN_model==true && CEB_model==true)
	{
		if (stunde != stunde_old || tag != tag_old || monat_old != monat || jahr_old != jahr)
		{
			// add temperature to stack
			double** temp_grid = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			for (int i = mpistartrow; i < mpiendrow; i++)
			{
				for (int j = this->FirstLastCol[i].First; j <= this->FirstLastCol[i].Last; j++)
				{
					if (MPIREG[i][j] == mpi_rank)
					{
						temp_grid[i][j] = TInput[i][j];
					}
				}
			}

			this->TempGrids->addToStack(temp_grid);
		}
		this->SCN_Load_stats->init();
		this->SCN_Age_stats->init();
		this->SCN_Alb_stats->init();
		this->SCN_Sub_stats->init();
		this->SCN_Out_solid_stats->init();
		this->SCN_Out_liq_stats->init();
		this->SCN_Temp_stats->init();
		this->SCN_RH_stats->init();
		this->SCN_RAD_stats->init();
		this->SCN_Wind_stats->init();
		if (CEB_model)
		{
			this->CEB_CANENBAL_stats->init();
			this->CEB_CANTEMP_stats->init();
			this->CEB_CANTEMP_OLD_stats->init();
			this->CEB_DUMMY_stats->init();
		}
	}

	for (i=mpistartrow;i<=mpiendrow;i++)
	{
		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
		{
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
			//      if (zn[i][j]!=nodata)
			{
				icevalue=ICE?ICE[i][j]:0; // valid range of icevalue now between 0...2
				if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
				else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
				if (icevalue<1.0) // performance is better when we can skip this section for pure ice cells 
				{
					double input_to_rainwater_interc = this->QInput[i][j];
					double lai = 0.;
					double z0 = 0.;
					double rsc = 3000;//high value
					double rsi = 3000;//high value
					double si = 0;
					for (lay=1;lay<=iNumVegetationLayers;lay++)
					{
						lai += LAI[lay][i][j];
						if(Z0[lay][i][j] > z0)
							z0 = Z0[lay][i][j]; // get maximum z0 of land use
						if(RSC[lay][i][j] < rsc)
							rsc = RSC[lay][i][j]; // mk: get stomata resistance
						if(RSI[lay][i][j] < rsi)
							rsi = RSI[lay][i][j]; // mk: get interception resistance 
						if(SI[lay][i][j] > si)
							si = SI[lay][i][j]; // mk: get intercepted water
						
					}
					if (methode > 20 && this->bSCN_model)
					{
						if (lai > this->SCN_min_lai && z0 > this->SCN_min_z0)
						{
							// $kf: 2016-06-15: compile meteo data
							// create SnowCanopyParam object for data exchange
							SnowCanopyParam* p = new SnowCanopyParam();
							// (1) Meteorological forcing data
							p->PRECIP = this->QInput[i][j];
							p->AIRTEMP = this->TInput[i][j] + 273.16;
							p->MEANTEMP = this->TempGrids->getAvg(i, j) + 273.16;
							p->HUM = this->RHInput[i][j] * 100.;
							p->SHORTWAVE = this->RInput[i][j];
							p->WIND = this->WSInput[i][j];
							p->CL = this->CLInput[i][j];
							p->DELTATEMP = (p->MEANTEMP - 273.16) / 3;
							// Korrektur von DELTATEMP, wie in SnowMod.f90
							if (p->DELTATEMP < -2.0)
								p->DELTATEMP = -2.0;
							else if (p->DELTATEMP > +2.0)
								p->DELTATEMP = +2.0;
							// (2) Parameters @todo: add albedo parameters (at present, default values are used)
							// canopy height is not used as the values cancel in the algorithm
							// @todo: calibration values for T-INDEX (Pelliciotti)
							p->INTERVALL = this->InV * 60;
							p->ALTITUDE = this->dhm[i][j];
							p->LAI = /* LAI[1][i][j]*/ lai * this->SCN_factor_lai;
							// (3) State variables
							p->ALBEDO = this->SCN_Albedo[i][j];
							p->CANOPINTLOADOLD = this->SCN_LOAD[i][j];
							p->SNOWAGE = this->SCN_AGE[i][j];

							// simple snow processes representation in the canopy
							// Phase detection
							// PhasePartitioning(p->AIRTEMP, p->HUM, p->PRECIP, p->ALTITUDE, p->SNOWPRECIP, p->RAINPRECIP);

							// $kf 2018-06-15: in order to avoid two different approaches for separating rain and snow in one model,
							// we use a simple threshold with transition zone (copied from snow model)
							// todo: define a function that does this!!!
							PhasePartitioningWaSiM(this->TInput[i][j], this->T0R, this->transzone, p->PRECIP, p->SNOWPRECIP, p->RAINPRECIP);

							// for mass balance check; remove later
							double last_ci = p->CANOPINTLOADOLD;

							// Call snow canpy model (if snowfall occurs or if a snow layer is prevailing on the trees):
							// Update: call the model for each time step to achieve consistent sub-canopy variables
							// ************************************************
							//if(p->SNOWPRECIP > 0. || p->CANOPINTLOADOLD > 0.0)
							//{
							CANOPY_INTERCEPTION(p);
							/* If the model has been called, set liquid water outflow as input
							 * for the subsequent rain water interception. Snow unload is not
							 * involved in any further interception computations. Instead, it
							 * is redirected to the snow model where it is acummulated to the
							 * snowpack (as solid input).*/
							input_to_rainwater_interc = p->RAINPRECIP;
							//}
							
							if (CEB_model)
							{
								//get Parameters from other Models
								double SNOWSTOR = this->SSTO[i][j]; //Include pointer!!!! --> Working!! --> with bisection Iteration no longer necessary!
								p->SURFTEMP = this->SURFTEMP[i][j]; //LW radiation
								p->GLOBRAD = this->RInput[i][j]; //SW Input
								p->ALBEDO_SURF = this->Albedo[i][j]; //SW --> CORRECT??
								p->ALBEDO_CAN = this->SCN_Albedo[i][j];//SW
								p->CANTEMP = p->AIRTEMP + 0.5;//this->CEB_CANTEMP[i][j];//BM ??
								p->CANTEMP_OLD = (this->CEB_CANTEMP_OLD[i][j] == 0) ? p->CANTEMP : this->CEB_CANTEMP_OLD[i][j];//BM ??
								//get leaf resistance + interception resistance
								p->RSC = rsc;
								p->RSI = rsi;
								//surface roughness from cntrl file
								p->Z0 = z0;
								//get intercepted water
								p->SI = si;
								p->SI_MAX = this->SImax;

								//DETERMINE SURFACE ROUGHNESS 
								
								if (SNOWSTOR > 0.0)
								{
									SOLVE_CEB(p);

									//CEB
									this->CEB_CANENBAL[i][j] = p->CANENBAL;
									this->CEB_CANTEMP[i][j] = p->CANTEMP;
									this->CEB_CANTEMP_OLD[i][j] = p->CANTEMP;
									this->CEB_DUMMY[i][j] = p->DUMMY;
								}
								else
								{
									this->CEB_CANENBAL[i][j] = 0;
									this->CEB_CANTEMP[i][j] = p->AIRTEMP;
									this->CEB_CANTEMP_OLD[i][j] = p->AIRTEMP;
									this->CEB_DUMMY[i][j] = p->DUMMY;
								}
								

//								//CALCULATION OF CANOPY ENERGY BALANCE
//								if (this->CEB_model == true)
//								{
//									SOLVE_CEB(p);
//								}
//								//CEB
//								this->CEB_CANENBAL[i][j] = p->CANENBAL;
//								this->CEB_CANTEMP[i][j] = p->CANTEMP;
//								this->CEB_CANTEMP_OLD[i][j] = p->CANTEMP;
//								this->CEB_DUMMY[i][j] = p->DUMMY;
							}

							// ************************************************
							//scn
							// update results
							// todo save original values in separate grids
							// (1) Meteorological forcing data (updated values in the canopy)
							this->SCN_Out_liquid[i][j] = p->RAINPRECIP; // save intermediate results to grid	
							this->SCN_Out_solid[i][j] = p->SNOWPRECIP;
							this->SCN_Temp[i][j] = p->AIRTEMP - 273.16;
							this->SCN_RH[i][j] = p->HUM / 100.;
							this->SCN_RAD[i][j] = p->SHORTWAVE;
							this->SCN_Wind[i][j] = p->WIND;
							this->SCN_IR[i][j] = p->IR; // intermediate results that might be relevant for energy balance snow models
							// (2) State variables
							this->SCN_Albedo[i][j] = p->ALBEDO;
							this->SCN_LOAD[i][j] = p->CANOPINTLOAD;
							this->SCN_AGE[i][j] = p->SNOWAGE;
							// (3) additional results (fluxes)
							this->SCN_SUBL[i][j] = p->SUBLIMSNOW;

//							// check mass balance
//							double watbal = p->PRECIP - p->SNOWPRECIP - p->RAINPRECIP - p->SUBLIMSNOW - (p->CANOPINTLOAD - last_ci);
//							if (watbal < 0)
//								watbal *= -1;
//							if (watbal > 1.E-9)
//								int z = 0;
//
//							// for testing
//							if (i == 50 && j == 50)
//								double TEST = p->CANTEMP - p->CANTEMP_OLD;

							// clean up
							delete p;
						}
						else
						{
							// no modified meteorological conditions when LAI and/or z0 were to small 
							this->SCN_Out_liquid[i][j] = this->QInput[i][j]; // input = output for reasons of consistency
							this->SCN_Temp[i][j] = this->TInput[i][j];
							this->SCN_RH[i][j] = this->RHInput[i][j];
							this->SCN_RAD[i][j] = this->RInput[i][j];
							this->SCN_Wind[i][j] = this->WSInput[i][j];
							this->SCN_IR[i][j] = 0.;	// $kf 2017-03-28: valid values for initialisation

							// write consistent values even if the snow canopy model isnot used
							this->SCN_Out_solid[i][j] = 0.;
							this->SCN_Albedo[i][j] = 0.2;
							this->SCN_LOAD[i][j] = 0.;
							this->SCN_AGE[i][j] = 0.;
							this->SCN_SUBL[i][j] = 0.;
							//ceb
							if (CEB_model)
							{
								this->CEB_CANENBAL[i][j] = 0;
								this->CEB_CANTEMP[i][j] = this->TInput[i][j] + 273.15;
								this->CEB_CANTEMP_OLD[i][j] = this->TInput[i][j] + 273.15;
								this->CEB_DUMMY[i][j] = 0;
							}
						}

						// redirect gridded data to time series output as well:
						this->SCN_Load_stats->wertzufuegen(i, j, zarray[(int)zn[i][j]], this->SCN_LOAD[i][j], 1 - icevalue);
						this->SCN_Age_stats->wertzufuegen(i, j, zarray[(int)zn[i][j]], this->SCN_AGE[i][j], 1 - icevalue);
						this->SCN_Alb_stats->wertzufuegen(i, j, zarray[(int)zn[i][j]], this->SCN_Albedo[i][j], 1 - icevalue);
						this->SCN_Sub_stats->wertzufuegen(i, j, zarray[(int)zn[i][j]], this->SCN_SUBL[i][j], 1 - icevalue);
						this->SCN_Out_solid_stats->wertzufuegen(i, j, zarray[(int)zn[i][j]], this->SCN_Out_solid[i][j], 1 - icevalue);
						this->SCN_Out_liq_stats->wertzufuegen(i, j, zarray[(int)zn[i][j]], this->SCN_Out_liquid[i][j], 1 - icevalue);
						this->SCN_Temp_stats->wertzufuegen(i, j, zarray[(int)zn[i][j]], this->SCN_Temp[i][j], 1 - icevalue);
						this->SCN_RH_stats->wertzufuegen(i, j, zarray[(int)zn[i][j]], this->SCN_RH[i][j], 1 - icevalue);
						this->SCN_RAD_stats->wertzufuegen(i, j, zarray[(int)zn[i][j]], this->SCN_RAD[i][j], 1 - icevalue);
						this->SCN_Wind_stats->wertzufuegen(i, j, zarray[(int)zn[i][j]], this->SCN_Wind[i][j], 1 - icevalue);
						//ceb
						if (CEB_model)
						{
							this->CEB_CANENBAL_stats->wertzufuegen(i, j, zarray[(int)zn[i][j]], this->CEB_CANENBAL[i][j], 1 - icevalue);
							this->CEB_CANTEMP_stats->wertzufuegen(i, j, zarray[(int)zn[i][j]], this->CEB_CANTEMP[i][j], 1 - icevalue);
							this->CEB_CANTEMP_OLD_stats->wertzufuegen(i, j, zarray[(int)zn[i][j]], this->CEB_CANTEMP_OLD[i][j], 1 - icevalue);
							this->CEB_DUMMY_stats->wertzufuegen(i, j, zarray[(int)zn[i][j]], this->CEB_DUMMY[i][j], 1 - icevalue);
						}
					}

					if (iNumVegetationLayers>1) // calculate interception for multiple layers
					{
						EI[0][i][j]=0.0;
						Q[0][i][j]=0.0;
						SI[0][i][j]=0.0;
						for (lay=1;lay<=iNumVegetationLayers;lay++)
						{
							if (LUSE[lay][i][j]>nodata)
							{
								double fIntercepStorage = LUtable[(int)(LUSE[lay][i][j])].getIntercepStor();
								if (fIntercepStorage<0) // if no "IntercepCap" is given in the landuse table or if the old landuse table is used, the fIntercepStorage will be set to -1 in LAtable
								{
									SImax = LAI[lay][i][j] * SIdicke; /* Beisppiel: 0.5 mm Wasserschichtdicke als SIdicke ergibt mit LAI = 10
														multipliziert den maximal moeglichen Inhalt des Speichers von 5 mm */
								}
								else
								{ // here, the interception storage content is calculated depending on the actual landuse as well as on the actual LAI
									SImax = LAI[lay][i][j]*fIntercepStorage;
								}
								
								/* z.B. only vegetation covered fraction can intercept water, the uncovered fraction will drain into lower layers or into the soil immediately */
								double input;
								if ((lay==1)||(do_ArcVeg>0)) // when arcVeg is active, all alyers are lateral fractions, so the input is the same for all layers
								{
									input = input_to_rainwater_interc;
								}
								else
								{
									input = Q[lay-1][i][j];
								}
								Q[lay][i][j]=(1.0-VBG[lay][i][j])*input;
								
								/* Bilanz des Speichers, QInput faellt in Speicher, was dort nicht unterkommt, faellt durch als Q */
								SI[lay][i][j]+=VBG[lay][i][j]*input;
								
								if (SI[lay][i][j]>SImax) 
								{
									/* outflow from interception storage */
									Q[lay][i][j]+=(SI[lay][i][j]-SImax); 
									SI[lay][i][j]=SImax;
								}
								
								/* changes in interflow storage due to evaporation 
									(concentrations of evaporating tracers remain as constant as 
									the abs. masses of non evaporating tracers) */
								double actFrac = 1.0;
								if (do_ArcVeg > 0)
								{
									actFrac = ARCVEGFRACTIONS[lay][i][j];
								}

								if ((methode == 1) || (methode == 11) || (methode == 21) || (methode == 31)) // traditional method using ETP for interception
								{
									if ((SI[lay][i][j]* actFrac) <=(ETP[lay][i][j]))
									{ /* alles kann verdunsten, Wassermenge ist limitierend, es bleibt was fuer ETP uebrig */
										EI[lay][i][j]   = SI[lay][i][j] * actFrac;
										ETP[lay][i][j] -= EI[lay][i][j];
										ETP[0][i][j]   -= EI[lay][i][j];    // index 0 holds the total amount of remaining ETP, so it must be subtracted there, too
										SI[lay][i][j]   = 0.0;
									}
									else /* SI > ETP Energiemenge ist limitierend, es bleibt also Wasser in SI drin */
									{
										EI[lay][i][j]   = ETP[lay][i][j];
										ETP[0][i][j]   -= ETP[lay][i][j];    // index 0 holds the total amount of remaining ETP, so it must be subtracted there, too
										ETP[lay][i][j]  = 0.0;
										SI[lay][i][j]  -= (EI[lay][i][j]/actFrac);
									}
								}
								else if ((methode == 2) || (methode == 12) || (methode == 22) || (methode == 32))
								{
									double reduamount=0;
									if ((SI[lay][i][j]* actFrac)<=EIP[lay][i][j])
									{ /* alles kann verdunsten, Wassermenge ist limitierend, es bleibt was fuer ETP uebrig */
										EI[lay][i][j]=SI[lay][i][j]* actFrac;
										if (ETP[lay][i][j]<EI[lay][i][j])
										{
											reduamount=ETP[lay][i][j];
										}
										else
										{
											reduamount=EI[lay][i][j];
										}
										ETP[lay][i][j]-=reduamount;
										ETP[0][i][j]-=reduamount;    // index 0 holds the total amount of remaining ETP, so it must be subtracted there, too
										SI[lay][i][j]=0.0;
									}
									else /* SI > ETP Energiemenge ist limitierend, es bleibt also Wasser in SI drin */
									{
										EI[lay][i][j]=EIP[lay][i][j];
										if (ETP[lay][i][j]<EI[lay][i][j])
										{
											reduamount=ETP[lay][i][j];
										}
										else
										{
											reduamount=EI[lay][i][j];
										}
										ETP[0][i][j]-=reduamount;    // index 0 holds the total amount of remaining ETP, so it must be subtracted there, too
										ETP[lay][i][j]-=reduamount;
										SI[lay][i][j]-=(EI[lay][i][j]/ actFrac);
									}
								}

								EI[0][i][j] += EI[lay][i][j]; // total amount of interception loss is the sum of all layers
								if (do_ArcVeg > 0)
								{
									Q[0][i][j] += Q[lay][i][j];    // for ArcVeg active, the total runoff is the fraction-weighed sum of all layers outflow (parallel!)
								}
								else
								{
									Q[0][i][j] = Q[lay][i][j];    // only the last layer gives its throughflow to the next module for stacked multilayer, this value will be overwritten for each layer, so only the last value is valid
								}
								if (do_ArcVeg>0)
								{ 
									SI[0][i][j] += SI[lay][i][j]* actFrac; // for ArcVeg, the fractions are beside each other, so each layers storage has to multiplied by the fraction to get a fraction-weighed sum for the entire cell
								}
								else
								{
									SI[0][i][j] += SI[lay][i][j]; // the total amount of intercepted water is the sum of the water stored in the seperate layers
								}
							}
							else
							{
								EI[lay][i][j]=0;
								Q[lay][i][j]=0;
								SI[lay][i][j]=0;


							}

							EIstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],EI[lay][i][j],1-icevalue);
							SIstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],SI[lay][i][j],1-icevalue);
							Qstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],Q[lay][i][j],1-icevalue);
						}	

					} // more than one vegetation layer
					else
					{ // original method with one layer only
						double fIntercepStorage = LUtable[(int)(LUSE[1][i][j])].getIntercepStor();
						if (fIntercepStorage<0) // if no "IntercepCap" is given in the landuse table or if the old landuse table is used, the fIntercepStorage will be set to -1 in LAtable
						{
							SImax = LAI[1][i][j] * SIdicke; /* Beisppiel: 0.5 mm Wasserschichtdicke als SIdicke ergibt mit LAI = 10
												multipliziert den maximal moeglichen Inhalt des Speichers von 5 mm */
						}
						else
						{ // here, the interception storage content is calculated depending on the actual landuse as well as on the actual LAI
							SImax = LAI[1][i][j]*fIntercepStorage;
						}
						
						/* z.B. only vegetation covered fraction can intercept water, the uncovered fraction will drain into lower layers or into the soil immediately */
						Q[1][i][j]=(1.0-VBG[1][i][j])*input_to_rainwater_interc; //$kf
						
						for (z=1;z<=tracer_count;z++)
						{
							if (tracer_info[z].radioactive)
							tracer_module[z].concstor[i][j]*=tracer_info[z].decfact;
							/* the uncovered fraction of soil gets the tracer immediately into the flux-grid 
				(for evaporating tracers the concentration has to be averaged by weighing it with prec. masses) */
							tracer_module[z].concflux[i][j] = tracer_info[z].evaporating?
							tracer_module[z].cinput[i][j] : /* if evaporating tracer: flux has input concentration */
							tracer_module[z].cinput[i][j]*Q[1][i][j]; /* salts: mass of salt is calculated */
							/* the other fraction: mixing or adding concentrations / masses with content of concstor */
							tracer_module[z].concstor[i][j] = tracer_info[z].evaporating ? 
							mixing(tracer_module[z].concstor[i][j],tracer_module[z].cinput[i][j],SI[1][i][j],VBG[1][i][j]*input_to_rainwater_interc) :
							/* QInput in mm/zeitintervall = Kg/m2, also ergibt sich Masse des
				nichtverdunstenden Tracers als Produkt aus QInput und Tracerconcentration */
							tracer_module[z].concstor[i][j]+tracer_module[z].cinput[i][j]*VBG[1][i][j]*input_to_rainwater_interc; 
						}
						/* Bilanz des Speichers, QInput faellt in Speicher, was dort nicht unterkommt, faellt durch als Q */
						SI[1][i][j]+=VBG[1][i][j]*input_to_rainwater_interc;
						
						if (SI[1][i][j]>SImax) 
						{
							/* tracer bilance: mixing for evaporating and summing masses for non evaporating tracers */
							for (z=1;z<=tracer_count;z++)
							{
								if (tracer_info[z].evaporating)
								{
									tracer_module[z].concflux[i][j]=
									mixing(tracer_module[z].concstor[i][j],tracer_module[z].concflux[i][j],SI[1][i][j]-SImax,Q[1][i][j]);
								}
								else
								{ /* subtracting massflux from storage and add it to the outflow */
									trms    = SI[1][i][j]>0.0 ? tracer_module[z].concstor[i][j]*((SI[1][i][j]-SImax)/SI[1][i][j]) : 0.0;
									trmsmax = tracer_info[z].max_conc*(SI[1][i][j]-SImax);
									if (trmsmax<trms) trms=trmsmax;
									tracer_module[z].concstor[i][j]-=trms;
									tracer_module[z].concflux[i][j]+=trms;
								}
							}
							/* outflow from interception storage */
							Q[1][i][j]+=(SI[1][i][j]-SImax); 
							SI[1][i][j]=SImax;
						}
						for (z=1;z<=tracer_count;z++)
						{
							if (!tracer_info[z].evaporating)
							tracer_module[z].concflux[i][j]=Q[1][i][j]>0?tracer_module[z].concflux[i][j]/Q[1][i][j]:0.0; /* make concentrations from the masses of non evaporating tracers */
						}
						
						/* changes in interception storage due to evaporation 
			(concentrations of evaporating tracers remain as constant as 
			the abs. masses of non evaporating tracers) */
						if ((methode == 1) || (methode == 11) || (methode == 21) || (methode == 31)) // traditional method using ETP for interception
						{
							if (SI[1][i][j]<=ETP[1][i][j]) 
							{ /* alles kann verdunsten, Wassermenge ist limitierend, es bleibt was fuer ETP uebrig */
								EI[1][i][j]=SI[1][i][j];
								ETP[1][i][j]-=EI[1][i][j];
								SI[1][i][j]=0.0;
							}
							else /* SI > ETP Energiemenge ist limitierend, es bleibt also Wasser in SI drin */
							{
								EI[1][i][j]=ETP[1][i][j];
								ETP[1][i][j]=0.0;
								SI[1][i][j]-=EI[1][i][j];
							}
						}
						else if ((methode == 2) || (methode == 12) || (methode == 22) || (methode == 32))
						{
							double reduamount=0;
							if (SI[1][i][j]<=EIP[1][i][j]) 
							{ /* alles kann verdunsten, Wassermenge ist limitierend, es bleibt was fuer ETP uebrig */
								EI[1][i][j]=SI[1][i][j];
								if (ETP[1][i][j]<EI[1][i][j])
								{
									reduamount=ETP[1][i][j];
								}
								else
								{
									reduamount=EI[1][i][j];
								}
								ETP[1][i][j]-=reduamount;
								SI[1][i][j]=0.0;
							}
							else /* SI > ETP Energiemenge ist limitierend, es bleibt also Wasser in SI drin */
							{
								EI[1][i][j]=EIP[1][i][j];
								if (ETP[1][i][j]<EI[1][i][j])
								{
									reduamount=ETP[1][i][j];
								}
								else
								{
									reduamount=EI[1][i][j];
								}
								ETP[1][i][j]-=reduamount;
								SI[1][i][j]-=EI[1][i][j];
							}
						}
						
						for (z=1;z<=tracer_count;z++)
						{
							if (tracer_info[z].evaporating)
							{
								if (RelCellArea!=NULL)
								{
									tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].sx+=tracer_module[z].concstor[i][j]*RelCellArea[i][j]*SI[1][i][j];
									tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].n+=SI[1][i][j];
									tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].RelCellSum+=RelCellArea[i][j]*SI[1][i][j];
								}
								else
								{
									tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].sx+=tracer_module[z].concstor[i][j]*SI[1][i][j];
									tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].n+=SI[1][i][j];
									tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].RelCellSum+=SI[1][i][j];
								}
							}
							else
							{
								if (RelCellArea!=NULL)
								{
									tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].sx+=RelCellArea[i][j]*tracer_module[z].concstor[i][j];
									tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].n++;
									tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].RelCellSum+=RelCellArea[i][j];
								}
								else
								{
									tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].sx+=tracer_module[z].concstor[i][j];
									tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].n++;
									tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].RelCellSum+=1.0;
								}
							}
						}
						// to have the same values for the totalisaation in index 0, the values have simply to be copied
						EI[0][i][j]  = EI[1][i][j];
						ETP[0][i][j] = ETP[1][i][j];
						SI[0][i][j]  = SI[1][i][j];
						Q[0][i][j]   = Q[1][i][j];


						EIstat[1]->wertzufuegen(i,j,zarray[(int)zn[i][j]],EI[0][i][j],1.0-icevalue);
						SIstat[1]->wertzufuegen(i,j,zarray[(int)zn[i][j]],SI[0][i][j],1.0-icevalue);
						Qstat[1]->wertzufuegen(i,j,zarray[(int)zn[i][j]],Q[0][i][j],1.0-icevalue);
					}

					EIstat[0]->wertzufuegen(i,j,zarray[(int)zn[i][j]],EI[0][i][j],1.0-icevalue);
					SIstat[0]->wertzufuegen(i,j,zarray[(int)zn[i][j]],SI[0][i][j],1.0-icevalue);
					Qstat[0]->wertzufuegen(i,j,zarray[(int)zn[i][j]],Q[0][i][j],1.0-icevalue);
				}
//				else
				if (icevalue>0) // add 0 for glacierized fractions of a cell in order to get weighed statistics for the entire basin
				{
					// change only pure ice cells' values
					for (lay=0;lay<=iNumVegetationLayers;lay++)
					{
						if (icevalue>=0.99999)
						{
							EI[lay][i][j] = 0;
							SI[lay][i][j] = 0;
							Q[0][i][j] = 0;
						}
					}
					Q[0][i][j] = Q[0][i][j]*(1.0-icevalue)+QInput[i][j]*icevalue;

					for (z=1;z<=tracer_count;z++)
					{
						if (tracer_info[z].evaporating)
						{
							if (RelCellArea!=NULL)
							{
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].sx+=tracer_module[z].concstor[i][j]*RelCellArea[i][j]*SI[0][i][j];
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].n+=SI[0][i][j];
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].RelCellSum+=RelCellArea[i][j]*SI[0][i][j];
							}
							else
							{
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].sx+=tracer_module[z].concstor[i][j]*SI[0][i][j];
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].n+=SI[0][i][j];
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].RelCellSum+=SI[0][i][j];
							}
						}
						else
						{
							if (RelCellArea!=NULL)
							{
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].sx+=RelCellArea[i][j]*tracer_module[z].concstor[i][j];
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].n++;
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].RelCellSum+=RelCellArea[i][j];
							}
							else
							{
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].sx+=tracer_module[z].concstor[i][j];
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].n++;
								tracer_module[z].stat->zone[zarray[(int)zn[i][j]]].RelCellSum+=1.0;
							}
						}
					}
					for (lay=0;lay<=iNumVegetationLayers;lay++)
					{

						EIstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],0,icevalue);
						SIstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],0,icevalue);
						Qstat[lay]->wertzufuegen(i,j,zarray[(int)zn[i][j]],0,icevalue);
					}
				}
			}
			else 
			{
				for (lay=0;lay<=iNumVegetationLayers;lay++)
				{	
					EI[lay][i][j] = nodata;
					Q[lay][i][j] = nodata;
				}
			}
		}
	}
	
	// $kf write results of the canopy model to grid files
	if (methode > 20)
	{
		writegridtofile(this->grd_idx_scn_temp);
		writegridtofile(this->grd_idx_scn_rh);
		writegridtofile(this->grd_idx_scn_rad);
		writegridtofile(this->grd_idx_scn_wind);
		writegridtofile(this->grd_idx_scn_age);
		writegridtofile(this->grd_idx_scn_alb);
		writegridtofile(this->grd_idx_scn_load);
		writegridtofile(this->grd_idx_scn_subl);
		writegridtofile(this->grd_idx_scn_out_s);
		writegridtofile(this->grd_idx_scn_out_l);
		//ceb
		if (CEB_model)
		{
			writegridtofile(this->grd_idx_ceb_canenbal);
			writegridtofile(this->grd_idx_ceb_cantemp);
			writegridtofile(this->grd_idx_ceb_cantemp_old);
			writegridtofile(this->grd_idx_ceb_dummy);
		}
		this->SCN_Temp_stats->calc_and_done(jahr, monat, tag, stunde, minute);
		this->SCN_RH_stats->calc_and_done(jahr, monat, tag, stunde, minute);
		this->SCN_RAD_stats->calc_and_done(jahr, monat, tag, stunde, minute);
		this->SCN_Wind_stats->calc_and_done(jahr, monat, tag, stunde,minute);
		this->SCN_Load_stats->calc_and_done(jahr, monat, tag, stunde,minute);
		this->SCN_Age_stats->calc_and_done(jahr, monat, tag, stunde,minute);
		this->SCN_Alb_stats->calc_and_done(jahr, monat, tag, stunde,minute);
		this->SCN_Sub_stats->calc_and_done(jahr, monat, tag, stunde,minute);
		this->SCN_Out_solid_stats->calc_and_done(jahr, monat, tag, stunde, minute);
		this->SCN_Out_liq_stats->calc_and_done(jahr, monat, tag, stunde, minute);
		//ceb
		if (CEB_model)
		{
			this->CEB_CANENBAL_stats->calc_and_done(jahr, monat, tag, stunde, minute);
			this->CEB_CANTEMP_stats->calc_and_done(jahr, monat, tag, stunde, minute);
			this->CEB_CANTEMP_OLD_stats->calc_and_done(jahr, monat, tag, stunde, minute);
			this->CEB_DUMMY_stats->calc_and_done(jahr, monat, tag, stunde, minute);
		}
	}
	
	for (lay=0;lay<=iNumVegetationLayers;lay++)
	{
		EIstat[lay]->calc_and_done(jahr,monat,tag,stunde, minute);
		SIstat[lay]->calc_and_done(jahr,monat,tag,stunde, minute);
		Qstat[lay]->calc_and_done(jahr,monat,tag,stunde, minute);
		writegridtofile(qindex[lay]);
		writegridtofile(eiindex[lay]);
		writegridtofile(siindex[lay]);
	}

	for (z=1;z<=tracer_count;z++)
	{
		tracer_module[z].stat->calc_and_done(jahr,monat,tag,stunde, minute);
		writegridtofile(tracer_module[z].storindex);
	}
	// for comparison with next time step:if values are identical, we are within the same hour and no new temperature grid needs to be added to the GridStack TempStack
	minute_old = minute;
	stunde_old = stunde;
	tag_old = tag;
	monat_old = monat;
	jahr_old = jahr;
}


/* 3.3.97: neu anstelle des Destructors wird der Speicher extra freigegeben.
C++ soll Standarddestructor erzeugen, da bei Borland C++ 4.5 oft beim
Freigeben (Aufruf des Destructors) unter Win32 Fehler in CW3215.dll erzeugt wurden. */
void interzeptionsmodellklasse::ClearMemory()
{
	if (do_it >0)
	{
		// $kf 2016-06-29: delete stack of temperature grids which was used
		// to calculate the daily average temperature
		// @todo: change conditional statement !!!
		if (methode > 20)
		{
			// delete temperature memory which was needed to calculate averages
			delete this->TempGrids;

			// clean up grids and stats associated to the snow canopy model
			this->SCN_Temp_stats->ClearMemory();
			this->SCN_RH_stats->ClearMemory();
			this->SCN_RAD_stats->ClearMemory();
			this->SCN_Wind_stats->ClearMemory();
			this->SCN_Load_stats->ClearMemory();
			this->SCN_Age_stats->ClearMemory();
			this->SCN_Alb_stats->ClearMemory();
			this->SCN_Sub_stats->ClearMemory();
			this->SCN_Out_solid_stats->ClearMemory();
			this->SCN_Out_liq_stats->ClearMemory();
			if (CEB_model)
			{
				this->CEB_CANENBAL_stats->ClearMemory();
				this->CEB_CANTEMP_stats->ClearMemory();
				this->CEB_CANTEMP_OLD_stats->ClearMemory();
				this->CEB_DUMMY_stats->ClearMemory();
			}
			if (this->SCN_Temp_stats) delete this->SCN_Temp_stats; this->SCN_Temp_stats = NULL;
			if (this->SCN_RH_stats) delete this->SCN_RH_stats; this->SCN_RH_stats = NULL;
			if (this->SCN_RAD_stats) delete this->SCN_RAD_stats; this->SCN_RAD_stats = NULL;
			if (this->SCN_Wind_stats) delete this->SCN_Wind_stats; this->SCN_Wind_stats = NULL;
			if (this->SCN_Load_stats) delete this->SCN_Load_stats; this->SCN_Load_stats = NULL;
			if (this->SCN_Age_stats) delete this->SCN_Age_stats; this->SCN_Age_stats = NULL;
			if (this->SCN_Alb_stats) delete this->SCN_Alb_stats; this->SCN_Alb_stats = NULL;
			if (this->SCN_Sub_stats) delete this->SCN_Sub_stats; this->SCN_Sub_stats = NULL;
			if (this->SCN_Out_solid_stats) delete this->SCN_Out_solid_stats; this->SCN_Out_solid_stats = NULL;
			if (this->SCN_Out_liq_stats) delete this->SCN_Out_liq_stats; this->SCN_Out_liq_stats = NULL;
			if (CEB_model)
			{
				if (this->CEB_CANENBAL_stats) delete this->CEB_CANENBAL_stats; this->CEB_CANENBAL_stats = NULL;
				if (this->CEB_CANTEMP_stats) delete this->CEB_CANTEMP_stats; this->CEB_CANTEMP_stats = NULL;
				if (this->CEB_CANTEMP_OLD_stats) delete this->CEB_CANTEMP_OLD_stats; this->CEB_CANTEMP_OLD_stats = NULL;
				if (this->CEB_DUMMY_stats) delete this->CEB_DUMMY_stats; this->CEB_DUMMY_stats = NULL;
			}
			if (this->SCN_Temp_zone) delete this->SCN_Temp_zone; this->SCN_Temp_zone = NULL;
			if (this->SCN_RH_zone) delete this->SCN_RH_zone; this->SCN_RH_zone = NULL;
			if (this->SCN_RAD_zone) delete this->SCN_RAD_zone; this->SCN_RAD_zone = NULL;
			if (this->SCN_Wind_zone) delete this->SCN_Wind_zone; this->SCN_Wind_zone = NULL;
			if (this->SCN_Load_zone) delete this->SCN_Load_zone; this->SCN_Load_zone = NULL;
			if (this->SCN_Age_zone) delete this->SCN_Age_zone; this->SCN_Age_zone = NULL;
			if (this->SCN_Alb_zone) delete this->SCN_Alb_zone; this->SCN_Alb_zone = NULL;
			if (this->SCN_Sub_zone) delete this->SCN_Sub_zone; this->SCN_Sub_zone = NULL;
			if (this->SCN_Out_solid_zone) delete this->SCN_Out_solid_zone; this->SCN_Out_solid_zone = NULL;
			if (this->SCN_Out_liq_zone) delete this->SCN_Out_liq_zone; this->SCN_Out_liq_zone = NULL;
			if (CEB_model)
			{
				if (this->CEB_CANENBAL_zone) delete this->CEB_CANENBAL_zone; this->CEB_CANENBAL_zone = NULL;
				if (this->CEB_CANTEMP_zone) delete this->CEB_CANTEMP_zone; this->CEB_CANTEMP_zone = NULL;
				if (this->CEB_CANTEMP_OLD_zone) delete this->CEB_CANTEMP_OLD_zone; this->CEB_CANTEMP_OLD_zone = NULL;
				if (this->CEB_DUMMY_zone) delete this->CEB_DUMMY_zone; this->CEB_DUMMY_zone = NULL;
			}
			if (this->SCN_Temp_filename) delete this->SCN_Temp_filename; this->SCN_Temp_filename = NULL;
			if (this->SCN_RH_filename) delete this->SCN_RH_filename; this->SCN_RH_filename = NULL;
			if (this->SCN_RAD_filename) delete this->SCN_RAD_filename; this->SCN_RAD_filename = NULL;
			if (this->SCN_Wind_filename) delete this->SCN_Wind_filename; this->SCN_Wind_filename = NULL;
			if (this->SCN_LOAD_filename) delete this->SCN_LOAD_filename; this->SCN_LOAD_filename = NULL;
			if (this->SCN_AGE_filename) delete this->SCN_AGE_filename; this->SCN_AGE_filename = NULL;
			if (this->SCN_ALB_filename) delete this->SCN_ALB_filename; this->SCN_ALB_filename = NULL;
			if (this->SCN_SUBL_filename) delete this->SCN_SUBL_filename; this->SCN_SUBL_filename = NULL;
			if (this->SCN_Out_solid_filename) delete this->SCN_Out_solid_filename; this->SCN_Out_solid_filename = NULL;
			if (this->SCN_Out_liq_filename) delete this->SCN_Out_liq_filename; this->SCN_Out_liq_filename = NULL;
			if (CEB_model)
			{
				if (this->CEB_CANENBAL_filename) delete this->CEB_CANENBAL_filename; this->CEB_CANENBAL_filename = NULL;
				if (this->CEB_CANTEMP_filename) delete this->CEB_CANTEMP_filename; this->CEB_CANTEMP_filename = NULL;
				if (this->CEB_CANTEMP_OLD_filename) delete this->CEB_CANTEMP_OLD_filename; this->CEB_CANTEMP_OLD_filename = NULL;
				if (this->CEB_DUMMY_filename) delete this->CEB_DUMMY_filename; this->CEB_DUMMY_filename = NULL;
			}
			free_matrix(this->SCN_Temp, mpistartrow, mpiendrow, 1, ncols, nrows); this->SCN_Temp = NULL;
			free_matrix(this->SCN_RH, mpistartrow, mpiendrow, 1, ncols, nrows); this->SCN_RH = NULL;
			free_matrix(this->SCN_RAD, mpistartrow, mpiendrow, 1, ncols, nrows); this->SCN_RAD = NULL;
			free_matrix(this->SCN_Wind, mpistartrow, mpiendrow, 1, ncols, nrows); this->SCN_Wind = NULL;
			free_matrix(this->SCN_LOAD, mpistartrow, mpiendrow, 1, ncols, nrows); this->SCN_LOAD = NULL;
			free_matrix(this->SCN_AGE, mpistartrow, mpiendrow, 1, ncols, nrows); this->SCN_AGE = NULL;
			free_matrix(this->SCN_Albedo, mpistartrow, mpiendrow, 1, ncols, nrows); this->SCN_Albedo = NULL;
			free_matrix(this->SCN_SUBL, mpistartrow, mpiendrow, 1, ncols, nrows); this->SCN_SUBL = NULL;
			free_matrix(this->SCN_Out_solid, mpistartrow, mpiendrow, 1, ncols, nrows); this->SCN_Out_solid = NULL;
			free_matrix(this->SCN_Out_liquid, mpistartrow, mpiendrow, 1, ncols, nrows); this->SCN_Out_liquid = NULL;
			free_matrix(this->SCN_IR, mpistartrow, mpiendrow, 1, ncols, nrows); this->SCN_IR = NULL;
			if (CEB_model)
			{
				free_matrix(this->CEB_CANENBAL, mpistartrow, mpiendrow, 1, ncols, nrows); this->CEB_CANENBAL = NULL;
				free_matrix(this->CEB_CANTEMP, mpistartrow, mpiendrow, 1, ncols, nrows); this->CEB_CANTEMP = NULL;
				free_matrix(this->CEB_CANTEMP_OLD, mpistartrow, mpiendrow, 1, ncols, nrows); this->CEB_CANTEMP_OLD = NULL;
				free_matrix(this->CEB_DUMMY, mpistartrow, mpiendrow, 1, ncols, nrows); this->CEB_DUMMY = NULL;
			}
			delete_grid(this->grd_idx_scn_temp);
			delete_grid(this->grd_idx_scn_rh);
			delete_grid(this->grd_idx_scn_rad);
			delete_grid(this->grd_idx_scn_wind);
			delete_grid(this->grd_idx_scn_load);
			delete_grid(this->grd_idx_scn_age);
			delete_grid(this->grd_idx_scn_alb);
			delete_grid(this->grd_idx_scn_subl);
			delete_grid(this->grd_idx_scn_out_s);
			delete_grid(this->grd_idx_scn_out_l);
			delete_grid(this->grd_idx_scn_ir);
			if (CEB_model)
			{
				delete_grid(this->grd_idx_ceb_canenbal);
				delete_grid(this->grd_idx_ceb_cantemp);
				delete_grid(this->grd_idx_ceb_cantemp_old);
				delete_grid(this->grd_idx_ceb_dummy);
			}
			if (this->SCN_Temp_gridname) free(this->SCN_Temp_gridname); this->SCN_Temp_gridname = NULL;
			if (this->SCN_RH_gridname) free(this->SCN_RH_gridname); this->SCN_RH_gridname = NULL;
			if (this->SCN_RAD_gridname) free(this->SCN_RAD_gridname); this->SCN_RAD_gridname = NULL;
			if (this->SCN_Wind_gridname) free(this->SCN_Wind_gridname); this->SCN_Wind_gridname = NULL;
			if (this->SCN_LOAD_gridname) free(this->SCN_LOAD_gridname); this->SCN_LOAD_gridname = NULL;
			if (this->SCN_AGE_gridname) free(this->SCN_AGE_gridname); this->SCN_AGE_gridname = NULL;
			if (this->SCN_ALB_gridname) free(this->SCN_ALB_gridname); this->SCN_ALB_gridname = NULL;
			if (this->SCN_SUBL_gridname) free(this->SCN_SUBL_gridname); this->SCN_SUBL_gridname = NULL;
			if (this->SCN_Out_solid_gridname) free(this->SCN_Out_solid_gridname); this->SCN_Out_solid_gridname = NULL;
			if (this->SCN_Out_liq_gridname) free(this->SCN_Out_liq_gridname); this->SCN_Out_liq_gridname = NULL;
			if (CEB_model)
			{
				if (this->CEB_CANENBAL_gridname) free(this->CEB_CANENBAL_gridname); this->CEB_CANENBAL_gridname = NULL;
				if (this->CEB_CANTEMP_gridname) free(this->CEB_CANTEMP_gridname); this->CEB_CANTEMP_gridname = NULL;
				if (this->CEB_CANTEMP_OLD_gridname) free(this->CEB_CANTEMP_OLD_gridname); this->CEB_CANTEMP_OLD_gridname = NULL;
				if (this->CEB_DUMMY_gridname) free(this->CEB_DUMMY_gridname); this->CEB_DUMMY_gridname = NULL;
			}
		}
		
		int lay;
		for (lay=0;lay<=iNumVegetationLayers;lay++)
		{
			Qstat[lay]->ClearMemory();
			EIstat[lay]->ClearMemory();
			SIstat[lay]->ClearMemory();

			if (Qstat[lay])  delete(Qstat[lay]);   Qstat[lay] = NULL;
			if (EIstat[lay]) delete(EIstat[lay]); EIstat[lay] = NULL;
			if (SIstat[lay]) delete(SIstat[lay]); SIstat[lay] = NULL;

			if (zoneEI[lay])  free(zoneEI[lay]); zoneEI[lay]=NULL;
			if (zoneSI[lay])  free(zoneSI[lay]); zoneSI[lay]=NULL;
			if (zoneQ[lay])   free(zoneQ[lay]);  zoneQ[lay]=NULL;

			if (Qstatfilename[lay]   ) free(Qstatfilename[lay]   ); Qstatfilename[lay]   = NULL;
			if (EIstatfilename[lay]  ) free(EIstatfilename[lay]  ); EIstatfilename[lay]  = NULL;
			if (SIstatfilename[lay]  ) free(SIstatfilename[lay]  ); SIstatfilename[lay]  = NULL;
		}

		for (lay=0;lay<=iNumVegetationLayers;lay++)
		{
			free_matrix(EI[lay], gridinfo[eiindex[lay]].nrl, gridinfo[eiindex[lay]].nrh, 1, ncols, nrows); EI[lay] = NULL;
			free_matrix(SI[lay], mpistartrow, mpiendrow, 1, ncols, nrows); SI[lay] = NULL;
			free_matrix(Q[lay], mpistartrow, mpiendrow, 1, ncols, nrows);  Q[lay] = NULL;

			delete_grid(eiindex[lay]);
			delete_grid(siindex[lay]);
			delete_grid(qindex[lay]);

			if (Qgridname[lay]       ) free(Qgridname[lay]       ); Qgridname[lay]       = NULL;
			if (EIgridname[lay]      ) free(EIgridname[lay]      ); EIgridname[lay]      = NULL;
			if (SIgridname[lay]      ) free(SIgridname[lay]      ); SIgridname[lay]      = NULL;
		}
                if (zoneQ) delete[] zoneQ; zoneQ=NULL;
                if (zoneEI) delete[] zoneEI; zoneEI=NULL;
                if (zoneSI) delete[] zoneSI; zoneSI=NULL;

		if (EIstat ) delete[] EIstat; EIstat = NULL;
		if (SIstat ) delete[] SIstat; SIstat = NULL;
		if (Qstat  ) delete[] Qstat;  Qstat = NULL;

		if (EIstatfilename ) delete[] EIstatfilename; EIstatfilename = NULL;
		if (SIstatfilename ) delete[] SIstatfilename; SIstatfilename = NULL;
		if (Qstatfilename  ) delete[] Qstatfilename;  Qstatfilename  = NULL;

		if (EIgridname )  delete[] EIgridname; EIgridname = NULL;
		if (SIgridname )  delete[] SIgridname; SIgridname = NULL;
		if (Qgridname  )  delete[] Qgridname;  Qgridname  = NULL;

		if (LUSE  ) delete[] LUSE; LUSE = NULL;
		if (WURZ  ) delete[] WURZ; WURZ = NULL;
		if (LAI   ) delete[] LAI;  LAI  = NULL;
		if (VBG   ) delete[] VBG;  VBG  = NULL;
		if (ETP   ) delete[] ETP;  ETP  = NULL;
		if (EIP   ) delete[] EIP;  EIP  = NULL;
		if (ETR   ) delete[] ETR;  ETR  = NULL;
		if (EI    ) delete[] EI;   EI   = NULL;
		if (SI    ) delete[] SI;   SI   = NULL;
		if (Q     ) delete[] Q;    Q    = NULL;
		if (Z0    ) delete[] Z0;   Z0   = NULL;
		if (RSC    ) delete[] RSC;   RSC   = NULL;
		if (RSI    ) delete[] RSI;   RSI   = NULL;
		if (ARCVEGFRACTIONS) delete[] ARCVEGFRACTIONS; ARCVEGFRACTIONS = NULL;

		if (siindex ) delete[] siindex; siindex = NULL;
		if (eiindex ) delete[] eiindex; eiindex = NULL;
		if (qindex  ) delete[] qindex;  qindex  = NULL;
	}
	modelobject::ClearMemory();

}

void interzeptionsmodellklasse::substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults)
{
	try
	{
		retvalue=0;
		int z;
		char *tempstr1,*tempstr2,*tempstr3;
		if (do_it==0) return;
		else
		{
			teilmodellklasse::substance_transport_ini(retvalue,p_tracer_count,p_tracer_info_table,LastResults,maxLastResults);
			if (retvalue!=0) throw retvalue;
			if (tracer_count<=0) return;
			else
			{
				tempstr1=newstring(100,sizeof(char));
				tempstr2=newstring(100,sizeof(char));
				tempstr3=newstring(100,sizeof(char));
				
				for (z=1;z<=tracer_count;z++)
				{
					sprintf(tempstr3,"concflux_tracer_%1d_%s",z,"input");
					/* if the tracer is evaporating, the input concentration is always the same as in the last storage, 
		but if it is not evaporating (salts), the input-concentration may not rise above the maximum 
		concentration, therefore another grid containing the concentration in the flux is allocated and 
		consequently searched for here in such cases */
					sprintf(tempstr2,"concflux_tracer_%1d_%s",z,_irrigation_parameters_); /* in irrigationmodellklasse the interpolated input function may be combined with tracer concentrations in irrigation water */
					sprintf(tempstr1,"concflux_tracer_%1d_%s",z,_Schneemodell_);
					tracer_module[z].cinput=gridinfo[getindex_first_of_3(tempstr1,tempstr2,tempstr3)].pointer;
				}
				if (tempstr1) free(tempstr1); tempstr1=NULL;
				if (tempstr2) free(tempstr2); tempstr2=NULL;
				if (tempstr3) free(tempstr3); tempstr3=NULL;
				return;
			}
		}
	}
	catch (int ex)
	{
		retvalue=ex;
	}
	catch (...)
	// catching all exceptions in order to allow the object to be constructed correct. This is required 
	// in order to allow a dextructor call in case of an error during construction
	{
		retvalue=103;
	}
}

void interzeptionsmodellklasse::substance_transport_delete()
{
	if (do_it==0) 
	{
		return;
	}
	else
	{
		teilmodellklasse::substance_transport_delete();
	}
}

void interzeptionsmodellklasse::flushbuffers()
{
	int z;
	if (do_it==0) return;
	else
	{
		int lay;
		for (lay=0;lay<=iNumVegetationLayers;lay++)
		{
			if (EIstat[lay]) EIstat[lay]->flushbuffer();
			if (SIstat[lay]) SIstat[lay]->flushbuffer();
			if (Qstat[lay]) Qstat[lay]->flushbuffer();
		}
		
		for (z=1;z<=tracer_count;z++)
		tracer_module[z].stat->flushbuffer();
		// $kf 2016-07-08: 
		if (methode > 20)
		{
			this->SCN_Temp_stats->flushbuffer();
			this->SCN_RH_stats->flushbuffer();
			this->SCN_RAD_stats->flushbuffer();
			this->SCN_Wind_stats->flushbuffer();
			this->SCN_Load_stats->flushbuffer();
			this->SCN_Age_stats->flushbuffer();
			this->SCN_Alb_stats->flushbuffer();
			this->SCN_Sub_stats->flushbuffer();
			this->SCN_Out_solid_stats->flushbuffer();
			this->SCN_Out_liq_stats->flushbuffer();
			//ceb
			if (CEB_model)
			{
				this->CEB_CANENBAL_stats->flushbuffer();
			}
		}
	}
}


/* ---------------- Ende der Definition der Interzeptionsmodellklasse ------------------ */

