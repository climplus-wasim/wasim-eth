/*   $Revision: 30 $ 
	$Date: 8.06.10 0:05 $  */

/*
(c) 2007 by Joerg Schulla

*/

#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <cstring>
#include <cmath>
#include <limits>
#include <signal.h>
#include "stringp.h"   /* String-funktionen */
#include "modconst.h"
#include "modef.h"     /* Definition des Ausgangsobjectes der Interpolationen und Teilmodelle, Code in modef.C */
#include "dataassim.h"
#include "wasim.h"
#ifdef _MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;


data_assimilation::data_assimilation(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults, int maxLastResults,
		double breitearg,double laengearg,double bezugsmeridianarg,char *bezeichnung,
		char *steuername,int LogLevel,int progaus,
		double p_min_to_int_end,int p_numsubb,int *p_listsubb,
		gridinfotype *p_gridinfo,
		int *p_zarray,statrecord *p_zone,firstlast *p_FirstLastCol ,
		int p_zonmax,int p_zonmaxwert,
		int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
		double p_xll,double p_yll,
		double p_cellsize,double p_nodata,
		long mina,long sta, long tga, long mta, long jra,
		long mine,long ste , long tge , long mte , long jre,
		double p_timekorr,
		LUSEtabletype *p_LUtable,
		SOILtabletype *p_SOILtable,
		MultiLayerLUSEType *p_MultiLayerLUTable,
		bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
		long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
		int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU, std::vector<SP_class*> &reservoirs, long *lakeCellsPerMPIRank,int substepcnt):
teilmodellklasse(Stream,breitearg,laengearg,bezugsmeridianarg,bezeichnung,steuername,LogLevel,progaus,p_min_to_int_end,p_numsubb,p_listsubb,
		p_gridinfo,p_zarray,p_zone,p_FirstLastCol,p_zonmax,p_zonmaxwert,
		p_ncols, p_nrows, p_mpistartrow, p_mpiendrow, p_xll, p_yll, p_cellsize, p_nodata, mine,ste, tge, mte, jre, p_timekorr,
		p_LUtable,p_SOILtable,p_MultiLayerLUTable,
		WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory, v_mpiregioncnt, p_sendHalo, p_recvHalo,
		p_doSpinUp, p_routingmode, p_vHRU,substepcnt,false),
	lakeCellsPerMPIRank(lakeCellsPerMPIRank), pondindex(0), lakeindex(0)
// initializing member variables here to 0/NULL or whatever...
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
		(*WSstream) << "file "<<__FILE__<<", compiled at "<<__DATE__<<", "<<__TIME__<<"[end-entry]";
#endif
		if (ist_vorhanden(_MPIRegions_))
		{
			MPIREG = gridinfo[getindex(_MPIRegions_)].pointer;
			MPIHLP = gridinfo[getindex(_MPIHelper_)].pointer;
		}
		else
		{
			MPIREG = NULL;
			MPIHLP = NULL;
		}

		// read control file, look for section [data_assimilation]
		std::string sectionname(_data_assimilation_);
		sectionname = std::string("[") + sectionname + "]";
		findeabschnitt(sectionname.c_str(), false);
		// overread first two lines as dummy parameters (were already read in in teilmodellklasse::teilmodellklasse)
		int dummyint;
		double dummydouble;
		fscanfln(steuerdatei, "%d",&dummyint,"",tstr);
		fscanfln(steuerdatei, "%lf",&dummydouble,"",tstr);
		// read number of entities
		fscanfln(steuerdatei, "%d", &numrules, "number of data-assimilation rules ", tstr);

		for (int n = 1; n <= numrules; n++)
		{
			// read entity string
			int depth = 1;
			tstr = fgets(tstr, _stringsize_, steuerdatei);
			std::string entityString = std::string(tstr);
			while (entityString.find("}") == std::string::npos || depth > 0)
			{
				tstr = fgets(tstr, _stringsize_, steuerdatei);
				if (std::string(tstr).find("{") != std::string::npos) depth++;
				if (std::string(tstr).find("}") != std::string::npos) depth--;
				entityString.append(std::string(tstr));
			}

			// add instance of data_assimilation_grid_entity if a grid entity was found or any other instance for other instances resp.
			if (entityString.find("grid") == 0)
			{
				vdaRules.push_back(new data_assimilation_grid_entity(this, WSstream, entityString, gridinfo, xll, yll, cellsize, ncols, nrows, mpistartrow, mpiendrow,mina,sta,tga,mta,jra, min_to_int_end,InV));
			}
			else if (entityString.find("reservoir") == 0)
			{
				if (pondindex == 0)
				{
					pondindex = getindex(_ponding_storage_);
					lakeindex = getindex(_LAKE_CODES_);
				}
				if (mpi_rank == 0)
				{
					vdaRules.push_back(new data_assimilation_reservoir_entity(this, WSstream, entityString, gridinfo, xll, yll, cellsize, ncols, nrows, mpistartrow, mpiendrow, mina,sta, tga, mta, jra, min_to_int_end, InV, reservoirs));
				}
			}
		}

		initialized=1; 

	}
	catch (int ex)
	{
		retvalue=ex;
	}
	catch (...)
	{
		retvalue=103;
	}
}



void data_assimilation::run(long me,long se,long tg,long mt,long jr)
{
	if (do_it == 0) return;
	// gibt lediglich eine Nachricht aus
	teilmodellklasse::run(me,se, tg, mt, jr,NULL,NULL);

	for (std::vector<data_assimilation_entity*>::iterator iter = vdaRules.begin(); iter != vdaRules.end(); iter++)
	{
		if ((*iter)->assimtype=="grid")
		{
			(*iter)->run(me,se, tg, mt, jr);
		}
	}


	if (vdaRules.size() > 0)
	{
//		if (mpi_rank > 0)
//		{
//			// MPIScatterGridWithMask will be called for rank 0 at the end of the ::run method in order to synchronize again with all ranks
//
//			MPIScatterGridWithMask(pondindex, lakeindex, lakeCellsPerMPIRank);
//			//		cout << "rank " << mpi_rank << " skipping routing" << endl;
//			return; // this module is currently not able to run on a cluster. Only the instance running mpi rank 0 will do channel routing in order to keep the model running at all
//		}

		if (mpi_rank == 0)
		{
			for (std::vector<data_assimilation_entity*>::iterator iter = vdaRules.begin(); iter != vdaRules.end(); iter++)
			{
				if ((*iter)->assimtype == "reservoir")
				{
					(*iter)->run(me,se, tg, mt, jr);
				}
			}
		}

	}


	if ((pondindex > 0) && (lakeindex > 0))
	{
		MPIScatterGridWithMask(pondindex, lakeindex, lakeCellsPerMPIRank); // will be called for rank >0 at the beginning of the ::run method
	}

}


/* 3.3.97: neu anstelle des Destructors wird der Speicher extra freigegeben.
C++ soll Standarddestructor erzeugen, da bei Borland C++ 4.5 oft beim
Freigeben (Aufruf des Destructors) unter Win32 Fehler in CW3215.dll erzeugt wurden. */
void data_assimilation::ClearMemory()
{
	if (do_it > 0)
	{
		for (std::vector<data_assimilation_entity*>::iterator iter = vdaRules.begin(); iter != vdaRules.end(); iter++)
		{
			delete (*iter);
		}
		vdaRules.clear();
	}
	modelobject::ClearMemory();
}

void data_assimilation::substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults)
{
	int z;
	char *tempstr1,*tempstr2,*tempstr3,*tempstr4;
	if (do_it==0) return;
	else
	{
		teilmodellklasse::substance_transport_ini(retvalue,p_tracer_count,p_tracer_info_table,LastResults,maxLastResults);
		if (retvalue!=0) throw retvalue;

		// initialize all objects that are part of substance gtransport or that need to be initialized after all other model components are initialized, like objects that deal with grids that are allocated after the construcotr of this class ran

		return;
	}
}

void data_assimilation::substance_transport_delete()
{
	if (do_it==0) return;
	else
	{
		// de-initialize all objects that were initialized in substance_tronsport_ini
		// fill in here
		teilmodellklasse::substance_transport_delete();
	}
}


void data_assimilation::flushbuffers()
{
	int z;
	if (do_it==0) return;
	else
	{
	}
}

/* ++++++++++++++++++++ end of data assimilation routines ++++++++++++++++++++++++++++++ */

/* ++++++++++++++++++++ start of implementation of the various entities liek grid and reservoir data assimilation types ++++++++++++++++++++++++++++++ */

/* ++++++++++++++++++++ abstract base class ++++++++++++++++++++++++++++++ */
data_assimilation_entity::data_assimilation_entity(data_assimilation *p_parent, wasim::mystream *Stream,
	std::string entityString, gridinfotype *gridinfo,
	double xll, double yll, double cs, int ncols, int nrows, int mpistartrow, int mpiendrow, int mina,int sta, int tga, int mta, int jra, double min_to_int_end,int InV) :
	parent(p_parent),
	WSstream(Stream),
	entityString(entityString),
	gridinfo(gridinfo),
	xll(xll), 
	yll(yll), 
	cs(cs), 
	ncols(ncols), 
	nrows(nrows),
	mpistartrow(mpistartrow), 
	mpiendrow(mpiendrow),
	mina(mina),
	sta(sta),
	tga(tga),
	mta(mta),
	jra(jra), 
	hour(0), 
	day(0), 
	month(0), 
	year(0), 
	min_to_int_end(min_to_int_end),
	InV(InV),
	oldmonth(0)
{
}


data_assimilation_entity::~data_assimilation_entity()
{
	// deleting all objects now that were allocated in this instance 
	if (infile.good())
	{
		infile.close();
	}
}

/* ++++++++++++++++++++ implementation for grid data assimilation type ++++++++++++++++++++++++++++++ */

data_assimilation_grid_entity::data_assimilation_grid_entity(data_assimilation *p_parent, wasim::mystream *Stream,
	std::string entityString, gridinfotype *gridinfo,
	double xll, double yll, double cs, int ncols, int nrows, int mpistartrow, int mpiendrow, int mina,int sta, int tga, int mta, int jra, double min_to_int_end,int InV) :
	data_assimilation_entity(p_parent,Stream,entityString,gridinfo,xll, yll, cs, ncols, nrows,mpistartrow, mpiendrow, mina,sta, tga, mta, jra, min_to_int_end,InV),
	entityName(""),defaultLiquidPerc(0),doSplitSWE(false), doWeighByGlacierizedFraction(false)
{
	try
	{
		assimtype = std::string("grid");
		// interpreting entityString --> reading parameters
		entityString = su.stripSectionBetweenBrackets(entityString);
		su.stripString(entityString, "grid-ID", entityName);

		char charEntityName[_stringsize_];
		sprintf(charEntityName, "%s", entityName.c_str());
		if (parent->ist_vorhanden(charEntityName))
		{
			gridindex = parent->getindex(charEntityName);
			grid = gridinfo[gridindex].pointer;
		}
		else
		{
			(*WSstream) << "FATAL: \n---------------------- error in data assimilation model ------------------- \n"
				"entity "<< entityName << " not available as grid in the actual model configuration\n"
				"-----------------------------------------------------------------------\n" << "[end-entry]";
			throw (1210);
		}
		su.stripString(entityString, "listfilename", infilename);

		FirstLastCol = parent->FirstLastCol;
		zn = parent->zn;
		this->MPIREG = parent->modelobject::MPIREG;
		mpi_rank = parent->mpi_rank;
		zarray = parent->zarray;
		nodata = parent->nodata;

		int pos = entityString.find("options");
		if (pos != std::string::npos)
		{
			entityString = entityString.substr(pos);
			entityString = su.stripSectionBetweenBrackets(entityString);
			
			std::string sSplit = "false";
			su.stripString(entityString, "split-into-liquid-and-solid", sSplit);
			if (sSplit == std::string("true"))
			{
				doSplitSWE = true;
			}

			std::string sWeighGlac="false";
			su.stripString(entityString, "weighing-by-glacierized-fraction", sWeighGlac);
			if (sWeighGlac == std::string("true"))
			{
				doWeighByGlacierizedFraction = true;
			}
			
			if (doSplitSWE)
			{
				su.stripDouble(entityString, "default-liquid-percentage", &defaultLiquidPerc);
				if (parent->ist_vorhanden(_snowstorage_liquid_))
				{
					gridindexliquid = parent->getindex(_snowstorage_liquid_);
					SWELIQUID = gridinfo[gridindexliquid].pointer;
				}
				else
				{
					(*WSstream) << "FATAL: \n---------------------- error in data assimilation model ------------------- \n"
						"entity " << _snowstorage_liquid_ << " not available as grid in the actual model configuration\n"
						"-----------------------------------------------------------------------\n" << "[end-entry]";
					throw (1211);
				}
				if (parent->ist_vorhanden(_snowstorage_solid_))
				{
					gridindexsolid = parent->getindex(_snowstorage_solid_);
					SWESOLID = gridinfo[gridindexsolid].pointer;
				}
				else
				{
					(*WSstream) << "FATAL: \n---------------------- error in data assimilation model ------------------- \n"
						"entity " << _snowstorage_solid_ << " not available as grid in the actual model configuration\n"
						"-----------------------------------------------------------------------\n" << "[end-entry]";
					throw (1212);
				}
			}

			if (doWeighByGlacierizedFraction)
			{
				if (parent->ist_vorhanden(_TotalSnowStorage_onGlacier_))
				{
					gridindexTotalSWEonglac = parent->getindex(_TotalSnowStorage_onGlacier_);
					SWETOTGL = gridinfo[gridindexTotalSWEonglac].pointer;
				}
				else
				{
					(*WSstream) << "FATAL: \n---------------------- error in data assimilation model ------------------- \n"
						"entity " << _TotalSnowStorage_onGlacier_ << " not available as grid in the actual model configuration\n"
						"set \"weighing-by-glacierized-fraction\" to false when glacier model is inactive \n"
						"-----------------------------------------------------------------------\n" << "[end-entry]";
					throw (1214);
				}

				if (parent->ist_vorhanden(_SolidSnowStorage_onGlacier_))
				{
					gridindexsolidSWEOnglac = parent->getindex(_SolidSnowStorage_onGlacier_);
					SWESOLGL = gridinfo[gridindexsolidSWEOnglac].pointer;
				}
				else
				{
					(*WSstream) << "FATAL: \n---------------------- error in data assimilation model ------------------- \n"
						"entity " << _SolidSnowStorage_onGlacier_ << " not available as grid in the actual model configuration\n"
						"set \"weighing-by-glacierized-fraction\" to false when glacier model is inactive \n"
						"-----------------------------------------------------------------------\n" << "[end-entry]";
					throw (1215);
				}

				if (parent->ist_vorhanden(_LiquidSnowStorage_onGlacier_))
				{
					gridindexliquidSWEOnglac = parent->getindex(_LiquidSnowStorage_onGlacier_);
					SWELIQGL = gridinfo[gridindexliquidSWEOnglac].pointer;
				}
				else
				{
					(*WSstream) << "FATAL: \n---------------------- error in data assimilation model ------------------- \n"
						"entity " << _LiquidSnowStorage_onGlacier_ << " not available as grid in the actual model configuration\n"
						"set \"weighing-by-glacierized-fraction\" to false when glacier model is inactive \n"
						"-----------------------------------------------------------------------\n" << "[end-entry]";
					throw (1216);
				}

				if (parent->ist_vorhanden(_GlacierizedCells_))
				{
					gridindexglacfrac = parent->getindex(_GlacierizedCells_);
					GLACFRAC = gridinfo[gridindexglacfrac].pointer;
				}
				else
				{
					(*WSstream) << "FATAL: \n---------------------- error in data assimilation model ------------------- \n"
						"entity " << _TotalSnowStorage_onGlacier_ << " not available as grid in the actual model configuration\n"
						"set \"weighing-by-glacierized-fraction\" to false when glacier model is inactive \n"
						"-----------------------------------------------------------------------\n" << "[end-entry]";
					throw (1214);
				}
			}
		}

		// open list file, reading x, y and elevation and seeking the first record with a file name
		infile.open(infilename.c_str());
		if (!infile.good())
		{
			(*WSstream) << "FATAL: \n---------------------- error in data assimilation model ------------------- \n"
				"error opening input list file " << infilename << " for grid-type data assimilation in entity " << entityName << "\n"
				"-----------------------------------------------------------------------\n" << "[end-entry]";
			throw (1200);

		}
		// winding forward until the actual or a later date is reached, whatever comes first
		tm *start_tm = new(tm);
		start_tm->tm_hour = sta - (int)((InV-min_to_int_end)/60);
		start_tm->tm_mday = tga;
		start_tm->tm_mon = mta - 1;
		start_tm->tm_year = jra - 1900;
		start_tm->tm_min = 0;
		start_tm->tm_isdst = 0;
		start_tm->tm_sec = 0;
#ifdef WIN64
		startdatetime = _mkgmtime(start_tm);
#else
		startdatetime = timegm(start_tm);
#endif
		//startdatetime = mktime(&start_tm);

		char *tstr = new char[2048];
		infile.getline(tstr, 2048); // reading first line with the comment
		std::string syear, smonth, sday, shour;
		infile >> syear >> smonth >> sday >> shour >> xcoord;
		infile >> syear >> smonth >> sday >> shour >> ycoord;
		infile >> syear >> smonth >> sday >> shour >> elevation;
		infile >> syear >> smonth >> sday >> shour >> entityNameInFile;
		// reading the first entry, i.e. the filename of the grid to be read when the given time step is reached
		infile >> year >> month >> day >> hour >> listFileEntry;
		
		tm *next_tm = new(tm);
		next_tm->tm_hour = hour - (int)((InV - min_to_int_end)/60);
		next_tm->tm_mday = day;
		next_tm->tm_mon = month - 1;
		next_tm->tm_year = year - 1900;
		next_tm->tm_min = 0;
		next_tm->tm_isdst = 0;
		next_tm->tm_sec = 0;
#ifdef WIN64
		nextdatetime = _mkgmtime(next_tm);
#else
		nextdatetime = timegm(next_tm);
#endif
		//nextdatetime = mktime(&next_tm);

		while ((nextdatetime < startdatetime) && (infile.good()))
		{
			if (infile.good())
			{
				infile >> year >> month >> day >> hour >> listFileEntry;
			}
			next_tm->tm_hour = hour - (int)((InV - min_to_int_end)/60);
			next_tm->tm_mday = day;
			next_tm->tm_mon = month - 1;
			next_tm->tm_year = year - 1900;
			next_tm->tm_min = 0;
			next_tm->tm_isdst = 0;
			next_tm->tm_sec = 0;
#ifdef WIN64
			nextdatetime = _mkgmtime(next_tm);
#else
			nextdatetime = timegm(next_tm);
#endif
			//nextdatetime = mktime(&next_tm);
		}
		delete start_tm;
		delete next_tm;
	}
	catch (int ex)
	{
		throw(ex);
	}
	catch (...)
	{
		throw(1299);
	}
}

data_assimilation_grid_entity::~data_assimilation_grid_entity()
{
	// deleting all objects now that were allocated in this instance 
}

void data_assimilation_grid_entity::run(long me, long se, long tg, long mt, long jr)
{
	tm *act_tm = new(tm);
	act_tm->tm_hour = se - (int)((InV - min_to_int_end) / 60);
	act_tm->tm_mday = tg;
	act_tm->tm_mon = mt - 1;
	act_tm->tm_year = jr - 1900;
	act_tm->tm_min = 0;
	act_tm->tm_isdst = 0;
	act_tm->tm_sec = 0;
#ifdef WIN64
	actdatetime = _mkgmtime(act_tm);
#else
	actdatetime = timegm(act_tm);
#endif
	//actdatetime = mktime(&act_tm);
	if (actdatetime == nextdatetime)
	{
		// reading the external grid and replace the existing grid with the read ones values (not the grid itself because of potentially invalid pointers)
		char *filename = new char[1024];
		sprintf(filename, "%s", listFileEntry.c_str());
		parent->readgridfromfile(filename, // filename
			&(gridinfo[gridindex].pointer), // pointer to the grid
			gridinfo[gridindex].fillcode, // fill missing values (0 or 1)
			false, // bAllocateGridMemory
			true, // bShowMessages
			gridinfo[gridindex].defaultValue, // default value
			false, // limitAutoFill
			true, // MPIFullGrid
			false, // MPIRank=Only
			false); // MPISyncHalo

// debugging only	parent->printgridfile("temp.grd", gridinfo[gridindex].pointer, false, true, false);

		// if SWE should be applied to liquid and solid storage split it now according to tha existing fractions or to the default fractioning
		// regard MPI ranks here, loop from mpistartrow to mpiendrow
		// if glacier model is active, the respective values must be set for the glacier grids of SWE_tot, SWE_liq and SWE_solid as well
		if (doSplitSWE)
		{
			for (int i = mpistartrow; i <= mpiendrow; i++)
			{
				for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
				{
					if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
					{
						if (SWESOLID[i][j] > 0)
						{ // simplest case: cell has already snow on it, so the distribution can be taken form the existing fractions
							double liquidfraction = SWELIQUID[i][j] / (SWELIQUID[i][j] + SWESOLID[i][j]);
							SWELIQUID[i][j] = grid[i][j] * liquidfraction;
							SWESOLID[i][j] = grid[i][j] * (1.0 - liquidfraction);
						}
						else
						{   // if no SWE is already there, we have two possibililtites: a) use a fix fraction, or b) look at the neighboring cells and take the nearest value from there, 
							// which may be more complex than expected, since the neighboring cell may be a cell in another elevation range
							// version a) first...
							SWELIQUID[i][j] = grid[i][j] * double(defaultLiquidPerc/100.0);
							SWESOLID[i][j] = grid[i][j] * double((100.0 - defaultLiquidPerc)/100.0);
						}
						if (doWeighByGlacierizedFraction)
						{
							if (SWESOLGL[i][j] > 0)
							{ // simplest case: cell has already snow on it, so the distribution can be taken form the existing fractions
								double liquidfraction = SWELIQGL[i][j] / (SWELIQGL[i][j] + SWESOLGL[i][j]);
								SWELIQGL[i][j] = grid[i][j] * liquidfraction;
								SWESOLGL[i][j] = grid[i][j] * (1.0 - liquidfraction);
							}
							else
							{   // if no SWE is already there, we have two possibililtites: a) use a fix fraction, or b) look at the neighboring cells and take the nearest value from there, 
								// which may be more complex than expected, since the neighboring cell may be a cell in another elevation range
								// version a) first...
								SWELIQGL[i][j] = grid[i][j] * double(defaultLiquidPerc / 100.0);
								SWESOLGL[i][j] = grid[i][j] * double((100.0 - defaultLiquidPerc) / 100.0);
							}
							SWETOTGL[i][j] = grid[i][j];

							double factorOffGlac = ((GLACFRAC[i][j] >= 1.0) ? 0.0 : 1.0); // if the cell is 100% glacierized, only the GL-grids need to hav valid vaklues to be compatible with the running model
							double factorOnGlac = ((GLACFRAC[i][j] <= 0.0) ? 0.0 : 1.0); // -9999 or 0 means: no glacier at all on tis cell, so the SWE values for total, solid and liquid SWE can be set to 0
							SWELIQUID[i][j] *= factorOffGlac;
							SWESOLID[i][j] *= factorOffGlac;
							SWETOTGL[i][j] *= factorOnGlac;
							SWELIQGL[i][j] *= factorOnGlac;
							SWESOLGL[i][j] *= factorOnGlac;
						}
					}
				}
			}
		}


		// reading next record from infile
//		signal(SIGINT, SIG_IGN);
		if (infile.good())
		{
			infile >> year >> month >> day >> hour >> listFileEntry;
		}
//		signal(SIGINT, SIG_DFL);
		tm *next_tm = new(tm);
		next_tm->tm_hour = hour - (int)((InV - min_to_int_end) / 60);
		next_tm->tm_mday = day;
		next_tm->tm_mon = month - 1;
		next_tm->tm_year = year - 1900;
		next_tm->tm_min = 0;
		next_tm->tm_isdst = 0;
		next_tm->tm_sec = 0;
#ifdef WIN64
		nextdatetime = _mkgmtime(next_tm);
#else
		nextdatetime = timegm(next_tm);
#endif
		//nextdatetime = mktime(&next_tm);
		delete next_tm;
	}
	delete act_tm;
}

/* ++++++++++++++++++++ implementation for reservoir data assimilation type ++++++++++++++++++++++++++++++ */

data_assimilation_reservoir_entity::data_assimilation_reservoir_entity(data_assimilation * p_parent, wasim::mystream * Stream,
	std::string entityString, gridinfotype * gridinfo, 
	double xll, double yll, double cs, int ncols, int nrows, int mpistartrow, int mpiendrow, int mina, int sta, int tga, int mta, int jra, double min_to_int_end, int InV, std::vector<SP_class*> &reservoirs) :
	data_assimilation_entity(p_parent, Stream, entityString, gridinfo, xll, yll, cs, ncols, nrows, mpistartrow, mpiendrow, mina, sta, tga, mta, jra, min_to_int_end, InV),
	lakeID(0),
	SP(NULL),
	AverageWaterTable(0),
	tstr(NULL),
	skiprows(0)
{
	assimtype = std::string("reservoir");

	try
	{
		// interpreting entityString --> reading parameters
		entityString = su.stripSectionBetweenBrackets(entityString);
		su.stripInt(entityString, "lake-ID", &lakeID);

		for (std::vector<SP_class*>::iterator iter = reservoirs.begin(); iter != reservoirs.end(); iter++)
		{
			if ((*iter)->SP_code == lakeID)
			{
				SP = (*iter);
				break;
			}

		}
		if (!SP)
		{
			(*WSstream) << "FATAL: \n---------------------- error in data assimilation model ------------------- \n"
				"reservoir " << lakeID << " as defined for data assimilation not found in the actual model configuration\n"
				"-----------------------------------------------------------------------\n" << "[end-entry]";
			throw (1220);
		}
		su.stripDouble(entityString, "AverageWaterTable", &AverageWaterTable);
		su.stripString(entityString, "listfilename", infilename);

		FirstLastCol = parent->FirstLastCol;
		zn = parent->zn;
		this->MPIREG = parent->modelobject::MPIREG;
		mpi_rank = parent->mpi_rank;
		zarray = parent->zarray;
		nodata = parent->nodata;

		int pos = entityString.find("options");
		if (pos != std::string::npos)
		{
			entityString = entityString.substr(pos);
			entityString = su.stripSectionBetweenBrackets(entityString);

			sFormat = "wasim";
			su.stripString(entityString, "format", sFormat);
			if ((sFormat.find("wasim")==std::string::npos) && (sFormat.find("bafu") == std::string::npos))
			{
				(*WSstream) << "FATAL: \n---------------------- error in data assimilation model ------------------- \n"
					"for reservoir " << lakeID << ": format "<< sFormat <<"not defined \n"
					"-----------------------------------------------------------------------\n" << "[end-entry]";
				throw (1220);
			}
			su.stripInt(entityString, "skiprows", &skiprows);

		}

		// open list file, reading x, y and elevation and seeking the first record with a file name
		infile.open(infilename.c_str());
		if (!infile.good())
		{
			(*WSstream) << "FATAL: \n---------------------- error in data assimilation model ------------------- \n"
				"error opening input list file " << infilename << " for reservoir-type data assimilation for reservoir " << lakeID << "\n"
				"-----------------------------------------------------------------------\n" << "[end-entry]";
			throw (1200);

		}
		// winding forward until the actual or a later date is reached, whatever comes first
		tm *start_tm = new(tm);
		start_tm->tm_hour = sta - (int)((InV - min_to_int_end) / 60);
		start_tm->tm_mday = tga;
		start_tm->tm_mon = mta - 1;
		start_tm->tm_year = jra - 1900;
		start_tm->tm_min = 0;
		start_tm->tm_isdst = 0;
		start_tm->tm_sec = 0;
#ifdef WIN64
		startdatetime = _mkgmtime(start_tm);
#else
		startdatetime = timegm(start_tm);
#endif
		//startdatetime = mktime(&start_tm);

		tstr = new char[2048];
		if (sFormat.find("bafu")!=std::string::npos)
		{

			infile.getline(tstr, 2048); // reading first line with the comment
			entityNameInFile = std::string(tstr);
			infile.getline(tstr, 2048); // reading second line with the comment
			infile.getline(tstr, 2048); // reading third line with the comment
			std::string syear, smonth, sday, shour;

			infile.getline(tstr, 2048);
			sLine = std::string(tstr);
			year = atoi(sLine.substr(22, 4).c_str());
			month = atoi(sLine.substr(27, 2).c_str());
			day = atoi(sLine.substr(30, 2).c_str());
			hour = atoi(sLine.substr(33, 2).c_str());

			tm *next_tm = new(tm);
			next_tm->tm_hour = hour - (int)((InV - min_to_int_end) / 60);
			next_tm->tm_mday = day;
			next_tm->tm_mon = month - 1;
			next_tm->tm_year = year - 1900;
			next_tm->tm_min = 0;
			next_tm->tm_isdst = 0;
			next_tm->tm_sec = 0;
#ifdef WIN64
			nextdatetime = _mkgmtime(next_tm);
#else
			nextdatetime = timegm(next_tm);
#endif
			//nextdatetime = mktime(&next_tm);

			if (nextdatetime == -1.0)
			{
				(*WSstream) << "FATAL: \n---------------------- error in data assimilation model ------------------- \n"
					"error reading datetime format in " << infilename << " for reservoir-type data assimilation for reservoir " << lakeID << "\n"
					"-----------------------------------------------------------------------\n" << "[end-entry]";
				throw (1205);

			}

			while ((nextdatetime < startdatetime)&&(infile.good()))
			{
				if (infile.good())
				{
					infile.getline(tstr, 2048);
					sLine = std::string(tstr);
					year = atoi(sLine.substr(22, 4).c_str());
					month = atoi(sLine.substr(27, 2).c_str());
					day = atoi(sLine.substr(30, 2).c_str());
					hour = atoi(sLine.substr(33, 2).c_str());
				}
				else
				{
					(*WSstream) << "WARNING: \n---------------------- possible error in data assimilation model ------------------- \n"
						"no valid date/time record fond in " << infilename << " for reservoir-type data assimilation for reservoir " << lakeID << "\n"
						"-----------------------------------------------------------------------\n" << "[end-entry]";
				}
				next_tm->tm_hour = hour - (int)((InV - min_to_int_end) / 60);
				next_tm->tm_mday = day;
				next_tm->tm_mon = month - 1;
				next_tm->tm_year = year - 1900;
				next_tm->tm_min = 0;
				next_tm->tm_isdst = 0;
				next_tm->tm_sec = 0;
#ifdef WIN64
				nextdatetime = _mkgmtime(next_tm);
#else
				nextdatetime = timegm(next_tm);
#endif
				//nextdatetime = mktime(&next_tm);
				if (oldmonth != month)
				{
					(*WSstream) << ".";
					oldmonth = month;
				}
			}
			if (infile.good()) waterlevel = atof(substr(tstr, 39, 9));
			(*WSstream) << "[end-entry]";
			delete next_tm;
		}
		else
		{ // format = "wasim" or omitted format specifier (i.e. wasim is default)

			infile.getline(tstr, 2048); // reading first line with the comment
			std::string syear, smonth, sday, shour;
			infile >> syear >> smonth >> sday >> shour >> xcoord;
			infile >> syear >> smonth >> sday >> shour >> ycoord;
			infile >> syear >> smonth >> sday >> shour >> elevation;
			infile >> syear >> smonth >> sday >> shour >> entityNameInFile;
			// reading the first entry, i.e. the waterlevel to be read when the given time step is reached
			infile >> year >> month >> day >> hour >> waterlevel;

			tm *next_tm = new(tm);
			next_tm->tm_hour = hour - (int)((InV - min_to_int_end) / 60);
			next_tm->tm_mday = day;
			next_tm->tm_mon = month - 1;
			next_tm->tm_year = year - 1900;
			next_tm->tm_min = 0;
			next_tm->tm_isdst = 0;
			next_tm->tm_sec = 0;
#ifdef WIN64
			nextdatetime = _mkgmtime(next_tm);
#else
			nextdatetime = timegm(next_tm);
#endif
			//nextdatetime = mktime(&next_tm);

			if (nextdatetime == -1.0)
			{
				(*WSstream) << "FATAL: \n---------------------- error in data assimilation model ------------------- \n"
					"error reading datetime format in " << infilename << " for reservoir-type data assimilation for reservoir " << lakeID << "\n"
					"-----------------------------------------------------------------------\n" << "[end-entry]";
				throw (1205);

			}

			while ((nextdatetime < startdatetime) && (infile.good()))
			{
				if (infile.good())
				{
					infile >> year >> month >> day >> hour >> waterlevel;
				}
				else
				{
					(*WSstream) << "WARNING: \n---------------------- possible error in data assimilation model ------------------- \n"
						"no valid date/time record fond in " << infilename << " for reservoir-type data assimilation for reservoir " << lakeID << "\n"
						"-----------------------------------------------------------------------\n" << "[end-entry]";
				}
				next_tm->tm_hour = hour - (int)((InV - min_to_int_end) / 60);
				next_tm->tm_mday = day;
				next_tm->tm_mon = month - 1;
				next_tm->tm_year = year - 1900;
				next_tm->tm_min = 0;
				next_tm->tm_isdst = 0;
				next_tm->tm_sec = 0;
#ifdef WIN64
				nextdatetime = _mkgmtime(next_tm);
#else
				nextdatetime = timegm(next_tm);
#endif
				//nextdatetime = mktime(&next_tm);
				if (oldmonth != month)
				{
					(*WSstream) << ".";
					oldmonth = month;
				}
			}
			(*WSstream) << "[end-entry]";
			delete next_tm;
		}
		delete start_tm;
	}
	catch (int ex)
	{
		throw(ex);
	}
	catch (...)
	{
		throw(1299);
	}

}

data_assimilation_reservoir_entity::~data_assimilation_reservoir_entity()
{
	delete[] tstr;
}

void data_assimilation_reservoir_entity::run(long me,long se, long tg, long mt, long jr)
{
	tm *act_tm = new(tm);
	act_tm->tm_hour = se - (int)((InV - min_to_int_end) / 60);
	act_tm->tm_mday = tg;
	act_tm->tm_mon = mt - 1;
	act_tm->tm_year = jr - 1900;
	act_tm->tm_min = 0;
	act_tm->tm_isdst = 0;
	act_tm->tm_sec = 0;
#ifdef WIN64
	actdatetime = _mkgmtime(act_tm);
#else
	actdatetime = timegm(act_tm);
#endif
	//actdatetime = mktime(&act_tm);
	if (actdatetime == nextdatetime)
	{

#ifdef _DEBUG
		if (waterlevel <= nodata)
		{
			(*WSstream) << "DEBUG: skipping data assimilation for reservoir " << lakeID << ": external watertable is nodata (" << waterlevel << ")" << "[end-entry]";
		}
		else
		{
			double wt = SP->GetWaterLevelByVolume(SP->Vakt[0]);
			(*WSstream) << "DEBUG: old value for watertable in reservoir " << lakeID << ": " << wt << " m" << "; Vakt=" << SP->Vakt[0] << " m3" << "[end-entry]";
		}
#endif
		// set the volume of the lake defined in lakeID to the value corresponding to waterlevel by calling the respective function
		if (waterlevel > nodata)
		{
			double V = SP->GetVolumeByWaterLevel(waterlevel);
			SP->SetPonds(V);
			SP->Vakt[0] = V;
			SP->actWaterLevel = waterlevel;
		}
#ifdef _DEBUG
		if (waterlevel > nodata)
		{
			(*WSstream) << "DEBUG: new value for watertable in reservoir " << lakeID << ": " << SP->actWaterLevel << " m" << "; Vakt=" << SP->Vakt[0] << " m3" << "[end-entry]";
		}
#endif

		// if configured, skip as many rows as defined in parameter skiprows
		for (int i=0;i<skiprows;i++) infile.getline(tstr, 2048);
		// reading next record from infile
		if (sFormat.find("bafu") != std::string::npos)
		{
			infile.getline(tstr, 2048);
			sLine = std::string(tstr);
			year = atoi(sLine.substr(22, 4).c_str());
			month = atoi(sLine.substr(27, 2).c_str());
			day = atoi(sLine.substr(30, 2).c_str());
			hour = atoi(sLine.substr(33, 2).c_str());
		}
		else
		{
			if (infile.good())
			{
				infile >> year >> month >> day >> hour >> waterlevel;
			}
		}
		tm *next_tm = new(tm);
		next_tm->tm_hour = hour - (int)((InV - min_to_int_end) / 60);
		next_tm->tm_mday = day;
		next_tm->tm_mon = month - 1;
		next_tm->tm_year = year - 1900;
		next_tm->tm_min = 0;
		next_tm->tm_isdst = 0;
		next_tm->tm_sec = 0;
#ifdef WIN64
		nextdatetime = _mkgmtime(next_tm);
#else
		nextdatetime = timegm(next_tm);
#endif
		//nextdatetime = mktime(&next_tm);
		delete next_tm;
	}
	delete act_tm;
}
