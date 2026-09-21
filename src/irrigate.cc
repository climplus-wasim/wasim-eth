/*   $Revision: 30 $ 
	$Date: 16.05.10 0:09 $  */

/*
(c) 1995-2007 by Joerg Schulla
more version history: see source safe comments

Teilmodell : irrigate
Zweck interpolieren der Bewaesserung und/oder optimieren der Wassergaben
in Abhaengigkeit von der Bodenfeuchte ium Wurzelbereich

date: 16.6.98, version 1
date       : march, 9 1999
version    : 1a: hour-entry can be given for start or for end of time step via min_to_int_end

date       : december, 11 1999
version    : stream output into all possible childs of baseclass 
				wasim::mystream now possible. A stream has to be given 
				to the Application object at initialization (constructor 
				parameter).
				Also new: Exception handling is used for exit from the program.
				All old calls of "exit(n);" are replaced by a "throw n;" in order
				to allow the modules to be used in a DLL under the WIN32-environment.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "stringp.h"
#include "modconst.h"
#include "modef.h"
#include "teilmod.h"
#include "wasim.h"

using namespace std;


irrigationclass::irrigationclass(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,int maxLastResults,
	double breitearg,double laengearg,double bezugsmeridianarg,const char *bezeichnung,
	char *steuername,int LogLevel,
	int progaus,double p_min_to_int_end,int p_numsubb,int *p_listsubb,
	gridinfotype *p_gridinfo,int *p_zarray,statrecord *p_zone,firstlast *p_FirstLastCol ,
	int p_zonmax, int p_zonmaxwert, int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow, double p_xll, double p_yll,
	double p_cellsize,double p_nodata,
	long mine,long ste , long tge , long mte , long jre, double p_timekorr,
	LUSEtabletype *p_LUtable,SOILtabletype *p_SOILtable,MultiLayerLUSEType *p_MultiLayerLUTable,
	irrigationtabletype *p_irrigationtable,
	bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
	long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
	int p_doSpinUp,bool p_routingmode,std::vector<HRUMapElement*> *p_vHRU):
teilmodellklasse(Stream,breitearg,laengearg,bezugsmeridianarg,bezeichnung,steuername,LogLevel,progaus,p_min_to_int_end,p_numsubb,p_listsubb,
	p_gridinfo,p_zarray,p_zone,p_FirstLastCol,p_zonmax,p_zonmaxwert,
	p_ncols, p_nrows, p_mpistartrow,p_mpiendrow, p_xll, p_yll, p_cellsize, p_nodata, mine,ste, tge, mte, jre, p_timekorr, p_LUtable, p_SOILtable, p_MultiLayerLUTable,
	WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory, v_mpiregioncnt, p_sendHalo, p_recvHalo,
	p_doSpinUp, p_routingmode,p_vHRU),
P(0),SOIL(0),WURZ(0),IRRIG(0),SB(0),SBLOW(0),SBHIGH(0),ETR2ETP(0),MAXPUMP(0),IRRIGCODE(0),ICE(0),PUMPCONC(0),
p(0),luse(0),soil(0),wurz(0),irrig(0),sb(0),maxpump(0),Qges(0),maxirrig(0),maxirrigFromReservoirs(0),IRRIGgwstat(0),IRRIGswstat(0),
riverflow_table(0),IRRIGgwstatcode(0),IRRIGswstatcode(0),IRRIGgwstatfilename(0),
IRRIGswstatfilename(0),TG(0),methode(0),rtarget(0),abfrel(0),code(0),zoneGW(0),zoneSW(0)
{
	try
	{
		retvalue=0;

		int i;
		if (doSpinUp == 1) return; // will not be processed in SpinUp mode

		if (routingmode)
		{
			do_it = 0; // wether the submodels code for do_it is set to 0 or 1: when routingmode is active, this submodel should never run
		}

		if (do_it == 0) return;
#ifdef __TIME__
		(*WSstream) << "file "<<__FILE__<<", compiled at "<<__DATE__<<", "<<__TIME__<<"[end-entry]";
		//  (*WSstream) << "file "__FILE__", last modified at "__TIMESTAMP__", compiled at "__DATE__", "__TIME__"\n");
#endif

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
		}
		if (MPIHLP == NULL)
		{
			if (ist_vorhanden(_MPIHelper_))
			{
				MPIHLP = gridinfo[getindex(_MPIHelper_)].pointer;
			}
			else
			{
				MPIHLP = NULL;
			}
		}

		irrigationtable = p_irrigationtable;

		if (MultiLayerLUTable)
		{
			iNumVegetationLayers=MultiLayerLUTable[1].getNumLayers();
		}
		else
		{
			iNumVegetationLayers=1;
		}


		IRRIGgwstatfilename = newstring(_stringsize_+1,1);
		IRRIGswstatfilename = newstring(_stringsize_+1,1);

		fscanfln(steuerdatei,"%s%d",IRRIGgwstatfilename,&IRRIGgwstatcode  ,"Name des GWSTANDstatfiles fuer "   ,tstr); /* groundwater level */
		fscanfln(steuerdatei,"%s%d",IRRIGswstatfilename,&IRRIGswstatcode  ,"Name des GWSTANDstatfiles fuer "   ,tstr); /* groundwater level */

		if ((fclose(steuerdatei))==EOF)
		{
			(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"[end-entry]";
			throw(11);
		}

		P         = gridinfo[getindex(_Niederschlag_)].pointer;
		IRRIGCODE = gridinfo[getindex(_irrigation_codes_)].pointer;
		MAXPUMP   = gridinfo[getindex(_Max_Pumping_)].pointer;
		IRRIG     = gridinfo[getindex(_irrigation_amounts_)].pointer;
		//  SOIL      = gridinfo[getindex(_Bodenarten_)].pointer;  /* soil types (texture) */
		SB        = gridinfo[getindex(_SB_1_grid_)].pointer;   /* mean theta within root zone */
		SBLOW     = gridinfo[getindex(_SBLOW_grid_)].pointer;  /* in unsatzonmodell during initialization the minimum soil moisture for irrigated plants id calculated from minimum hydraulic head */
		SBHIGH    = gridinfo[getindex(_SBHIGH_grid_)].pointer; /* initialized in unsatzonmodell during initialization -> soil moisture to which soil is irrigated if SBLOW was reached */
		ETR2ETP   = gridinfo[getindex(_ETR2ETP_)].pointer; // initialized in unsatzonmodell --> relation between reall and potential evapotranspiration used in method 4

		char name[500];
		sprintf(name,"%s%d",_Wurzeltiefe_,1);
		WURZ	= gridinfo[getindex(name)].pointer; /* root depth */

		/* initializing a statistc register set for groundwater level [m below surface] */
		zoneGW = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneGW) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{ zoneGW[i].n         = zoneGW[i].RelCellSum = zoneGW[i].sx        = 0;
			zoneGW[i].zonenwert = zone[i].zonenwert;
			zoneGW[i].rel       = zone[i].rel;
		}

		/* initializing a statistc register set for soil moisturte [0...0.5 m below surface] */
		zoneSW = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneSW) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{ zoneSW[i].n         = zoneSW[i].RelCellSum = zoneSW[i].sx        = 0;
			zoneSW[i].zonenwert = zone[i].zonenwert;
			zoneSW[i].rel       = zone[i].rel;
		}

		IRRIGgwstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
		&zoneGW,(int)zonmax,(int)IRRIGgwstatcode,(int)4,(char*)IRRIGgwstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, InV < 60);
		if (retvalue!=0) throw retvalue;
		IRRIGgwstat->schreibzonenwerte("irrigation water amount from groundwater [mm]",instanzname,LastResults,maxLastResults);

		IRRIGswstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
		&zoneSW,(int)zonmax,(int)IRRIGswstatcode,(int)4,(char*)IRRIGswstatfilename,DefaultOutputDirectory,(int)numsubb,(int*)listsubb, InV < 60);
		if (retvalue!=0) throw retvalue;
		IRRIGswstat->schreibzonenwerte("irrigation water amount from  surface water [mm]",instanzname,LastResults,maxLastResults);
		initialized=0; // darf erst nach SetAbfrelAdress() auf 1 gesetzt werden
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

/* the irrigation demand has to be checked against the available water from surface water (rivers, reservoirs). 
So the adresses of the flow-arrays must be known within the irrigation moduler */
void irrigationclass::SetQgesAdresse(double **p_Qges)
{
	if (do_it==0) 
	return;
	else  
	Qges=p_Qges;
}

/* the routing module needs the informations about the irrigation demands */
statrecord *irrigationclass::GetZoneSWAdress() 
{
	if (do_it==0)
	return NULL;
	else
	return zoneSW;
}

void  irrigationclass::SetMaxIrrigAdress(double **p_maxirrig,double **p_maxirrigFromReservoirs)
{
	if (do_it==0) 
	{
		return;
	}
	else  
	{
		maxirrig=p_maxirrig;
		maxirrigFromReservoirs = p_maxirrigFromReservoirs;
	}
	if (!maxirrig)
	{
		(*WSstream) << "FATAL: if irrigation is modelled, also routing has to be modelled!"
		<< "\n(Adress of array maxirrig[] is NULL\n" << "[end-entry]";
		throw(60);
	}
	if (!maxirrigFromReservoirs)
	{
		(*WSstream) << "FATAL: if irrigation is modelled, also routing has to be modelled!"
		<< "\n(Adress of array maxirrigFromReservoirs[] is NULL\n" << "[end-entry]";
		throw(60);
	}
}

void  irrigationclass::SetRtargetAdress(int **p_rtarget)
{
	if (do_it==0) return;
	else  rtarget=p_rtarget;
	if (!rtarget)
	{
		(*WSstream) << "FATAL: if irrigation is modelled, also routing has to be modelled! "
		<< "\n(Adress of array routing-target is NULL)\n" << "[end-entry]";
		throw(60);
	}
}

void  irrigationclass::SetAbfrelAdress(double **p_abfrel)
{
	if (do_it==0) return;
	else  abfrel=p_abfrel;
	if (!abfrel)
	{
		(*WSstream) << "FATAL: if irrigation is modelled, also routing has to be modelled!"
		<< " \n(Adress of array abfrel is NULL)\n" << "[end-entry]";
		throw(60);
	}
}

void irrigationclass::SetRiverflow_tableAdress(riverflow_table_type *adresse)
{
	if (do_it==0) return; else riverflow_table=adresse;
}

void irrigationclass::run(long me,long se,long tg,long mt,long jr)
{
	int i,j,z,ic,tc,flag;
	int actindex,routindex;
	double irrigamount;
	if (do_it == 0) return;
	if (doSpinUp == 1) return; // will not be processed in SpinUp mode
	/* Ausgabe einer Nachricht */
	teilmodellklasse::run(me,se,tg,mt,jr,NULL,NULL);

	// MPI parallelizing is not easily possible for processes which are sharing a common reservoir
	// this is the case here: the source of irrigation water may be a river or a lake (reservoir, i.e.)
	// so all irrigation will be done in instance 0 only
	// to speed the model up, all irrigated cells may be included into MPIREG-value 0 (the master thread).
	// however, this will be a problem when the fraction of irrigated areas is higher than the numbers of processors are
	// because the overall performance would go down with an unproportionally growing region for the master instance
	// at the end of irrigation, precipitation and irrigsw/irriggw must be broadcasted to all instances (or scattered, depends on performance)

	for (z=1;z<=tracer_count;z++) tracer_module[z].stat->init();
	IRRIGgwstat->init();
	IRRIGswstat->init();

	//char *aktname;
	//aktname=new char[200];
	//sprintf(aktname,"Prec_%d",mpi_rank);
	//printgridfile(aktname,P);	

	int precgrid=getindex(_Niederschlag_);
	MPIGatherGrid(precgrid);
	for (z = 1; z <= tracer_count; z++)
	{
		MPIGatherGrid(tracer_module[z].fluxindex);
	}

	//if (mpi_rank==0)
	//{
	//	sprintf(aktname,"PrecGathered_%d",mpi_rank);
	//	printgridfile(aktname,P);	
	//}

	//if (aktname) free(aktname); aktname=NULL;

	// for test only: delete P and see, if it will be set to valid values in MPIScatterGrid again
	//  if (mpi_rank>0)
	//  {
	//	  for (i=1;i<=nrows;i++)
	//	  {
	//		for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
	//		{
	//		  if (zn[i][j]>nodata)
	//		  {
	//			  P[i][j]=0;
	//		  }
	//		}
	//	  }
	//	char *aktname;
	//	aktname=new char[200];
	//	sprintf(aktname,"Prec0gesetzt_%d",mpi_rank);
	//	printgridfile(aktname,P);	
	//	if (aktname) free(aktname); aktname=NULL;
	//  }

	if (mpi_rank==0)
	{
		for (i=1;i<=nrows;i++)
		{
			for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			{
				if (zn[i][j]>nodata)
				{
					TG=zarray[(int)zn[i][j]];
					ic=(int)IRRIGCODE[i][j];
					/* now calculating irrigation water amount for the actual cell, depending on the settings in the irrigation table or on the actual soil moisture */
					if ((ic>0)&&(ic<=irrigationtable[0].Code)) /* only if this type of land use is present within the grid containing irrigation codes */
					{
						if (irrigationtable[ic].methode>0) /* only if this crop type should be irrigated within this model run (method 1: irrig. table, 2: demand optimizing irrigation) */
						{
							/* calculating the actual irrigation amount */
							if (irrigationtable[ic].methode==4) { // according to demand goven by shortfall of the relation ETR/ETP below 0.8
								if (ETR2ETP[i][j]<irrigationtable[ic].MinEtr2Etp) 
								{
									if (IRRIG[i][j]==0) /* otherwise the value from the last interval is still valid until SBHIGH is reached */
									{
										IRRIG[i][j]=(SBHIGH[i][j]-SB[i][j])*WURZ[i][j]*1000.0; /* irrigation water demand in mm */
										IRRIG[i][j]*=(InV/1440); /* if timestep == 1d -> InV=1440 min -> factor = 1; if timeste==1h -> InV=60min, factor=1/24 */
									}
								}
								else
								{
									if (InV<1440)
									{
										if (SB[i][j]>=SBHIGH[i][j]) IRRIG[i][j]=0.0;  /* end of irrigation, if maximum level for irrigation is reached */
									}
									else IRRIG[i][j]=0.0; /* in resolution 1 day all irrigation water is applied within the same intervall, IRRIG can be set to zero now */
								}
							} else if (irrigationtable[ic].methode==3) /* according to demand */ {
								if (SB[i][j]<SBLOW[i][j]) /* if the soil within the root zone is below lower limit for irrigation */
								{
									/* if IRRIG[i][j] is != 0 this means, that there is actual a filling in progress (needs 24 hours if running in hourly resolution) */
									if (IRRIG[i][j]==0) /* otherwise the value from the last interval is still valid until SBHIGH is reached */
									{
										IRRIG[i][j]=(SBHIGH[i][j]-SB[i][j])*WURZ[i][j]*1000.0; /* irrigation water demand in mm */
										IRRIG[i][j]*=(InV/1440); /* if timestep == 1d -> InV=1440 min -> factor = 1; if timeste==1h -> InV=60min, factor=1/24 */
									}
								}
								else
								if (InV<1440)
								{
									if (SB[i][j]>=SBHIGH[i][j]) IRRIG[i][j]=0.0;  /* end of irrigation, if maximum level for irrigation is reached */
								}
								else IRRIG[i][j]=0.0; /* in resolution 1 day all irrigation water is applied within the same intervall, IRRIG can be set to zero now */
							} else if ((irrigationtable[ic].methode==1)||(irrigationtable[ic].methode==2)) { /* data from irrigation table */
								IRRIG[i][j]=0; /* if irrigation is demanded, the volume for irigation is calculated in the following rows */
								for (tc=1;tc<=irrigationtable[ic].tabcnt;tc++) /* tc: table-count */
								{
									if ((irrigationtable[ic].month[tc]==mt)&&(irrigationtable[ic].day[tc]==tg)) /* an irrigation day is reached */
									{
										IRRIG[i][j]=irrigationtable[ic].amount[tc]; /* amount from irrigation table in control file in mm/day */
										IRRIG[i][j]*=(InV/1440); /* if timestep == 1d -> InV=1440 min -> factor = 1; if timeste==1h -> InV=60min, factor=1/24 */
										tc=irrigationtable[ic].tabcnt+1;
									}
								}
							}

							if (irrigationtable[ic].origin==1) /* irrigation from groundwater */
							{
								/* comparing with maxpump, if irrigation water source is groundwater */
								if (IRRIG[i][j]>MAXPUMP[i][j])
								IRRIG[i][j]=MAXPUMP[i][j]; /* reducing irrigation to the maximum amount of pumpable water */
								IRRIGgwstat->wertzufuegen(i,j,TG,IRRIG[i][j]); /* Statistics for irrigation water demand from groundwater */
								IRRIGswstat->wertzufuegen(i,j,TG,0.0); /* Statistics for irrigation water demand from surface water */
							}
							else 
							if (irrigationtable[ic].origin==2) /* irrigation from surface water */
							{
								IRRIGgwstat->wertzufuegen(i,j,TG,0.0); /* Statistics for irrigation water demand from groundwater */
								IRRIGswstat->wertzufuegen(i,j,TG,IRRIG[i][j]); /* Statistics for irrigation water demand from surface water */
							}
						}
						else 
						{
							IRRIG[i][j]=0.0; /* if this land use type should not be irrigated (method=0 in control file, but code is present in irrigation codes grid) */
							IRRIGgwstat->wertzufuegen(i,j,TG,0.0); /* Statistics for irrigation water demand from groundwater */
							IRRIGswstat->wertzufuegen(i,j,TG,0.0); /* Statistics for irrigation water demand from surface water */
						}
					}
					else 
					{
						IRRIG[i][j]=0.0; /* if no irrigationcode is present for this cell */
						IRRIGgwstat->wertzufuegen(i,j,TG,0.0); /* Statistics for irrigation water demand from groundwater */  
						IRRIGswstat->wertzufuegen(i,j,TG,0.0); /* Statistics for irrigation water demand from surface water */
					}
				}
				else IRRIG[i][j]=nodata;
			}
		}

		/*
		calculating the fraction of irrigation water coming from reservoirs and from rivers.
		Reservoirs have first priority, rivers have second priority. Only if a reservoir
		falls dry, then the rivers will be taken into account (which are feed by inflowing 
		tributaries then, only.
	*/

		double *usedIrrigationFromReservoirs;
		usedIrrigationFromReservoirs = new double[zonmax+1];
		for (i=1;i<=zonmax;i++) 
		{

			irrigamount=zoneSW[i].sx/zoneSW[i].RelCellSum;
			if (maxirrigFromReservoirs[0][i]>0)
			{
				if (irrigamount<=maxirrigFromReservoirs[0][i])
				{
					// there was enough water for irrigation stored in the reservoir 
					usedIrrigationFromReservoirs[i]=irrigamount;
				}
				else
				{   
					// all water in the reservoir was used for irrigation
					usedIrrigationFromReservoirs[i]=maxirrigFromReservoirs[0][i];
				}
			}
			else
			{
				usedIrrigationFromReservoirs[i]=0;
			}
		}

		/* subtracting irrigation amount from maxirrig values in routing subbasins 
		This must be done first before calculating the relative amount of irrigation water
		available for each subbasin, because maxirrig in routing basins depends on the 
		irrigation amount in the tributaries, too
	*/
		for (i=1;i<=zonmax;i++) 
		{
			routindex=rtarget[1][i];
			actindex=routindex;
			if (routindex>0)
			{
				irrigamount=zoneSW[routindex].sx/zoneSW[routindex].RelCellSum-usedIrrigationFromReservoirs[routindex];
				if (irrigamount>maxirrig[0][actindex]) irrigamount=maxirrig[0][actindex]; /* the larger value of both variables limits the real irrigation amount */
				if (irrigamount<0) irrigamount=0;

				if (riverflow_table)
				{
					if (riverflow_table[actindex].ROUT_index!=0)
					{
						riverflow_table[actindex].flowm3s-=irrigamount*zoneSW[0].rel*zoneSW[actindex].rel/(0.06*InV); /* reduction by m3/s */
						riverflow_table[actindex].flow_mm[0]=riverflow_table[actindex].flowm3s / riverflow_table[actindex].upstream_area * InV * 0.06;
						if (riverflow_table[actindex].flow_mm[0]<0.0)
						{
							riverflow_table[actindex].flow_mm[0]=0.0;
							riverflow_table[actindex].flowm3s=0.0;
						}
					}
				}

				while (rtarget[1][routindex]>0) /* as long as the subbasin has a target subbasin for routing, the irrigation amount is subtracted from the maxiirig value of the target subbasin */
				{
					routindex=rtarget[1][routindex];
					if (routindex>0)
					{
						maxirrig[0][routindex]-=irrigamount*zoneSW[actindex].rel/zoneSW[routindex].rel;
						if (maxirrig[0][routindex]<0) maxirrig[0][routindex]=0;
						if (riverflow_table)
						{
							if (riverflow_table[routindex].ROUT_index!=0)
							{
								riverflow_table[routindex].flowm3s-=(zoneSW[actindex].sx/zoneSW[actindex].RelCellSum-usedIrrigationFromReservoirs[actindex])*zoneSW[0].rel*zoneSW[actindex].rel/(0.06*InV); /* reduction by m3/s */
								riverflow_table[routindex].flow_mm[0]=riverflow_table[routindex].flowm3s / riverflow_table[routindex].upstream_area * InV * 0.06;
								if (riverflow_table[routindex].flow_mm[0]<0.0)
								{
									riverflow_table[routindex].flow_mm[0]=0.0;
									riverflow_table[routindex].flowm3s=0.0;
								}
							}
						}
					}
				}
			}
		}

		/* comparing irrigation water demand from surface water with maxirrig (in routing module calculated max. available irrigation water from rivers/reservoirs) */
		flag=0; /* if flag is set to 1, than irrigation from surface water could'nt be satisfied */
		for (i=1;i<=zonmax;i++)
		{
			irrigamount=zoneSW[i].sx/zoneSW[i].RelCellSum;
			if ((zoneSW[i].sx>0.0)&&(((double)(zoneSW[i].sx/zoneSW[i].RelCellSum)-(double)(usedIrrigationFromReservoirs[i]))>1e-8))
			{
				maxirrig[0][i]= maxirrig[0][i]/((double)(zoneSW[i].sx/zoneSW[i].RelCellSum)-(double)(usedIrrigationFromReservoirs[i])); /* maxirrig[0] contains now for each subbaisn the reduction factor for all irrigated cells */
			}
			else
			{
				maxirrig[0][i]=1.0; /* reduction factor is set to 1.0 if no irrigation is required */
			}

			if (maxirrig[0][i]>=1.0) 
			{
				maxirrig[0][i]=1.0; 
			}
			else 
			{
				flag=1;
			}
			/* if flag is still 0 we have to check if the reservoirs may be the limiting factor for irrigation */
			if (irrigamount>maxirrigFromReservoirs[0][i])
			{
				flag+=2;
			}

		}

		// only if irrgation amount is limited by surface water (rivers as well as reservoirs) 
		// availability and therefore a reduction has to be applied */
		if (flag>0) // 1: river water was limiting; 2: reservoir water was limiting; 3: both was limiting 
		{
			for (i=1;i<=zonmax;i++) { zoneSW[i].sx=0; zoneSW[i].n=0; zoneSW[i].RelCellSum=0;}
			/* reducing irrigation from surface water by mean reduction factors for subbasins */
			double dFromReservoir;
			double dMaxFromResrvoir;
			double dAvailableFractionFromRes;
			double dFractionFromRes;
			double dFractionFromRiv;
			for (i=1;i<=nrows;i++)
			{
				for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
				{
					if (zn[i][j]>0)
					{
						TG=zarray[(int)zn[i][j]];
						dFromReservoir=usedIrrigationFromReservoirs[TG];
						dMaxFromResrvoir=maxirrigFromReservoirs[0][TG];
						if (maxirrigFromReservoirs[0][TG]>0)
						{
							dAvailableFractionFromRes=usedIrrigationFromReservoirs[TG]/maxirrigFromReservoirs[0][TG];
						}
						else
						{
							dAvailableFractionFromRes=0;
						}
						if (dAvailableFractionFromRes>1.0)
						{
							dAvailableFractionFromRes=1.0;
						}
						if ((maxirrigFromReservoirs[0][TG]+maxirrig[0][TG])>0)
						{
							dFractionFromRes=maxirrigFromReservoirs[0][TG]/(maxirrigFromReservoirs[0][TG]+maxirrig[0][TG]);
						}
						else
						{
							dFractionFromRes = 0;
						}
						dFractionFromRiv=1.0-dFractionFromRes;

						IRRIG[i][j] = dFractionFromRiv*IRRIG[i][j]*maxirrig[0][TG] + IRRIG[i][j]*dFractionFromRes*dAvailableFractionFromRes;

						IRRIGswstat->wertzufuegen(i,j,TG,IRRIG[i][j]);
					}
				}
			}
		}
		for (i=1;i<=zonmax;i++)
		{
			maxirrigFromReservoirs[0][i]=usedIrrigationFromReservoirs[i];
		}
		/* now adding the irrigation amount to the precipitation and mixing tracer concentrations(or adding salt masses) */
		for (i=1;i<=nrows;i++)
		{
			for (j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++)
			{
				if (zn[i][j]>0)
				{
					TG=zarray[(int)zn[i][j]];
					ic=(int)IRRIGCODE[i][j];
					if ((ic>=0)&&(ic<=irrigationtable[0].Code))
					{
						for (z=1;z<=tracer_count;z++)
						{
							if ((IRRIG[i][j]>0)&&(ic>=0))
							{
								if (irrigationtable[ic].origin==1) /* from groundwater */
								tracer_module[z].concflux[i][j]=mixing(tracer_module[z].cinput[i][j],PUMPCONC[z][i][j],P[i][j],IRRIG[i][j]);
								else
								if (irrigationtable[ic].origin==2) /* from surface water */
								tracer_module[z].concflux[i][j]=mixing(tracer_module[z].cinput[i][j],maxirrig[z][TG],P[i][j],IRRIG[i][j]);
							}
							else 
							{
								//              tracer_module[z].concflux[i][j]=tracer_module[z].cinput[i][j]; /* if only precipitation is falling (no irrigation): tracer concentration is get from precipitation only */
								if (P[i][j]>0)
								tracer_module[z].concflux[i][j]=tracer_module[z].cinput[i][j]; /* if only precipitation is falling (no irrigation): tracer concentration is get from precipitation only */
								else
								tracer_module[z].concflux[i][j]=0.0; /* if no irrigation and no precipitation is falling: concentration is set to zero */
							}
							if (RelCellArea!=NULL)
							{
								tracer_module[z].stat->zone[TG].sx+=tracer_module[z].concflux[i][j]*RelCellArea[i][j]*(P[i][j]+IRRIG[i][j]);
								tracer_module[z].stat->zone[TG].n+=(P[i][j]+IRRIG[i][j]);
								tracer_module[z].stat->zone[TG].RelCellSum+=(P[i][j]+IRRIG[i][j])*RelCellArea[i][j];
							}
							else
							{
								tracer_module[z].stat->zone[TG].sx+=tracer_module[z].concflux[i][j]*(P[i][j]+IRRIG[i][j]);
								tracer_module[z].stat->zone[TG].n+=(P[i][j]+IRRIG[i][j]);
								tracer_module[z].stat->zone[TG].RelCellSum+=(P[i][j]+IRRIG[i][j]);
							}
						}
						P[i][j]+=IRRIG[i][j];
					}
				}
			}
		}
		delete[] usedIrrigationFromReservoirs;
	} // mpi_rank==0

	// now P and IRRIG must be send to all instances in order to be able to use correct input values in 
	// snow-, interception-, and soilmodel

	int irrigindex=getindex(_irrigation_amounts_);
	MPIScatterGrid(precgrid);
	MPIScatterGrid(irrigindex);
	for (z = 1; z <= tracer_count; z++)
	{
		MPIScatterGrid(tracer_module[z].fluxindex);
	}
	writegridtofile(irrigindex);

	//  aktname=new char[200];
	//  if (mpi_rank>=0)
	//  {
	//    aktname=new char[200];
	//	  sprintf(aktname,"PrecScattered_%d",mpi_rank);
	//	  printgridfile(aktname,P);	
	//	  sprintf(aktname,"IrrigScattered_%d",mpi_rank);
	//	  printgridfile(aktname,IRRIG);	
	//  }
	//  if (aktname) free(aktname); aktname=NULL;


	IRRIGgwstat->calc_and_done(jahr,monat,tag,stunde, minute);
	IRRIGswstat->calc_and_done(jahr,monat,tag,stunde, minute);
	for (z=1;z<=tracer_count;z++)
	{
		tracer_module[z].stat->calc_and_done(jahr,monat,tag,stunde, minute);
		//    writegridtofile(tracer_module[z].storindex);
		writegridtofile(tracer_module[z].fluxindex);
	}

}

void irrigationclass::flushbuffers()
{
	int z;
	if (doSpinUp == 1) return; // will not be processed in SpinUp mode
	if (do_it==0) return;
	else
	{
		if (IRRIGgwstat) IRRIGgwstat->flushbuffer();
		if (IRRIGswstat) IRRIGswstat->flushbuffer();
		for (z=1;z<=tracer_count;z++)
		tracer_module[z].stat->flushbuffer();
	}
}


void irrigationclass::ClearMemory()
{
	if (doSpinUp == 1) return; // will not be processed in SpinUp mode
	if (do_it>0)
	{
		IRRIGgwstat->ClearMemory(); if (IRRIGgwstat ) delete(IRRIGgwstat ); IRRIGgwstat = NULL;
		IRRIGswstat->ClearMemory(); if (IRRIGswstat ) delete(IRRIGswstat ); IRRIGswstat = NULL;
		if (zoneGW) free(zoneGW); zoneGW = NULL;
		if (zoneSW) free(zoneSW); zoneSW = NULL;
		if (IRRIGgwstatfilename) free(IRRIGgwstatfilename); IRRIGgwstatfilename  = NULL;
		if (IRRIGswstatfilename) free(IRRIGswstatfilename); IRRIGswstatfilename  = NULL;
	}
	modelobject::ClearMemory();
}

void irrigationclass::substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults)
{
	int z;
	char *tempstr;
	if (doSpinUp == 1) return; // will not be processed in SpinUp mode
	if (do_it==0) return;
	else
	{
		teilmodellklasse::substance_transport_ini(retvalue,p_tracer_count,p_tracer_info_table,LastResults,maxLastResults);
		if (tracer_count<=0) return;
		else
		{
			tempstr=newstring(200,sizeof(char));
			for (z=1;z<=tracer_count;z++)
			{
				sprintf(tempstr,"concflux_tracer_%1d_%s",z,"input"); /* if irrigation is considered, tracer concentrations are get from irrigation submodel */
				tracer_module[z].cinput=gridinfo[getindex(tempstr)].pointer;
			}
			if (tempstr) free(tempstr); tempstr=NULL;
			return;
		}
	}
}

void irrigationclass::GetPumpingConcentrations()
{
	int z;
	char *identify;
	if (doSpinUp == 1) return; // will not be processed in SpinUp mode
	if (do_it==0) return;
	else
	{
		if (tracer_count<=0) return;
		else
		{
			PUMPCONC=(double***)calloc(tracer_count+1,sizeof(double**));
			identify=(char*)calloc(200,sizeof(char));
			for (z=1;z<=tracer_count;z++)
			{
				sprintf(identify,"concstor_tracer_%d_%s",z,_Unsatzonmodell_);
				PUMPCONC[z]=gridinfo[getindex(identify)].pointer;
			}
			if (identify) free(identify); identify=NULL;
		}
	}
}

void irrigationclass::substance_transport_delete()
{
	if (doSpinUp == 1) return; // will not be processed in SpinUp mode
	if (do_it==0) return;
	else
	{
		teilmodellklasse::substance_transport_delete();
		if (tracer_count<=0) return;
		else
		{
			if (PUMPCONC) free(PUMPCONC); PUMPCONC=NULL;
		}
	}
}


