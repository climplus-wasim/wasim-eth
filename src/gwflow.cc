/*   $Revision: 49 $ 
	$Date: 5.02.11 22:53 $  */

/*
(c) 1995-2007 by Joerg Schulla

groundwater flow and transport module
finite differences with equal grid cells, layered version (max. 3 layers)
no dry-falling of cells allowed, 
The module gets the recharge/evapotranspiration loss from unsaturated zone as
well as the changed groundwater levels and calculates new flow and heads

	more version history: see source safe comments

*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <csignal>
#include <sstream>
#include <stddef.h>
#include "stringp.h"
#include "modconst.h"
#include "modef.h"
#include "teilmod.h"
#include "wasim.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _MPI
#include <mpi.h>
#endif

using namespace std;


void gwflowclass::flushbuffer(FILE *datei,char *dateiname,unsigned long *bufcnt,buffertype *buffer)
/* flushing the file buffer */
{
	int i;
	unsigned long j;

	if (*bufcnt>0)
	{
		signal(SIGINT,SIG_IGN);
#if defined _BCC_
		signal(SIGBREAK,SIG_IGN);
#endif
		(*WSstream) << "write "<<dateiname;
		if ((datei=fopen(dateiname,"a+"))==NULL)
		{
			(*WSstream) << "FATAL: error opening file "<<dateiname<<"[end-entry]";
			throw(65);
		}
		for (int j=1;j<=*bufcnt;j++)
		{
			fprintf(datei,"%ld\t%ld\t%ld\t%ld",buffer[j].j,buffer[j].m,buffer[j].t,buffer[j].s);
			for (i=1;i<=layers;i++)
			fprintf(datei,"\t%g\t%g\t%g\t%g\t%g\t%g",buffer[j].wert[6*(i-1)+1],buffer[j].wert[6*(i-1)+2],buffer[j].wert[6*(i-1)+3]
			,buffer[j].wert[6*(i-1)+4],buffer[j].wert[6*(i-1)+5],buffer[j].wert[6*(i-1)+6]);
			fprintf(datei,"\n");
		}
		fclose(datei);
		(*WSstream) << "("<<(long)((*bufcnt)*((long)4*sizeof(long)+(6*layers+1)*sizeof(buffer[1].wert[1])))
		<<" of "<<(long)buffersize<<" buffer-bytes)"<<"[end-entry]";
		*bufcnt=0;
		SetSignalAbfangfunktion();
	}
}

// initialize lakes is mainly for makeing sure, the POND grid is used as it was possible before, but not required.
// the routine itself seems a little bit over dimensioned, but it's easier to extend later, when additions may be required
void gwflowclass::initialize_lakes(std::string instanzname)
{
	char *tstr;

	if (!ist_vorhanden(_LAKE_CODES_))
	{
		(*WSstream) << "INFO: no LAKE-grid was found in "<<this->instanzname<<"\n"<<
		"Lake modelling is switched off for "<<this->instanzname<<".\n" << "[end-entry]";
		//	throw(223);
	}
	else
	{

		tstr = new char[_stringsize_];

		if ((steuerdatei=fopen(steuerdateiname,"rt"))==NULL)
		{
			(*WSstream) << "FATAL: error opening "<<steuerdateiname<<"[end-entry]";
                        fclose(steuerdatei);
			throw(28);
		}

		strcpy(tstr,instanzname.c_str());

		/* suchen des richtigen Abschnittes in der Steuerdatei anhand des Namens der Meteogroesse */
		findeabschnitt(tstr);

		fscanfln(steuerdatei,"%d",&lake_active,"flow control parameter for ",tstr);
		// switch off lake model when in Spinup-mode
		if (doSpinUp) {
			lake_active=0;
		} 

		if (lake_active != 0)
		{
			if (!iUsePond)
			{
				(*WSstream) << "FATAL: \n=======================================\n"
				<< "!!!! ERROR !!!! \n\n"
				<< "Error in initialization method for "<<instanzname<<".\n"
				<< "pond_use must be set in "<<this->instanzname<<"\n"
				<< "for using the "<<instanzname<<
				"\n\n=======================================\n"<< "[end-entry]";
                                fclose(steuerdatei);
				throw(214);
			}

			fscanfln(steuerdatei,"%d",&lake_ini_method,"initialization method for ",tstr);
			if ((fclose(steuerdatei))==EOF)
			{
				(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"after reading initialization method for "<<instanzname<<"[end-entry]";
                                fclose(steuerdatei);
				throw(11);
			}

			LAKE=gridinfo[getindex(_LAKE_CODES_)].pointer;
			if (ist_vorhanden(_ponding_storage_))
			{
				POND=gridinfo[getindex(_ponding_storage_)].pointer;
			}
			else
			{
				(*WSstream) << "FATAL: \n=======================================\n"<<
				"!!!! ERROR !!!! \n\n"<<
				"\nLAKE-grid was found, but no POND-Grid. For using lakes,"<<
				"\nboth grids must be present."<<
				"\n\n=======================================\n"<< "[end-entry]";
                                fclose(steuerdatei);
				throw(28);
			}

			if ((lake_ini_method<1)||(lake_ini_method>2))
			{
				(*WSstream) << "FATAL: \n=======================================\n"
				<< "!!!! ERROR !!!! \n\n "
				<< "Error in initialization method for " << instanzname << ".\n"
				<< "Lake initialization method (" << lake_ini_method << ") not valid.\n"
				<< "Allowed values are 1 (no change of grids) or 2 (adjust DEM, SLOPE, DRAIN).\n"
				<< "In either case, theta_soil will be set to saturation,\n"
				<< "if the cell is marked as lake and POND is > 0. Also, the groudwater level\n"
				<< "will be set to 0." 
				<< "\n\n=======================================\n"<< "[end-entry]";
                                fclose(steuerdatei);
				throw(214);
			}

			// change some values as initialization for lake use
			if (lake_ini_method==1) // only change saturation values
			{
				for (int i=mpistartrow;i<=mpiendrow;i++)
				{
					for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
					{
						if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
						{
							if (LAKE[i][j]>=0)
							{
								// no special initializations required for now
							}
						}
					}
				}
			}
			else if (lake_ini_method==2)
			{
				for (int i = mpistartrow; i <= mpiendrow; i++)
				{
					for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
					{
						if ((zn[i][j] != nodata) && (MPIREG[i][j] == mpi_rank))
						{
							if (LAKE[i][j] >= 0)
							{
								// no special initializations required for now
							}
						}
					}
				}
			}
		}
		delete[] tstr;
//                fclose(steuerdatei);
	}
}


//#pragma optimize ("p",on)
gwflowclass::gwflowclass(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults,int maxLastResults,double breitearg,
	double laengearg,double bezugsmeridianarg,const char *bezeichnung,char *steuername,int LogLevel,
	int progaus,double p_min_to_int_end,int p_numsubb,int *p_listsubb,
	gridinfotype *p_gridinfo,int *p_zarray,statrecord *p_zone,firstlast *p_FirstLastCol ,
	int p_zonmax, int p_zonmaxwert, int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow, double p_xll, double p_yll,
	double p_cellsize,double p_nodata,
	long mine,long ste , long tge , long mte , long jre, double p_timekorr,
	LUSEtabletype *p_LUtable,SOILtabletype *p_SOILtable,MultiLayerLUSEType *p_MultiLayerLUTable,irrigationtabletype *p_irrigationtable,
	bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
	long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
	int p_doSpinUp,bool p_routingmode,std::vector<HRUMapElement*> *p_vHRU,int substepcnt):
teilmodellklasse(Stream,breitearg,laengearg,bezugsmeridianarg,bezeichnung,steuername,LogLevel,progaus,p_min_to_int_end,p_numsubb,p_listsubb,
	p_gridinfo,p_zarray,p_zone,p_FirstLastCol,p_zonmax,p_zonmaxwert,
	p_ncols, p_nrows, p_mpistartrow, p_mpiendrow, p_xll, p_yll, p_cellsize, p_nodata, mine,ste, tge, mte, jre, p_timekorr, p_LUtable, p_SOILtable, p_MultiLayerLUTable,
	WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory, v_mpiregioncnt, p_sendHalo, p_recvHalo,
	p_doSpinUp, p_routingmode,p_vHRU,substepcnt),
DHM(0),GWN(0),GWSTAND(0),SOIL   (0),DELTAGW(0),BOUNDH(0),BOUNDQ(0),HEAD(0),HEADold(0),HEADolditer(0),BALANCE(0),THICK(0),
KX(0),KY(0),TRANSX(0),TRANSY(0),FLOWX(0),FLOWY(0),S0(0),LEAKAGE(0),UPLIM(0),KOLM(0),HEADoldold(0),ICE(0),POND(0),
glogbuffer(0),glogbufmax(0),glogbufcnt(0),glogfile(0),glogfilename(0),control_row(0),control_col(0),
layers(0),auto_sor(0),s(0),num_oscill(0),dx(0),dy(0),alpha(0),sor(0),sor_sum(0),sor_anz(0),sum_diff(0),
maxsum_diff(0),maxerror(0),maxerrorf(0),errorsum(0),errorcount(0),meanerror(0),meanerrorf(0),mass0(0),
massflux1(0),massflux2(0),massflux3(0),massflux4(0),massflux5(0),massflux6(0),iter(0),maxiter(0),FLOWXgridname(0),
FLOWYgridname(0),
HEADgridname(0),HEADoldgridname(0),HEADoldoldgridname(0),HEADolditergridname(0),
FLOWXwritegrid(0),FLOWYwritegrid(0),
HEADwritegrid(0),HEADoldwritegrid(0),HEADoldoldwritegrid(0),HEADolditerwritegrid(0),
transxindex(0),transyindex(0),
flowxindex(0),flowyindex(0),leakageindex(0),uplimindex(0),
headindex(0),headoldindex(0),headolditerindex(0),headoldoldindex(0),
Qges(0),riverflow_table(0),TG(0),methode(0),iUsePond(0),iPondAvailable(0),
balanceindex(0), BALANCEwritegrid(0), QDSU(0), THAWDEPTH(0), FROZENDEPTH(0), thawdepthindex(-1), frozendepthindex(-1),frozenDepthScaling(1.0)
{
	try 
	{
		retvalue=0;
//		int i,j,k,lr;
		char *tempstr;
		unsigned long bufrow;


		if (routingmode)
		{
			do_it = 0; // wether the submodels code for do_it is set to 0 or 1: when routingmode is active, this submodel should never run
		}

		if (do_it == 0) return;

		// gw model should be inactive when the model is used in spinup mode since the communication between 
		// unsatzon model and groundwater model will be as expected onm the HRU master cells only and all other cells 
		// will not exchange any information.
		// This will inevitably lead to balance problems (e.g. high contant direct flwo for saturated cells ot other phenomens)
		// so it is the only option to switch off groundwater model for spinup mode.
		// it is, however, required to kind of calibrate the baseflow in [unsatzonmodel] in order to get 
		// at least some realistic groundwater level.
		if (doSpinUp!=0 && HRU==NULL)
		{
			do_it=0;
			return;
		}

		if (MPIREG==NULL)
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

		int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow; // sr = start row for mpi with halo
		int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow; // er = end row for mpi with halo

		irrigationtable = p_irrigationtable;
#ifdef __TIME__
		(*WSstream) << "file "<<__FILE__<<", compiled at "<<__DATE__<<", "<<__TIME__<<"[end-entry]";
#endif
		/* read method Parameter methode */
		fscanfln(steuerdatei,"%d",&methode       ,"method in "                  ,tstr);
		/* read maximum iteration steps */
		fscanfln(steuerdatei,"%d",&maxiter       ,"maxiter in "                 ,tstr);
		/* read maxerror_allowed */
		fscanfln(steuerdatei,"%lf",&maxerrorf    ,"maxerror_allowed in "        ,tstr);
		/* read Crank_Nicholson Parameter alpha */
		fscanfln(steuerdatei,"%lf",&alpha        ,"alpha in "                   ,tstr);
		/* read SOR Parameter sor */
		fscanfln(steuerdatei,"%lf",&sor          ,"sor-parameter in "           ,tstr);
		if (sor<=0) auto_sor=1; else auto_sor=0;
		/* other parameters */
		fscanfln(steuerdatei,"%d",&gridslesen   ,"grid read parameter in "      ,tstr);

		/* read number of layers */
		fscanfln(steuerdatei,"%d",&layers        ,"number of layers in "        ,tstr);
		/* here the row/column indices of a control location is given.
	For this location, the three following files are generated */

		glogfilename       = (char*)calloc(_stringsize_,sizeof(char));

		fscanfln(steuerdatei,"%d%d",&control_row,&control_col   ,"control_row/column for control point fluxes in ",tstr);
		fscanfln(steuerdatei,"%s",glogfilename        ,"name of the control point file  ",tstr);

		ExtendFileName(DefaultOutputDirectory,glogfilename);

		/* row/column-indizes are wrong */
		if ((control_row<0) || (control_row>nrows))
		{
			(*WSstream) << "FATAL: in modul " << bezeichnung << ":\n"
				<< "row index (" << control_row << ") for output-cell must be between " << 1 << " and " << nrows << " (Input as row col) or 0 0 " << "[end-entry]";
			throw(8);
		}
		if ((control_col<0) || (control_col>ncols))
		{
			(*WSstream) << "FATAL: in modul " << bezeichnung << ":\n"
				<< "column index (" << control_col << ") for output-cell must be between " << 1 << " and " << ncols << " (Input as row col) " << "[end-entry]";
			throw(9);
		}

		// try to read iUsePond: if next line does not start with a number 0 or 1, then the flag is not available
		long fposition=ftell(steuerdatei);
		iUsePond=fgetc(steuerdatei);
		if ((iUsePond!=int('0'))&&(iUsePond!=int('1')))
		{
			iUsePond=0;
			iPondAvailable=0;
			fseek(steuerdatei,fposition,0);
			(*WSstream) << "INFO: in modul "<<bezeichnung<<": no POND-Grid will be used because Use-Pond-Flag is not present" << "[end-entry]";
			pondindex=-1;
		}
		else
		{
			if (iUsePond==int('1'))
			{
				iUsePond=1;
				tstr=fgets(tstr,_stringsize_,steuerdatei);
				if (ist_vorhanden(_ponding_storage_))
				{
					iPondAvailable=1;
					pondindex=getindex(_ponding_storage_);
					POND=gridinfo[pondindex].pointer;
					(*WSstream) << "INFO: in modul "<<bezeichnung<<": Use-Pond-Flag is " << iUsePond << "\n POND-Grid was found --> PONDs will be used in GW-Model" << "[end-entry]";
				}
				else
				{
					iPondAvailable=0;
					pondindex=-1;
					POND=NULL;
					(*WSstream) << "INFO: in modul "<<bezeichnung<<": Use-Pond-Flag is " << iUsePond << "\n POND-Grid was not found --> no PONDs will be used in GW-Model" << "[end-entry]";
				}
			}
			else
			{
				iUsePond=0;
				iPondAvailable=0;
				pondindex=0;
				(*WSstream) << "INFO: in modul "<<bezeichnung<<": Use-Pond-Flag is " << iUsePond << "; --> no PONDs will be used in GW-Model" << "[end-entry]";
			}
		}

		HEADgridname = (char**)calloc(layers + 1, sizeof(char*));
		HEADoldgridname     = (char**)calloc(layers+1,sizeof(char*));
		HEADoldoldgridname  = (char**)calloc(layers+1,sizeof(char*));
		HEADolditergridname = (char**)calloc(layers+1,sizeof(char*));
		FLOWXgridname		= (char**)calloc(layers+1,sizeof(char*));
		FLOWYgridname		= (char**)calloc(layers+1,sizeof(char*));
		BALANCEgridname		= (char**)calloc(layers+1,sizeof(char*));
		headindex			= (int*)calloc(layers+1,sizeof(int));
		headoldindex		= (int*)calloc(layers+1,sizeof(int));
		headoldoldindex		= (int*)calloc(layers+1,sizeof(int));
		headolditerindex	= (int*)calloc(layers+1,sizeof(int));
		balanceindex		= (int*)calloc(layers+1,sizeof(int));

		transxindex     = (int*)calloc(layers+1,sizeof(int));
		transyindex     = (int*)calloc(layers+1,sizeof(int));
		flowxindex      = (int*)calloc(layers+1,sizeof(int));
		flowyindex      = (int*)calloc(layers+1,sizeof(int));
		leakageindex    = (int*)calloc(layers+1,sizeof(int));
		uplimindex      = (int*)calloc(layers+1,sizeof(int));
		FLOWXwritegrid  = (int*)calloc(layers+1,sizeof(int));
		FLOWYwritegrid  = (int*)calloc(layers+1,sizeof(int));
		HEADwritegrid   = (int*)calloc(layers+1,sizeof(int));
		HEADoldwritegrid   = (int*)calloc(layers+1,sizeof(int));
		HEADoldoldwritegrid   = (int*)calloc(layers+1,sizeof(int));
		HEADolditerwritegrid   = (int*)calloc(layers+1,sizeof(int));
		BALANCEwritegrid= (int*)calloc(layers+1,sizeof(int));

		TRANSX         = (double***)calloc(layers+1,sizeof(double**));
		TRANSY         = (double***)calloc(layers+1,sizeof(double**));
		BOUNDH         = (double***)calloc(layers+1,sizeof(double**));
		BOUNDQ         = (double***)calloc(layers+1,sizeof(double**));
		HEAD           = (double***)calloc(layers+1,sizeof(double**));
		HEADold        = (double***)calloc(layers+1,sizeof(double**));
		HEADoldold     = (double***)calloc(layers+1,sizeof(double**));
		HEADolditer    = (double***)calloc(layers+1,sizeof(double**));
		BALANCE        = (double***)calloc(layers+1,sizeof(double**));
		
		THICK          = (double***)calloc(layers+1,sizeof(double**));
		KX             = (double***)calloc(layers+1,sizeof(double**));
		KY             = (double***)calloc(layers+1,sizeof(double**));
		FLOWX          = (double***)calloc(layers+1,sizeof(double**));
		FLOWY          = (double***)calloc(layers+1,sizeof(double**));
		S0             = (double***)calloc(layers+1,sizeof(double**));
		LEAKAGE        = (double***)calloc(layers+1,sizeof(double**));
		UPLIM          = (double***)calloc(layers+1,sizeof(double**));
		KOLM           = (double***)calloc(layers+1,sizeof(double**));
		/* read names for each layer */


		for (int lr=1;lr<=layers;lr++)
		{
			HEADgridname[lr]         = newstring(_stringsize_+1,1);
			HEADoldgridname[lr]      = newstring(_stringsize_+1,1);
			HEADoldoldgridname[lr]   = newstring(_stringsize_+1,1);
			HEADolditergridname[lr]  = newstring(_stringsize_+1,1);
			FLOWXgridname[lr]        = newstring(_stringsize_+1,1);
			FLOWYgridname[lr]        = newstring(_stringsize_+1,1);
			BALANCEgridname[lr]      = newstring(_stringsize_+1,1);
			fscanfln(steuerdatei,"%s",HEADgridname[lr]         ,"name of grid with hydraulic heads in "    ,tstr);

			std::string temp(HEADgridname[lr]);
			if (temp.find_last_of(".")!=std::string::npos)
			{
				sprintf(HEADoldgridname[lr],"%s_old%s",temp.substr(0,temp.find_last_of(".")).c_str(),temp.substr(temp.find_last_of(".")).c_str());
				sprintf(HEADoldoldgridname[lr],"%s_oldold%s",temp.substr(0,temp.find_last_of(".")).c_str(),temp.substr(temp.find_last_of(".")).c_str());
				sprintf(HEADolditergridname[lr],"%s_olditer%s",temp.substr(0,temp.find_last_of(".")).c_str(),temp.substr(temp.find_last_of(".")).c_str());
			}
			else
			{
				sprintf(HEADoldgridname[lr],"%s_old",temp.c_str());
				sprintf(HEADoldoldgridname[lr],"%s_oldold",temp.c_str());
				sprintf(HEADolditergridname[lr],"%s_olditer",temp.c_str());
			}

			fscanfln(steuerdatei,"%d",&(HEADwritegrid[lr])     ,"Writegrid (hydraulic heads) in "          ,tstr);
			HEADoldwritegrid[lr]=HEADwritegrid[lr];
			HEADoldoldwritegrid[lr]=HEADwritegrid[lr];
			HEADolditerwritegrid[lr]=HEADwritegrid[lr];
			fscanfln(steuerdatei,"%s",FLOWXgridname[lr]        ,"name of grid with x-transmissivities in " ,tstr);
			fscanfln(steuerdatei,"%d",&(FLOWXwritegrid[lr])    ,"Writegrid x-transmissivities in "         ,tstr);
			fscanfln(steuerdatei,"%s",FLOWYgridname[lr]        ,"name of grid with y-transmissivities in " ,tstr);
			fscanfln(steuerdatei,"%d",&(FLOWYwritegrid[lr])    ,"Writegrid y-transmissivities in "         ,tstr);
			fscanfln(steuerdatei,"%s",BALANCEgridname[lr]      ,"name of grid with balances in " ,tstr);
			fscanfln(steuerdatei,"%d",&(BALANCEwritegrid[lr])    ,"Writegrid balances in "         ,tstr);
		}

		if ((fclose(steuerdatei))==EOF)
		{
			(*WSstream) << "FATAL: error closing "<<steuerdateiname<<"[end-entry]";
			throw(11);
		}

		if ((steuerdatei=fopen(steuerdateiname,"rt"))==NULL)
		{
			(*WSstream) << "FATAL: error opening "<<steuerdateiname<<"[end-entry]";
			throw(28);
		}
		strcpy(tstr,"[");
		strcat(tstr,_groundwater_tuning_);
		strcat(tstr,"]");
	
		/* suchen des richtigen Abschnittes in der Steuerdatei anhand des Namens der Meteogroesse */
		if (findeabschnitt(tstr,false))
		{
			bool DoContinue=true;
			do
			{
				fgets(tstr,_stringsize_,steuerdatei);
				std::string s(tstr);
				s=s.substr(s.find_first_not_of(" \t"));
				sprintf(tstr,"%s",s.c_str());
				if ((!feof(steuerdatei)) && (index(tstr, "\n")>0) && (index(tstr, "\r")>0)) // skip all empty lines
				{
					DoContinue=true;
					std::string parameter(tstr);
					if (parameter.find("frozenDepthScaling")!=std::string::npos) {
						std::stringstream ss;
						parameter=parameter.substr(parameter.find("=")+1);
						ss << parameter.c_str();
						ss >> frozenDepthScaling;
					}
				}
				else 
				{
					DoContinue=false;
				}
			} while (DoContinue);

		}

		if (fclose(steuerdatei)==EOF)
		{
			(*WSstream) << "FATAL: control file \""<<steuerdateiname<<"\" cannot be closed"<<"[end-entry]";
			throw(11);
		}

		/* getting adresses of matrices already stored in memory */
		DHM     = gridinfo[getindex(_Hoehenmodell_)].pointer;
		GWN     = gridinfo[getindex(_Grundwasserneubildung_)].pointer;
		
		gwstandindex = getindex(_Grundwasserstand_);
		GWSTAND = gridinfo[gwstandindex].pointer;
		
		SOIL    = gridinfo[getindex(_Bodenarten_)].pointer;                        /* soil types (texture) */
		
		deltagwindex = getindex(_Basisabfluss_);
		DELTAGW = gridinfo[deltagwindex].pointer; /* balance of lateral inflows into and lateral outflows out of a grid cell */

		if (RelCellArea!=NULL)
		{
			CellSizeX=gridinfo[getindex(_CellSizeX_)].pointer;
			CellSizeY=gridinfo[getindex(_CellSizeY_)].pointer;
		}
		/* looking for ICE grid */
		if (ist_vorhanden(_ice_firn_)||ist_vorhanden(_GlacierizedCells_))
		{
			if (ist_vorhanden(_GlacierizedCells_))
			{
				ICE = gridinfo[getindex(_GlacierizedCells_)].pointer;
			}
			else
			{
				ICE = gridinfo[getindex(_ice_firn_)].pointer;
			}
		}
		else
		{
			ICE=NULL;
		}

		tempstr=newstring(_stringsize_+1,1);
		for (int lr=1;lr<=layers;lr++)
		{
			sprintf(tempstr,"%s%d",_aquifer_thickness_,lr);
			THICK[lr]  = gridinfo[getindex(tempstr)].pointer;               /* thickness of aquifer lr */
			sprintf(tempstr,"%s%d",_gw_k_x_,lr);
			KX[lr]     = gridinfo[getindex(tempstr)].pointer;            /* hydraulic conductivity in x direction for layer lr */
			sprintf(tempstr,"%s%d",_gw_k_y_,lr);
			KY[lr]     = gridinfo[getindex(tempstr)].pointer;            /* hydraulic conductivity in y direction for layer lr */
			sprintf(tempstr,"%s%d",_gw_storage_coeff_,lr);
			S0[lr]     = gridinfo[getindex(tempstr)].pointer;  /* storage coefficients for layer lr */
			if (lr==1) s0index=getindex(tempstr);
			sprintf(tempstr,"%s%d",_gw_boundary_fix_h_,lr); 
			BOUNDH[lr] = gridinfo[getindex(tempstr)].pointer; /* Dirichlet-boundary conditions (h=fix) for layer lr */
			sprintf(tempstr,"%s%d",_gw_boundary_fix_q_,lr); 
			BOUNDQ[lr] = gridinfo[getindex(tempstr)].pointer; /* Neumann-boundary condition dh/dx or dh/dy = fix for layer lr */
		}
		for (int lr=1;lr<layers;lr++) /* getting kolmation only for layers-1 levels */
		{
			sprintf(tempstr,"%s%d",_gw_kolmation_,lr); KOLM[lr] = gridinfo[getindex(tempstr)].pointer; /* Neumann-boundary condition dh/dx or dh/dy = fix for layer 1 */
		}

		// now setting the first aquifers thickness to the thickness of the soil including the new scaling factor
		double **SOIL = gridinfo[getindex(_Bodenarten_)].pointer;
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
				{
					if ((THICK[1][k][j] <= nodata) || (SOIL[k][j] <= nodata))
					{
						THICK[1][k][j] = nodata;
					}
					else
					{
						THICK[1][k][j] = this->SOILtable[int(SOIL[k][j])].getTotalThickness();
					}
				}
			}
		}
		MPISynchronizeGridHalo(0, THICK[1]);
		//MPISynchronizeGridHaloJustRecv(0, THICK[1]); //LAURA/JOERG
		//MPISynchronizeGridHaloJustSend(0, THICK[1]); //LAURA/JOERG
		//MPISynchronizeGridHaloJustWait(0, THICK[1]); //LAURA/JOERG



		// searching for the grid containing the thawing depth for permafrost soils.
		// if found, permafrost must be regarded at points where THAWDEPTH hat valid values
		if (ist_vorhanden(_PermafrostThawDepth_)&&ist_vorhanden(_SnowFreeDaysGrid_)&&ist_vorhanden(_ThawCoeffPermaFrost_))
		{
			thawdepthindex = getindex(_PermafrostThawDepth_);
			THAWDEPTH=gridinfo[thawdepthindex].pointer;
		}
		else
		{
			if (ist_vorhanden(_ThawDepthTMod_))
			{
				thawdepthindex = getindex(_ThawDepthTMod_);
				THAWDEPTH=gridinfo[thawdepthindex].pointer;
				frozendepthindex = getindex(_frozenSoilInGroundwater_);
				FROZENDEPTH=gridinfo[frozendepthindex].pointer;
			}
			else
			{
				thawdepthindex = -1;
				THAWDEPTH=NULL;
				frozendepthindex=-1;
				FROZENDEPTH=NULL;
			}
		}

		// if silting up model is switched off, then surface storage from surface routing (which has the same name as surface storage from silting up) should be considered for hydraulic heads when the soil is saturated
		// Otherwise, set QDSU to NULL to indicate that no additional
		// input from silting up nor surface routing has to be regarded in unsatzon model
		if ((ist_vorhanden(_SURFSTORSILTINGUP_)>0) && (ist_vorhanden(_DIRECTRUNOFFSILTINGUP_)<=0))
		{
			QDSU = gridinfo[getindex(_SURFSTORSILTINGUP_)].pointer;
//			QDSU=NULL;
		}
		else
		{
			QDSU=NULL;
		}

		/* where is memory space for new grids in gridinfo? looking for the first free pointer */
		int gridindex=1;
		while (gridinfo[gridindex].identification!=NULL)
		{
			gridindex++;
			if (gridindex==gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no free grid (gridmaxanzahl must be larger than "<<gridmaxanzahl<<")!\n"<< "[end-entry]";
				throw(12);
			}
		}

		/* reading stored hydraulic heads or generating it new */
		if (gridslesen!=1)
		{
			/* grids are not read from disk */
			for (int lr=1;lr<=layers;lr++)
			{
				HEAD[lr]  = matrix(sr,er,1,ncols,nrows);
				sprintf(tempstr,"%s%d",_gw_head_,lr);
				initialize_grid(gridindex, HEADgridname[lr], HEAD[lr], NULL, 0, tempstr, HEADwritegrid[lr], sr, er);
				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int j=FirstLastCol[k].First;j<=FirstLastCol[k].Last;j++) 
					{
						if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
						{
							double icevalue = ICE ? ICE[k][j] : 0; // valid range of icevalue now between 0...2
							if (icevalue > 1.0) icevalue = 1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
							else if (icevalue < 0.0) icevalue = 0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
							if (iPondAvailable)
							{
								if (fabs(GWSTAND[k][j]) < 0.001)
								{
									// when the groundwater stands at the soil surface, then the POND content counts for HEAD-initialization, too
									HEAD[lr][k][j] = zn[k][j] == nodata ? nodata : DHM[k][j] + GWSTAND[k][j] + POND[k][j] / 1000.0;
								}
								else
								{
									// no saturation, so only the real GWSTAND value is needed for initializing HEAD
									HEAD[lr][k][j] = zn[k][j] == nodata ? nodata : DHM[k][j] + GWSTAND[k][j];
								}
							}
							else
							{
								if ((QDSU) && (QDSU[k][j] > 0) && (fabs(GWSTAND[k][j]) < 0.001))
								{
									HEAD[lr][k][j] = zn[k][j] == nodata ? nodata : DHM[k][j] + GWSTAND[k][j] + QDSU[k][j] / 1000.0;
								}
								else
								{
									HEAD[lr][k][j] = zn[k][j] == nodata ? nodata : DHM[k][j] + GWSTAND[k][j];
								}
							}
						}
					}
				}
				headindex[lr] = gridindex;
				MPISynchronizeGridHalo(headindex[lr], HEAD[lr]);
				//MPISynchronizeGridHaloJustRecv(headindex[lr], HEAD[lr]); //LAURA/JOERG
				//MPISynchronizeGridHaloJustSend(headindex[lr], HEAD[lr]); //LAURA/JOERG
				//MPISynchronizeGridHaloJustWait(headindex[lr], HEAD[lr]); //LAURA/JOERG
				gridindex++;

				// moved here from creating new grids
				HEADold[lr]  = matrix(sr,er,1,ncols,nrows);
				sprintf(tempstr,"%s%d",_gw_head_old_,lr);
				initialize_grid(gridindex, HEADoldgridname[lr], HEADold[lr], NULL, 0, tempstr, HEADoldwritegrid[lr], sr,er);
				for (int k = mpistartrow; k <= mpiendrow; k++)
				{
					for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
					{
						if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
						{
							HEADold[lr][k][j] = HEAD[lr][k][j];
						}
						else
						{
							HEADold[lr][k][j] = nodata;
						}
					}
				}
				headoldindex[lr] = gridindex;
				MPISynchronizeGridHalo(headoldindex[lr], HEADold[lr]);
				//MPISynchronizeGridHaloJustRecv(headoldindex[lr], HEADold[lr]); //LAURA
				//MPISynchronizeGridHaloJustSend(headoldindex[lr], HEADold[lr]); //LAURA
				//MPISynchronizeGridHaloJustWait(headoldindex[lr], HEADold[lr]); //LAURA
				gridindex++;

/*				HEADoldold[lr]  = matrix(1,nrows,1,ncols);
				sprintf(tempstr,"%s%d",_gw_head_old_"old_",lr);
				initialize_grid(gridindex,HEADoldoldgridname[lr],HEADoldold[lr],NULL,0,tempstr,HEADoldoldwritegrid[lr]);
				for (int k=1;k<=nrows;k++) 
				{
					for (int j=1;j<=ncols;j++) 
					{
						HEADoldold[lr][k][j]=zn[k][j]==nodata ? nodata : HEAD[lr][k][j];
					}
				}
				headoldoldindex[lr] = gridindex;
				gridindex++;

				HEADolditer[lr]  = matrix(1,nrows,1,ncols);
				sprintf(tempstr,"%s%d",_gw_head_old_"iter_",lr);
				initialize_grid(gridindex,HEADolditergridname[lr],HEADolditer[lr],NULL,0,tempstr,HEADolditerwritegrid[lr]);
				for (int k=1;k<=nrows;k++) 
				{
					for (int j=1;j<=ncols;j++) 
					{
						HEADolditer[lr][k][j]=zn[k][j]==nodata ? nodata : HEAD[lr][k][j];
					}
				}
				headolditerindex[lr] = gridindex;
				gridindex++;
*/
			}
		}
		else
		{
			/* alte Dateien einlesen */
			for (int lr=1;lr<=layers;lr++)
			{
				readgridfromfile(HEADgridname[lr], &HEAD[lr], 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOn_);
				sprintf(tempstr,"%s%d",_gw_head_,lr);
				initialize_grid(gridindex, HEADgridname[lr], HEAD[lr], NULL, 0, tempstr, HEADwritegrid[lr], sr,er);
				headindex[lr] = gridindex;
				gridindex++;

				readgridfromfile(HEADoldgridname[lr], &HEADold[lr], 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOn_);
				sprintf(tempstr,"%s%d",_gw_head_old_,lr);
				initialize_grid(gridindex, HEADoldgridname[lr], HEADold[lr], NULL, 0, tempstr, HEADoldwritegrid[lr], sr,er);
				headoldindex[lr] = gridindex;
				gridindex++;

/*				
				readgridfromfile(HEADoldoldgridname[lr],&HEADoldold[lr],1);
				sprintf(tempstr,"%s%d",_gw_head_old_"old_",lr);
				initialize_grid(gridindex,HEADoldoldgridname[lr],HEADoldold[lr],NULL,0,tempstr,HEADoldoldwritegrid[lr]);
				headoldoldindex[lr] = gridindex;
				gridindex++;

				readgridfromfile(HEADolditergridname[lr],&HEADolditer[lr],1);
				sprintf(tempstr,"%s%d",_gw_head_old_"iter_",lr);
				initialize_grid(gridindex,HEADolditergridname[lr],HEADolditer[lr],NULL,0,tempstr,HEADolditerwritegrid[lr]);
				headolditerindex[lr] = gridindex;
				gridindex++;
*/
			}
		}

		/* now generating new grids */

		for (int lr=1;lr<=layers;lr++)
		{
/*			HEADold[lr]  = matrix(1,nrows,1,ncols);
			sprintf(tempstr,"%s%d",_gw_head_old_,lr);
			initialize_grid(gridindex,tempstr,HEADold[lr],NULL,0,tempstr,0);
			for (int k=1;k<=nrows;k++) 
			{
				for (int j=1;j<=ncols;j++) 
				{
					HEADold[lr][k][j]=zn[k][j]==nodata ? nodata : HEAD[lr][k][j];
				}
			}
			headoldindex[lr] = gridindex;
			gridindex++;
*/
			HEADoldold[lr]  = matrix(sr,er,1,ncols,nrows);
			sprintf(tempstr,"%s%d",_gw_head_old_"old_",lr);
			initialize_grid(gridindex, tempstr, HEADold[lr], NULL, 0, tempstr, 0, sr,er);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						HEADoldold[lr][k][j] = HEAD[lr][k][j];
					}
					else
					{
						HEADoldold[lr][k][j] = nodata;
					}
				}
			}
			headoldoldindex[lr] = gridindex;
			MPISynchronizeGridHalo(headoldoldindex[lr], HEADoldold[lr]);
			//MPISynchronizeGridHaloJustRecv(headoldoldindex[lr], HEADoldold[lr]); //LAURA
			//MPISynchronizeGridHaloJustSend(headoldoldindex[lr], HEADoldold[lr]); //LAURA
			//MPISynchronizeGridHaloJustWait(headoldoldindex[lr], HEADoldold[lr]); //LAURA
			gridindex++;

			HEADolditer[lr] = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			sprintf(tempstr,"%s%d",_gw_head_old_"old_",lr);
			initialize_grid(gridindex, tempstr, HEADolditer[lr], NULL, 0, tempstr, 0, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						HEADolditer[lr][k][j] = HEAD[lr][k][j];
					}
					else
					{
						HEADolditer[lr][k][j] = nodata;
					}
				}
			}
			headolditerindex[lr] = gridindex;
//			MPISynchronizeGridHalo(headolditerindex[lr], HEADolditer[lr]);
			gridindex++;

			TRANSX[lr] = matrix(sr, er, 1, ncols, nrows);
			sprintf(tempstr,"%s%d",_gw_transmissivity_x_,lr);
			initialize_grid(gridindex, tempstr, TRANSX[lr], NULL, 0, tempstr, 0, sr, er);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						TRANSX[lr][k][j] = KX[lr][k][j] * (THICK[lr][k][j] + GWSTAND[k][j]);
					}
					else
					{
						TRANSX[lr][k][j] = nodata;
					}
				}
			}
			transxindex[lr] = gridindex;
			MPISynchronizeGridHalo(transxindex[lr], TRANSX[lr]);
			//MPISynchronizeGridHaloJustRecv(transxindex[lr], TRANSX[lr]); //LAURA
			//MPISynchronizeGridHaloJustSend(transxindex[lr], TRANSX[lr]); //LAURA
			//MPISynchronizeGridHaloJustWait(transxindex[lr], TRANSX[lr]); //LAURA
			gridindex++;

			TRANSY[lr] = matrix(sr, er, 1, ncols, nrows);
			sprintf(tempstr,"%s%d",_gw_transmissivity_y_,lr);
			initialize_grid(gridindex, tempstr, TRANSY[lr], NULL, 0, tempstr, 0, sr, er);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						TRANSY[lr][k][j] = KY[lr][k][j] * (THICK[lr][k][j] + GWSTAND[k][j]);
					}
					else
					{
						TRANSY[lr][k][j] = nodata;
					}
				}
			}
			transyindex[lr] = gridindex;
			MPISynchronizeGridHalo(transyindex[lr], TRANSY[lr]);
			//MPISynchronizeGridHaloJustRecv(transyindex[lr], TRANSY[lr]); //LAURA
			//MPISynchronizeGridHaloJustSend(transyindex[lr], TRANSY[lr]); //LAURA
			//MPISynchronizeGridHaloJustWait(transyindex[lr], TRANSY[lr]); //LAURA
			gridindex++;

			FLOWX[lr] = matrix(sr, er, 1, ncols, nrows);
			sprintf(tempstr,"%s%d",_gw_flow_x_,lr);
			initialize_grid(gridindex, FLOWXgridname[lr], FLOWX[lr], NULL, 0, tempstr, FLOWXwritegrid[lr], sr, er);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						FLOWX[lr][k][j] = 0;
					}
					else
					{
						FLOWX[lr][k][j] = nodata;
					}
				}
			}
			flowxindex[lr] = gridindex;
			MPISynchronizeGridHalo(flowxindex[lr], FLOWX[lr]);
			//MPISynchronizeGridHaloJustRecv(flowxindex[lr], FLOWX[lr]); //LAURA
			//MPISynchronizeGridHaloJustSend(flowxindex[lr], FLOWX[lr]); //LAURA
			//MPISynchronizeGridHaloJustWait(flowxindex[lr], FLOWX[lr]); //LAURA
			gridindex++;

			FLOWY[lr] = matrix(sr, er, 1, ncols, nrows);
			sprintf(tempstr,"%s%d",_gw_flow_y_,lr);
			initialize_grid(gridindex, FLOWYgridname[lr], FLOWY[lr], NULL, 0, tempstr, FLOWYwritegrid[lr], sr,er);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						FLOWY[lr][k][j] = 0;
					}
					else
					{
						FLOWY[lr][k][j] = nodata;
					}
				}
			}
			flowyindex[lr] = gridindex;
			MPISynchronizeGridHalo(flowyindex[lr], FLOWY[lr]);
			//MPISynchronizeGridHaloJustRecv(flowyindex[lr], FLOWY[lr]); //LAURA
			//MPISynchronizeGridHaloJustSend(flowyindex[lr], FLOWY[lr]); //LAURA
			//MPISynchronizeGridHaloJustWait(flowyindex[lr], FLOWY[lr]); //LAURA
			gridindex++;

			BALANCE[lr] = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			sprintf(tempstr,"%s%d",_gw_balance_,lr);
			initialize_grid(gridindex, BALANCEgridname[lr], BALANCE[lr], NULL, 0, tempstr, BALANCEwritegrid[lr], mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						BALANCE[lr][k][j] = 0;
					}
					else
					{
						BALANCE[lr][k][j] = nodata;
					}
				}
			}
			balanceindex[lr] = gridindex;
			gridindex++;

			LEAKAGE[lr] = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			sprintf(tempstr,"%s%d",_gw_leakage_,lr);
			initialize_grid(gridindex, tempstr, LEAKAGE[lr], NULL, 0, tempstr, 0, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						LEAKAGE[lr][k][j] = 0;
					}
					else
					{
						LEAKAGE[lr][k][j] = nodata;
					}
				}
			}
			leakageindex[lr] = gridindex;
			//MPISynchronizeGridHalo(leakageindex[lr], LEAKAGE[lr]);
			gridindex++;

			UPLIM[lr] = matrix(sr,er, 1, ncols, nrows);
			sprintf(tempstr,"%s%d",_aquifer_upper_limit_,lr);
			initialize_grid(gridindex, tempstr, UPLIM[lr], NULL, 0, tempstr, 0, sr,er);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					if ((zn[k][j] != nodata) && (MPIREG[k][j] == mpi_rank))
					{
						UPLIM[lr][k][j] = (lr == 1) ? DHM[k][j] : UPLIM[lr - 1][k][j] - THICK[lr - 1][k][j];
					}
					else
					{
						UPLIM[lr][k][j] = nodata;
					}
				}
			}
			uplimindex[lr] = gridindex;
			MPISynchronizeGridHalo(uplimindex[lr], UPLIM[lr]);
			//MPISynchronizeGridHaloJustRecv(uplimindex[lr], UPLIM[lr]); //LAURA
			//MPISynchronizeGridHaloJustSend(uplimindex[lr], UPLIM[lr]); //LAURA
			//MPISynchronizeGridHaloJustWait(uplimindex[lr], UPLIM[lr]); //LAURA
			gridindex++;
		}

		dx=dy=cellsize;

		if ((control_row>0)&&(MPIREG[control_row]!=NULL)&&(mpi_rank==MPIREG[control_row][control_col])) /* writing the files only, if control_row is greater than zero */
		{
			/* file buffer for control point fluxes */
			glogbufmax      = (unsigned long)((long)buffersize/(4*sizeof(long)+(6*layers+1)*sizeof(double)));
			if (glogbufmax<1) glogbufmax=1;
			glogbufcnt      = 0; /* initialising the buffer index counter with 0 */
			glogbuffer      = (buffertype*)calloc(glogbufmax+1,sizeof(buffertype));
			for (bufrow=1;bufrow<=glogbufmax;bufrow++) glogbuffer[bufrow].wert=(double*)calloc(6*layers+1,sizeof(double));

			if ((glogfile=fopen(glogfilename,"wt"))==NULL)
			{
				(*WSstream) << "FATAL: error opening file for gw-fluxes (single point) "
				"\n-> check disk space or user rights\n"<<"[end-entry]";
				throw(62);
			}
			fprintf(glogfile,"YY\tMM\tDD\tHH\tcontrol_point_fluxes[%d][%d]",control_row,control_col);
			for (int i = 1; i <= layers; i++)
			{
				fprintf(glogfile, "[lr: %d, 1]\t[lr: %d, 2]\t[lr: %d, 3]\t[lr: %d, 4]\t[lr: %d, up]\t[lr: %d, down]\t", i, i, i, i, i, i);
			}
			fprintf(glogfile,"\n");
			fclose(glogfile);
		}

		if (tempstr) free(tempstr); tempstr=NULL;
		initialized=0; // darf erst nach SetQgesADresse auf 1 gesetzt werden

		this->initialize_lakes(_lake_model_);
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
//#pragma optimize ("p",off)

void gwflowclass::setSubStepCnt(int step)
{
	substepcnt=step;
}

//#pragma optimize ("p",on)
void gwflowclass::SetRiverflow_tableAdress(riverflow_table_type *adresse)
{
	if (do_it==0) return; else riverflow_table=adresse;
}
//#pragma optimize ("p",off)


/* the irrigation demand has to be checked against the available water from surface water (rivers, reservoirs). 
So the adresses of the flow-arrays must be known within the irrigation moduler */
//#pragma optimize ("p",on)
void gwflowclass::SetQgesAdresse(double **p_Qges)
{
	if (do_it==0) 
	return;
	else  
	Qges=p_Qges;
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void gwflowclass::transmissivities()
{
	/* transmissivities for the first aquifer (other aquifers are assumed to be 
	confined and so they have conctant transmissivities which have already 
	been calculated during initialization) */

	for (int lr=1;lr<=layers;lr++)
	{
#pragma omp parallel for schedule(runtime) shared(lr)
		for (int i=mpistartrow;i<=mpiendrow;i++)
		{
			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				double T0,T1,TS;
				double icevalue2=0;
				double icevalue=ICE?ICE[i][j]:0; // valid range of icevalue now between 0...2
				if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
				else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
				if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank)&&(icevalue<1.0))
				{

					double	thick0; // thickness of actual cell
					double	thick1; // thickness of the cell east to the actual cell (j+1)
					double thick2; // thickness of the cell south to the actual cell (i+1)

					thick0=THICK[lr][i][j]; 
					if (j<ncols) 
					{
						icevalue2 = ICE?ICE[i][j+1]:0; // icevalue may be 0...1 for dynamic glaciers and 1 or 2 for static glaciers. In any case, values below 1 mean that there is at least a fraction of the target cell unglacierized, thus allowing some groundwater flow to the actual cell
						if (icevalue2>1.0) icevalue2=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
						else if (icevalue2<0.0) icevalue2=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
						if (icevalue2<1.0)
						{
							thick1=THICK[lr][i][j+1];
						}
						else
						{
							thick1=0;
						}
					} 
					else 
					{
						thick1=0;
						icevalue2=0;
					}
					
					if (i<nrows) 
					{
						icevalue2 = ICE?ICE[i+1][j]:0; // icevalue may be 0...1 for dynamic glaciers and 1 or 2 for static glaciers. In any case, values below 1 mean that there is at least a fraction of the target cell unglacierized, thus allowing some groundwater flow to the actual cell
						if (icevalue2>1.0) icevalue2=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
						else if (icevalue2<0.0) icevalue2=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
						if (icevalue2<1.0)
						{
							thick2=THICK[lr][i+1][j];
						}
						else
						{
							thick2=0;
						}
					} 
					else 
					{
						thick2=0;
						icevalue2=0;
					}


					// here, we will manipulate the thicknesses for the first aquifer when a THAWDEPTH-grid is available
					if ((FROZENDEPTH)&&(lr==1)) {
						if ((FROZENDEPTH[i][j]>nodata))
						{
							thick0 -= FROZENDEPTH[i][j]*frozenDepthScaling;
							if (thick0<0) thick0=0;
						}
						if ((j<ncols)&&(FROZENDEPTH[i][j+1]>nodata))
						{
							thick1 -= FROZENDEPTH[i][j+1]*frozenDepthScaling;
							if (thick1<0) thick1=0;
						}
						if ((i<nrows)&&(FROZENDEPTH[i+1][j]>nodata))
						{
							thick2 -= FROZENDEPTH[i+1][j]*frozenDepthScaling;
							if (thick2<0) thick2=0;
						}
					} 
					// old code. obsolete since froznedepth now takes care of active layer.
					//else if ((THAWDEPTH)&&(lr==1)) {
						// if ((THAWDEPTH[i][j]>nodata)&&(-THAWDEPTH[i][j]<thick0))
						// {
							// thick0=-THAWDEPTH[i][j];
						// }
						// if ((j<ncols)&&(THAWDEPTH[i][j+1]>nodata)&&(-THAWDEPTH[i][j+1]<thick1))
						// {
							// thick1=-THAWDEPTH[i][j+1];
						// }
						// if ((i<nrows)&&(THAWDEPTH[i+1][j]>nodata)&&(-THAWDEPTH[i+1][j]<thick2))
						// {
							// thick2=-THAWDEPTH[i+1][j];
						// }
					//}


					/* transmissivity in x-direction */
					if (HEAD[lr][i][j]>(UPLIM[lr][i][j]+0.01)) 
					// only if the hydraulic head is significantly higher than the upper limit of the aquifer 
					// --> is important for permafrost modelling and saturated aquifers to prevent confined conditions 
					// when numerical instabilities may cause slight overestimations of hydraulic head
					{
						T0  = THICK[lr][i][j]*KX[lr][i][j];  /* confined layer */
					}
					else
					{
						T0  = (HEAD[lr][i][j]-UPLIM[lr][i][j]+thick0)*KX[lr][i][j];
					}
					if (icevalue>0) T0=0;
					// different versions of reducing the flow from partly glacierized cells
					T0 *= ((1.0-icevalue)*(1.0-icevalue)*(1.0-icevalue)*(1.0-icevalue));

					if ((j<ncols)&&(zn[i][j+1]!=nodata) && !((this->do_it==3)&&(zn[i][j]!=zn[i][j+1])))
					{
						icevalue2 = ICE?ICE[i][j+1]:0; // icevalue may be 0...1 for dynamic glaciers and 1 or 2 for static glaciers. In any case, values below 1 mean that there is at least a fraction of the target cell unglacierized, thus allowing some groundwater flow to the actual cell
						if (icevalue2>1.0) icevalue2=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
						else if (icevalue2<0.0) icevalue2=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
						if (icevalue2<1.0)
						{
							if (HEAD[lr][i][j+1]>=(UPLIM[lr][i][j+1]+0.01)) 
							// only if the hydraulic head is significantly higher than the upper limit of the aquifer 
							// --> is important for permafrost modelling and saturated aquifers to prevent confined conditions 
							// when numerical instabilities may cause slight overestimations of hydraulic head
							{
								T1  = THICK[lr][i][j+1]*KX[lr][i][j+1];  /* confined layer */
							}
							else
							{
								T1  = (HEAD[lr][i][j+1]-UPLIM[lr][i][j+1]+thick1)*KX[lr][i][j+1];
							}
							if (icevalue2>0) T1=0;
							// different versions of reducing the flow from partly glacierized cells
							T1 *= ((1.0-icevalue2)*(1.0-icevalue2)*(1.0-icevalue2)*(1.0-icevalue2));
						}
						else
						{
							T1=0.0;
						}
					}
					else
					{
						T1=0.0;
						icevalue2=0;
					}

					/* if the cell is not a no-flow-boundary (marked by K==0),
			it should be possible for water to enter the cell even if the cell is dry.
			in such a case, the transmissivities are build only from the neigboring
			cell from which water com3es into the actual dry cell.
			This should be true also for groundwater tables below the bottom of an aquifer
			(then is T<0 because of THICK+GWSTAND <0) */
					// making sure that a dry cell can be filled up but an empty cell cannot drain more water: actual cell
					if (T0<=0.0) 
					{
						if ((T1>0)&&(HEAD[lr][i][j+1]>HEAD[lr][i][j])) // but only if the gradient points towards the actual cell
						{
							if ((KX[lr][i][j  ]>0.0)&&(T1>0.0)&&(icevalue==0.0))  // only if it is not glacierized
							{
								T0=T1; 
							} 
							else 
							{
								T0=0.0;
							}
						}
						else
						{
							T0=0.0;
						}
					}
					// making sure that a dry cell can be filled up but an empty cell cannot drain more water: eastern cell
					if ((T1<=0.0)&&(j<ncols)&&(zn[i][j+1]!=nodata) && !((this->do_it==3)&&(zn[i][j]!=zn[i][j+1])))
					{
						if ((T0>0)&&(HEAD[lr][i][j]>HEAD[lr][i][j+1]))
						{
							icevalue2 = ICE?ICE[i][j+1]:0; 
							if (icevalue2>1.0) icevalue2=1.0; 
							else if (icevalue2<0.0) icevalue2=0.0; 
							if ((KX[lr][i][j+1]>0.0)&&(T0>0.0)&&(icevalue2==0.0))  
							{
								T1=T0; 
							} 
							else 
							{
								T1=0.0;
							}
						}
						else
						{
							T1=0.0;
						}
					}

					TS=T0+T1;
					TRANSX[lr][i][j]=TS>0.0?2*T0*T1/TS:0.0;

					/* transmissivity in y-direction */
					if (HEAD[lr][i][j]>=(UPLIM[lr][i][j]+0.01)) 
					// only if the hydraulic head is significantly higher than the upper limit of the aquifer 
					// --> is important for permafrost modelling and saturated aquifers to prevent confined conditions 
					// when numerical instabilities may cause slight overestimations of hydraulic head
					{ 
						T0  = THICK[lr][i][j]*KY[lr][i][j];  /* confined layer */
					} 
					else 
					{
						T0  = (HEAD[lr][i][j]-UPLIM[lr][i][j]+thick0)*KY[lr][i][j];
					}
					if (icevalue>0) T0=0;
					// different versions of reducing the flow from partly glacierized cells
					T0 *= ((1.0-icevalue)*(1.0-icevalue)*(1.0-icevalue)*(1.0-icevalue));

					if ((i<nrows)&&(zn[i+1][j]!=nodata) && !((this->do_it==3)&&(zn[i][j]!=zn[i+1][j])))
					{
						icevalue2 = ICE?ICE[i+1][j]:0; // icevalue may be 0...1 for dynamic glaciers and 1 or 2 for static glaciers. In any case, values below 1 mean that there is at least a fraction of the target cell unglacierized, thus allowing some groundwater flow to the actual cell
						if (icevalue2>1.0) icevalue2=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
						else if (icevalue2<0.0) icevalue2=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
						if (icevalue2<1.0)
						{
							if (HEAD[lr][i+1][j]>=(UPLIM[lr][i+1][j]+0.01)) 
							// only if the hydraulic head is significantly higher than the upper limit of the aquifer 
							// --> is important for permafrost modelling and saturated aquifers to prevent confined conditions 
							// when numerical instabilities may cause slight overestimations of hydraulic head
							{
								T1  = THICK[lr][i+1][j]*KY[lr][i+1][j];  /* confined layer */
							} 
							else 
							{
								T1  = (HEAD[lr][i+1][j]-UPLIM[lr][i+1][j]+thick2)*KY[lr][i+1][j];
							}
							if (icevalue2>0) T1=0;
							// different versions of reducing the flow from partly glacierized cells
							T1 *= ((1.0-icevalue2)*(1.0-icevalue2)*(1.0-icevalue2)*(1.0-icevalue2));
						}
						else
						{
							T1 = 0.0;
						}
					}
					else
					{
						T1=0.0;
						icevalue2=0;
					}
					/* if the cell is not a no-flow-boundary (marked by K==0),
			it should be possible for water to enter the cell even if the cell is dry.
			in such a case, the transmissivities are biuld only from the neigboring
			cell from which water com3es into the actual dry cell.
			This should be true also for groundwater table below the bottom of aquifer
			(then is T<0 because of THICK+GWSTAND <0) */
					// making sure that a dry cell can be filled up but an empty cell cannot drain more water: actual cell
					if ((T0<=0.0))
					{
						if ((T1>0)&&(HEAD[lr][i+1][j]>HEAD[lr][i][j]))
						{
							if ((KY[lr][i  ][j]>0.0)&&(T1>0.0)&&(icevalue==0.0)) 
							{
								T0=T1; 
							}
							else
							{
								T0=0.0;
							}
						}
						else
						{
							T0=0.0;
						}
					}
					// making sure that a dry cell can be filled up but an empty cell cannot drain more water: southern cell
					if ((T1<=0.0)&&(i<nrows)&&(zn[i+1][j]!=nodata) && !((this->do_it==3)&&(zn[i][j]!=zn[i+1][j])))
					{
						if ((T0>0)&&(HEAD[lr][i][j]>HEAD[lr][i+1][j]))
						{
							icevalue2 = ICE?ICE[i+1][j]:0; 
							if (icevalue2>1.0) icevalue2=1.0; 
							else if (icevalue2<0.0) icevalue2=0.0; 
							if ((KY[lr][i+1][j]>0.0)&&(T0>0.0)&&(icevalue2==0.0))
							{
								T1=T0; 
							}
							else
							{
								T1=0.0;
							}
						}
						else
						{
							T1=0.0;
						}
					}
					TS=T0+T1;
					TRANSY[lr][i][j]=TS>0.0?2*T0*T1/TS:0.0;
				}
				else
				{
					// glacierized cells should be handled as having frozen ground, so the transmissivity is 0, as should be set for nodata cells in zn, since the algorithm is taking into account cells from the nodata neighborhood as well
					// do NOT call this line for MPIREG[i][j]==mpi_rank only, since in gauss_seidel_iteration the values are assumed to have valid values even for cells outside the actual zones (shoudl then be 0)
					TRANSX[lr][i][j]=TRANSY[lr][i][j]=0.0;
				}
			}
		}
	}
}
//#pragma optimize ("p",off)


double gwflowclass::HorizontalEquationSystem(int lr)
{
	double DT = 60.0*InV / substepcnt; /* timestep in seconds */
	double maxerror = 0; /* initializing maxerror with 0 */
#ifdef _OPENMP
#ifndef _OMP3
	int nthreads = omp_get_max_threads();
	double *maxerrorarray = new double[nthreads];
	for (int n = 0; n < nthreads; n++) maxerrorarray[n] = 0;

#pragma omp parallel for schedule(runtime) shared(lr,DT,maxerrorarray) 
#else
#pragma omp parallel for schedule(runtime) shared(lr,DT) reduction(max:maxerror)
#endif
#endif
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
#ifdef _OPENMP
		int th_id = omp_get_thread_num();
#endif
		double a, b, c, d, e, r; // temprary variables
		double err, headold, headoldold;

		double *coeff_a = new double[ncols];
		double *coeff_b = new double[ncols];
		double *coeff_c = new double[ncols];
		double *right_v = new double[ncols];
		double *result_x = new double[ncols];

		int start = FirstLastCol[i].First; // initializing cell counter, condition e) is initialized here
		int end = start - 1;
		int act_zn = nodata;
		double icevalue; 

		// the method applied here solves the equation system for rows only (head for time n+1), taking the heads of neighboring rows at time n, 
		// so they are known entities and put on the rght hand side of the equation (in right_v vector). By iteration or combining with the Gauss-Seidel-Algorithm the soution should converge to the same result as a real 2D solver
		// this technique is used for better applicabilty in MPI environments
		// continous cell regions must all fullfill the following conditions:
		// a) share the same MPI region
		// b) no nodata value in the zone grid
		// c) no glacier cell, which has a frozen ground and as such also no flow
		// d) share the same sub basion code (if do_it==3, otherwise this parameter is not regarded)
		// e) between FirstLastCol[i].First and FirstLastCol[i].Last
		// f) no constant head boundary
		// g) no no-flow boundary (i.e. T=0 bcause k=0 or h=0, i.e. like a threshold the water cannot pass yet)
		
		while (start <= FirstLastCol[i].Last) // checking condition e). There may be multiple sections to be solved duw to interrupts by nodata cells or boundary condition cells, so a while loope is needed here es encompassing element
		{
			// checking conditions a), b), c), d), e)
			// looking for first variable head cell in the region of interest to construct the equation system that will be solved with tri-diagonal matrix solver for each continuous cell range
			icevalue = ICE ? ICE[i][start] : 0; // valid range of icevalue now between 0...2
			if (icevalue > 1.0) icevalue = 1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
			else if (icevalue < 0.0) icevalue = 0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
			act_zn = zn[i][start];

			while (
				(start < FirstLastCol[i].Last) &&
				(
				(MPIREG[i][start] != mpi_rank) ||
				(zn[i][start] == nodata) ||
				(icevalue == 1.0) ||
				((act_zn != nodata) && (do_it == 3) && (zn[i][start] != act_zn)) ||
				(BOUNDH[lr][i][start] != nodata) ||
				((TRANSX[lr][i][start] <= 0) &&
				(TRANSY[lr][i][start] <= 0) &&
				((start>1) ? (TRANSX[lr][i][start - 1] <= 0) : true) &&
				((i>1) ? (TRANSY[lr][i - 1][start] <= 0) : true)
				)
				)
				)
			{
				// checking condition f): if conditions a to 3 are met, we need to make sure that the actual cell is no constant head boundary. 
				// If so, set the result to the boundary condition and increment the first-cell counter
				if (BOUNDH[lr][i][start] != nodata)
				{
					HEAD[lr][i][start] = BOUNDH[lr][i][start];
				}
				if ((BOUNDQ[lr][i][start] != nodata) &&
					(
					(TRANSX[lr][i][start] <= 0) &&
					(TRANSY[lr][i][start] <= 0) &&
					((start>1) ? (TRANSX[lr][i][start - 1] <= 0) : true) &&
					((i > 1) ? (TRANSY[lr][i - 1][start] <= 0) : true)
					)
					)
				{	// noflow cell that can still be recharged or drained be constant flux boundary condition
					HEAD[lr][i][start] += BOUNDQ[lr][i][start] * DT / S0[lr][i][start];
				}
				start++;
				act_zn = zn[i][start];
				icevalue = ICE ? ICE[i][start] : 0; // valid range of icevalue now between 0...2
				if (icevalue > 1.0) icevalue = 1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
				else if (icevalue < 0.0) icevalue = 0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells

			}
			end = start;
			if (end < FirstLastCol[i].Last) // at the end of a i we need to skip the solver
			{
				act_zn = zn[i][start];
				// otherwise search the next boundary condition or nodata cell etc (conditions a) to g) above) to get the last cell that can be used for the tridiag-solver
				// checking conditions a), b), c), d), e), f) and g) that will all break the actual continuous series of depending cells
				icevalue = ICE ? ICE[i][end+1] : 0; // valid range of icevalue now between 0...2
				if (icevalue > 1.0) icevalue = 1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
				else if (icevalue < 0.0) icevalue = 0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
				while ((end < FirstLastCol[i].Last) &&   // must be within allowed range of cells for this i
						(BOUNDH[lr][i  ][end + 1] == nodata) && // must not be a boundary condition of fixed head
						((TRANSX[lr][i][end + 1] > 0) || (TRANSY[lr][i][end + 1] > 0) || 
						 (TRANSX[lr][i  ][end    ] > 0) || ((i>1)?(TRANSY[lr][i-1][end + 1] > 0):true)) &&		// transmissivity in either direction x or y musn't be 0
						(MPIREG[i][end + 1] == mpi_rank) &&	// cell must belong to the actual mpi-rank --> halo cells are only regarded as neigboring cells for getting the head values from.
						(zn[i][end+1] != nodata) &&			// cell must be an active cell
						(icevalue < 1.0) &&					// musn't be a glacier cell (with 100% glacier on it)
						(((do_it == 3) && (zn[i][end+1] == act_zn)) || (do_it != 3)) // must belong to the same sub basin when do_it==3
					  )
				{
					end++;
					icevalue = ICE ? ICE[i][end + 1] : 0; // valid range of icevalue now between 0...2
					if (icevalue > 1.0) icevalue = 1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
					else if (icevalue < 0.0) icevalue = 0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
				}
			}
			// now running the solution algorithm
			int j = start;
			while (j <= end)
			{
				// compute coefficients
				// checking for no-flow boundary contiions (either at the edge of the grid or the neighbor is in another zone or in no zone or a completely glacierized cell


				// a, i.e. to the left of the actual cell
				if (j == 1)	{
					a = 0;
				} else if (((zn[i][j - 1] != zn[i][j]) && (do_it == 3)) || (zn[i][j - 1]==nodata) || (ICE?(ICE[i][j-1]==1.0):false)) {
					a = 0;
				} else {
					a = TRANSX[lr][i][j - 1] / CellSizeX[i][j - 1];
				}


				// c, i.e. to the right od the actual cell
				if (j == ncols)	{
					c = 0;
				}
				else if (((zn[i][j + 1] != zn[i][j]) && (do_it == 3)) || (zn[i][j + 1] == nodata) || (ICE ? (ICE[i][j + 1] == 1.0) : false)) {
					c = 0;
				}
				else {
					// TRANSX and TRANSY are already computed as transmissivities between two cells, always assuming the flow between tha actual cell and the cell to the east or to the south, so index j is the flow to the 
					// cell j+1 and index i is the flow to the cell i+1, while the flow j-1 is the flow from the cell to the west into the actual cell resp. i-1 is the flow from the cell to the north into the actual cell
					c = TRANSX[lr][i][j] / CellSizeX[i][j + 1]; 
				}


				// d (upper cell, heads are used from the actual time step, so the coefficient moves to the right hand side)
				if (i == 1) {
					d = 0;
				} else if (((zn[i-1][j] != zn[i][j]) && (do_it == 3)) || (zn[i-1][j] == nodata) || (ICE ? (ICE[i-1][j] == 1.0) : false)) {
					d = 0;
				} else {
					d = TRANSY[lr][i - 1][j] / CellSizeY[i][j] * (HEADold[lr][i][j]-HEADold[lr][i-1][j]) / CellSizeY[i - 1][j];
				}


				// e (lower cell, heads are used from the actual time step, so the coefficient moves to the right hand side)
				if (i == nrows) {
					e = 0;
				}
				else if (((zn[i + 1][j] != zn[i][j]) && (do_it == 3)) || (zn[i + 1][j] == nodata) || (ICE ? (ICE[i + 1][j] == 1.0) : false)) {
					e = 0;
				}
				else {
					e = TRANSY[lr][i][j] / CellSizeY[i][j] * (HEADold[lr][i+1][j] - HEADold[lr][i][j]) / CellSizeY[i + 1][j];
				}


				b = -CellSizeX[i][j] * S0[lr][i][j] / DT - a - c;

				if (BOUNDQ[lr][i][j] != nodata)
				{
					r = (-S0[lr][i][j] * HEAD[lr][i][j] / DT - e + d - BOUNDQ[lr][i][j])*CellSizeX[i][j];
				}
				else
				{
					r = (-S0[lr][i][j] * HEAD[lr][i][j] / DT - e + d)*CellSizeX[i][j];
				}
				
				// constant head boundary conditions
				if ((j > 1) && ((BOUNDH[lr][i][j - 1] != nodata) || (MPIREG[i][j- 1] != mpi_rank)))
				{
					// move term with a to the right hand side and set a to 0
					r = r - a*HEAD[lr][i][j - 1];
					a = 0;
				}
				if ((j < ncols) && ((BOUNDH[lr][i][j + 1] != nodata) || (MPIREG[i][j + 1] != mpi_rank)))
				{
					// move term with c to the right hand side and set a to 0
					r = r - c*HEAD[lr][i][j + 1];
					c = 0;
				}
				// put parameters into the vectors
				coeff_a[j - start] = a;
				coeff_b[j - start] = b;
				coeff_c[j - start] = c;
				right_v[j - start] = r;
				// move to the next cell to the right 
				j++;
			}

			// solving the tri-diagonal matrices
			int cnt = end - start; // if end==start then count is 0 and this is still valid
			for (int j = 1; j <= cnt; j++)
			{
				double m = coeff_a[j] / coeff_b[j - 1];
				coeff_b[j] -= (m*coeff_c[j - 1]);
				right_v[j] -= (m*right_v[j - 1]);
			}

			result_x[cnt] = right_v[cnt] / coeff_b[cnt];

			for (int j = cnt - 1; j>=0; j--)
			{
				result_x[j] = (right_v[j] - coeff_c[j] * result_x[j + 1]) / coeff_b[j];
			}

			for (int j = start; j <= end; j++)
			{
				headold = HEAD[lr][i][j]; // old head is used to calculate changes (errors) between iteration steps
				HEAD[lr][i][j] = result_x[j - start];
				err = fabs(HEAD[lr][i][j] - headold); 
#ifndef _OMP3
				//#pragma omp critical
#ifdef _OPENMP
				if (err > maxerrorarray[th_id])
				{
					maxerrorarray[th_id] = err;
				}
#else
				if (err > maxerror)
				{
					maxerror = err;
				}
#endif
#else
				if (err > maxerror)
				{
					maxerror = err;
				}
#endif
			}

			// increment start to the next cell 
			start = end + 1;
		}
		delete[] coeff_a;
		delete[] coeff_b;
		delete[] coeff_c;
		delete[] right_v;
		delete[] result_x;
	}
#ifdef _OPENMP
#ifndef _OMP3
	for (int n = 0; n < nthreads; n++)
	{
		if (maxerrorarray[n] > maxerror)
		{
			maxerror = maxerrorarray[n];
		}
	}
	delete[] maxerrorarray;
#endif
#endif

	return maxerror;
}

//#pragma optimize ("p",on)
double gwflowclass::gauss_seidel_iteration(int lr)
/* performing the iteration for the aquifer given by the parameter lr 
using the Gauss-Seidel-Algorithm with a Crank-Nicholson-Approximation */
{
	double sor_sum_local,sor_anz_local,DT;
	
	sor_sum_local=0;
	sor_anz_local=0;
	DT=60.0*InV/substepcnt; /* timestep in seconds */
	double maxerror=0; /* initializing maxerror with 0 */
#ifdef _OPENMP
#ifndef _OMP3
	int nthreads = omp_get_max_threads();
	double *maxerrorarray = new double[nthreads];
	for (int n = 0; n < nthreads; n++) maxerrorarray[n] = 0;

#pragma omp parallel for schedule(runtime) shared(lr,DT,maxerrorarray) reduction(+:sor_sum_local,sor_anz_local)
#else
#pragma omp parallel for schedule(runtime) shared(lr,DT) reduction(+:sor_sum_local,sor_anz_local), reduction(max:maxerror)
#endif
#endif
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
#ifdef _OPENMP
		int th_id = omp_get_thread_num();
#endif
		double err, headold, headoldold, const_boundary_q; /* for calculating iteration error */
		double h1,h2,h3,h4, /* heads  (1) back, (2) right, (3) front and (4) left cell if looking on the cell from a higer row number */
		T1,T2,T3,T4, /* Transmissivities for the flow from (1) back, (2) right, (3) front and (4) left cell to the centered cell */
		f1,f2,f3,f4; /* flows from (1) back to center, (2) center to right, (3) center to front, (4) left to center*/
		double fsum, Tsum;
		double atemp,btemp,ctemp,storage_coefficient,sor_local,dx,dy;
		for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
		{
			double icevalue = ICE ? ICE[i][j] : 0; // valid range of icevalue now between 0...2
			if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
			else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
			if ((zn[i][j]>nodata)&&(MPIREG[i][j]==mpi_rank)&&(icevalue<1.0)) // glacierized cells should not be taken into account (even partly glacierized)
			{
				if (BOUNDH[lr][i][j]==nodata)
				{
					if (RelCellArea!=NULL)
					{
						dx=CellSizeX[i][j];
						dy=CellSizeY[i][j];
					}
					else
					{
						dx=dy=cellsize;
					}
					headold    = HEAD[lr][i][j]; // old head is used to calculate changes (errors) between iteration steps 
					// heads for t' in centerd coordinates 
					// note that h1 and h4 are newly calculated values, already, if the row index is 2 or larger. This may be a problem when using many cores in openMP or when using many ranks in MPI
					// since a growing number of rows has then no precalculated values as starting points and get more and more instable
					// depending on alpha, it is an explicite (alpha = 0), an fully implicite (alpha = 1) or an inbetween method (0 < alpha < 1) 
					h1=i>1?     (1-alpha)*HEADold[lr][i-1][j]+alpha*HEAD[lr][i-1][j]:HEAD[lr][i][j];
					h2=j<ncols? (1-alpha)*HEADold[lr][i][j+1]+alpha*HEAD[lr][i][j+1]:HEAD[lr][i][j];
					h3=i<nrows? (1-alpha)*HEADold[lr][i+1][j]+alpha*HEAD[lr][i+1][j]:HEAD[lr][i][j];
					h4=j>1?     (1-alpha)*HEADold[lr][i][j-1]+alpha*HEAD[lr][i][j-1]:HEAD[lr][i][j];
					/* transmissivities in centered coordinates */
					T1=i>1?      dx/dy*TRANSY[lr][i-1][j]:0;
					T2=j<=ncols? dy/dx*TRANSX[lr][i][j]  :0;
					T3=i<=nrows? dx/dy*TRANSY[lr][i][j]  :0;
					T4=j>1?      dy/dx*TRANSX[lr][i][j-1]:0;
					/* lokal fluxes in centered coordinates */
					f1=T1*h1;
					f2=T2*h2;
					f3=T3*h3;
					f4=T4*h4;
					Tsum=T1+T3+T2+T4;
					fsum=f1+f2+f3+f4;
					/* if any layer is unconfined, the storage coefficient of the
						uppermost aquifer is used which is assumed to be the one for unconfined conditions */
					storage_coefficient=HEADold[lr][i][j]>=UPLIM[lr][i][j]?S0[lr][i][j]:S0[1][i][j];
					//          storage_coefficient=S0[lr][i][j];
					const_boundary_q=BOUNDQ[lr][i][j]>nodata?BOUNDQ[lr][i][j]:0.0;  /* constant in/outflows */
					if (lr==1)
						atemp=storage_coefficient/DT*HEADold[lr][i][j]+LEAKAGE[lr][i][j]+const_boundary_q;
					else
						atemp=storage_coefficient/DT*HEADold[lr][i][j]+LEAKAGE[lr][i][j]-LEAKAGE[lr-1][i][j]+const_boundary_q;
					atemp*=(dx*dy);
					btemp=fsum-(1-alpha)*HEADold[lr][i][j]*Tsum;
					ctemp=alpha*Tsum+dx*dy*storage_coefficient/DT;

					HEAD[lr][i][j]=(atemp+btemp)/ctemp;
					err=fabs(HEAD[lr][i][j]-headold);
					/* due to the much smaller storage coefficient in confined conditions,
			it is possible to use larger errors (since the same error of heads between
			succesive iteration steps means much less difference in flux-balance than
			for unconfined conditions) */
					if (HEADold[lr][i][j]>=UPLIM[lr][i][j]) err/=5.0;
#ifndef _OMP3
//#pragma omp critical
#ifdef _OPENMP
					if (err>maxerrorarray[th_id])
					{
						maxerrorarray[th_id] = err;
					}
#else
					if (err > maxerror)
					{
						maxerror = err;
					}
#endif
#else
					if (err>maxerror)
					{
						maxerror = err;
					}
#endif
					/* now applying successive over relaxation SOR */
					if (auto_sor)
					{
						headoldold = HEADolditer[lr][i][j]; /* forelast old head is used to calculate sor if auto_sor==1 */
						if (headoldold!=headold)
						sor_local=1+(HEAD[lr][i][j]-headold)/(headold-headoldold);
						else sor_local=1.0;
						if (sor_local<=0) sor_local=1.0; else if (sor_local>2) sor_local=2.0;
						HEAD[lr][i][j]=headold+sor_local*(HEAD[lr][i][j]-headold);
						{
							sor_sum_local+=sor_local;
							sor_anz_local++;
						}
						HEADolditer[lr][i][j]=headold;
					}
					else
					if (sor!=1.0)
					HEAD[lr][i][j]=headold+sor*(HEAD[lr][i][j]-headold);
				}
				else
				{
					// if a constant head consition is true for this cell, then all HEADs should be set to this constant head
					// in order to avoid oscillation indications
					if (BOUNDH[lr][i][j]!=nodata) {
						HEAD[lr][i][j]       = BOUNDH[lr][i][j];
						HEADold[lr][i][j]    = BOUNDH[lr][i][j];
						HEADoldold[lr][i][j] = BOUNDH[lr][i][j];
					}
				}
			}
		}
	}
	sor_sum=sor_sum_local;
	sor_anz=sor_anz_local;

#ifdef _OPENMP
#ifndef _OMP3
	for (int n = 0; n < nthreads; n++)
	{
		if (maxerrorarray[n] > maxerror)
		{
			maxerror = maxerrorarray[n];
		}
	}
	delete[] maxerrorarray;
#endif
#endif


	return maxerror;
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void gwflowclass::fluxfield()
/* calculates the velocities of flow in x- and y-direction in m/s */
{
	for (int lr=1;lr<=layers;lr++)
	{
#pragma omp parallel for schedule(runtime) shared(lr)
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				double dx, dy;
				double icevalue=ICE?ICE[i][j]:0; // valid range of icevalue now between 0...2
				if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
				else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells

				if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank)&&(icevalue<1.0))
				{
					if (RelCellArea!=NULL)
					{
						dx=CellSizeX[i][j];
						dy=CellSizeY[i][j];
					}
					else
					{
						dx=dy=cellsize;
					}
					/* Flux in x-direction in m/s */
					if (j<ncols)
					{
						double icevalue2 = ICE?ICE[i][j+1]:0; // icevalue may be 0...1 for dynamic glaciers and 1 or 2 for static glaciers. In any case, values below 1 mean that there is at least a fraction of the target cell unglacierized, thus allowing some groundwater flow to the actual cell
						if (icevalue2>1.0) icevalue2=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
						else if (icevalue2<0.0) icevalue2=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
						if ((zn[i][j]>nodata)&&(zn[i][j+1]>nodata))
						{
							if (icevalue2<1.0)
							{
								FLOWX[lr][i][j]=KX[lr][i][j]*(HEAD[lr][i][j]-HEAD[lr][i][j+1])/dx;
							}
							else
							{
								FLOWX[lr][i][j]=0;
							}
						}
						else
						{
							FLOWX[lr][i][j]=nodata;
						}
					}
					else
					{
						if (zn[i][j]>nodata)
						{
							FLOWX[lr][i][j]=0;
						}
						else
						{
							FLOWX[lr][i][j]=nodata;
						}
					}

					/* Flux in y-direction in m/s */
					if (i<nrows)
					{
						double icevalue2 = ICE?ICE[i+1][j]:0; // icevalue may be 0...1 for dynamic glaciers and 1 or 2 for static glaciers. In any case, values below 1 mean that there is at least a fraction of the target cell unglacierized, thus allowing some groundwater flow to the actual cell
						if (icevalue2>1.0) icevalue2=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
						else if (icevalue2<0.0) icevalue2=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
						if ((zn[i][j]>nodata)&&(zn[i+1][j]>nodata))
						{
							if (icevalue2<1.0)
							{
								FLOWY[lr][i][j]=KY[lr][i][j]*(HEAD[lr][i][j]-HEAD[lr][i+1][j])/dy;
							}
							else
							{
								FLOWY[lr][i][j]=0;
							}
						}
						else
						{
							FLOWY[lr][i][j]=nodata;
						}
					}
					else
					{
						if (zn[i][j]>nodata)
						{
							FLOWY[lr][i][j]=0;
						}
						else
						{
							FLOWY[lr][i][j]=nodata;
						}
					}
				}
				else
				{
					if (MPIREG[i][j]==mpi_rank)
					{
						FLOWX[lr][i][j]=0;
						FLOWY[lr][i][j]=0;
					}
				}
			}
		}
	}
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
double gwflowclass::get_oscillations(double **H,double **Hold,double **Holdold)
{
	int i,j;
	double k,diff1,diff2;
	k=0;
	sum_diff=0;
	maxsum_diff=0;
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
		{
			double icevalue = ICE ? ICE[i][j] : 0; // valid range of icevalue now between 0...2
			if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
			else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
			{
				// if a constant head boundary condition applies to this cell, then H, Hold and Holdold are set to the BOUNDH-value
				// resulting in diff1=diff2=0 with no oscillation at all.
				diff1=Hold[i][j]-Holdold[i][j];
				diff2=H[i][j]-Hold[i][j];
				if (((diff1>maxerrorf)&&(diff2<-maxerrorf))||((diff1<-maxerrorf)&&(diff2>maxerrorf)))
				{
					k++;
					sum_diff+=fabs(diff2);
					if (maxsum_diff<fabs(diff2)) maxsum_diff=fabs(diff2);
				}
			}
		}
	}
	return k;
}
//#pragma optimize ("p",off)

void gwflowclass::pre_run(long me,long se,long tg,long mt,long jr, long actsubstep)
{
	if (do_it == 0) return;
	/* Ausgabe einer Nachricht */
#ifdef _OPENMP
	if (!omp_in_parallel())
	{
		teilmodellklasse::run(me,se,tg,mt,jr,NULL,NULL);
	}
	else
	{
		minute = me;
		stunde = se;
		tag    = tg;
		monat  = mt;
		jahr   = jr;
	}
#else
	teilmodellklasse::run(me,se,tg,mt,jr,NULL,NULL);
#endif

}

//#pragma optimize ("p",on)
void gwflowclass::run(long me,long se,long tg,long mt,long jr, long actsubstep)
{
	double DT;

	if (do_it == 0) return;

	DT=60.0*InV/substepcnt; /* timestep in seconds */
	for (int z=1;z<=tracer_count;z++)
	{
		for (int lr=1;lr<=layers;lr++)
		{
			tracer_module[z].stat_stack[lr]->init();
		}
	}

	/* calculating leakage between 1st and 2nd layer, 2nd and 3rd layer ...*/
	for (int lr=1;lr<layers;lr++) /* for all layers (except the lowermost) */
	{
#pragma omp parallel for schedule(runtime) shared(lr)
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				double MaxHeadDifference, HeadDifference;

				double icevalue=ICE?ICE[i][j]:0; // valid range of icevalue now between 0...2
				if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
				else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
				if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank)&&(icevalue<1.0))
				{
					MaxHeadDifference=1.0;
					// check if there is a saturated hydraulic connection between the layers at all, 
					// otherwise the pressure head difference must be calculated differently
					if ((HEAD[lr+1][i][j]<=HEAD[lr][i][j])&&(HEAD[lr+1][i][j]<UPLIM[lr+1][i][j]))
					{
						HeadDifference=UPLIM[lr+1][i][j]-HEAD[lr][i][j];
						if (HeadDifference>0) HeadDifference=0; 
					}
					else
					{
						HeadDifference=HEAD[lr+1][i][j]-HEAD[lr][i][j];
					}
					if (HeadDifference>MaxHeadDifference) HeadDifference=MaxHeadDifference;
					else
					if (HeadDifference<-MaxHeadDifference) HeadDifference=-MaxHeadDifference;
					LEAKAGE[lr][i][j]=KOLM[lr][i][j]*HeadDifference;
				}
				else
				{
					if (MPIREG[i][j]==mpi_rank)	LEAKAGE[lr][i][j]=0.0;
				}
			}
		}
		//MPISynchronizeGridHalo(leakageindex[lr]);
		//MPISynchronizeGridHaloJustRecv(leakageindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustSend(leakageindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustWait(leakageindex[lr]); //LAURA
	}
// removed again 20160920 (was only in here to debug stripes)	if (thawdepthindex>0) MPISynchronizeGridHalo(thawdepthindex);
// removed again 20160920 (was only in here to debug stripes)	if (frozendepthindex>0) MPISynchronizeGridHalo(frozendepthindex);
// removed again 20160920 (was only in here to debug stripes)	if (iPondAvailable)	MPISynchronizeGridHalo(pondindex);

	for (int lr=1;lr<=layers;lr++) /* for all layers (except the lowermost) */
	{

#pragma omp parallel for schedule(runtime) shared(lr)
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				double icevalue = ICE ? ICE[i][j] : 0; // valid range of icevalue now between 0...2
				if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
				else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
				if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank)&&(icevalue<1.0))
				{
					HEADoldold[lr][i][j]=HEADold[lr][i][j];   /* not only the result from the last but also the forelast timestep is needed for recognizing oscillations */
					if (lr==1)
					{
						if ((iPondAvailable)&&(fabs(GWSTAND[i][j])<0.001))
						{
							HEADold[lr][i][j]=HEADolditer[lr][i][j]=HEAD[lr][i][j]=DHM[i][j]+GWSTAND[i][j]+POND[i][j]/1000.0; /* setting the new head to the interpolated groundwater level from the unsaturated zone model */
						}
						else
						{
							if ((QDSU)&&(QDSU[i][j]>0)&&(fabs(GWSTAND[i][j])<0.001))
							{
								HEADold[lr][i][j]=HEADolditer[lr][i][j]=HEAD[lr][i][j]=DHM[i][j]+GWSTAND[i][j]+QDSU[i][j]/1000.0; /* setting the new head to the interpolated groundwater level from the unsaturated zone model */
							}
							else
							{
								HEADold[lr][i][j]=HEADolditer[lr][i][j]=HEAD[lr][i][j]=DHM[i][j]+GWSTAND[i][j]; /* setting the new head to the interpolated groundwater level from the unsaturated zone model */
							}
						}
					}
					else
					HEADold[lr][i][j]=HEADolditer[lr][i][j]=HEAD[lr][i][j]; /* setting the comparison-values head to the calculated value from the last interval */
					for (int z=1;z<=tracer_count;z++) /* radioactive decay for radioactive tracers in the storage */
					if (tracer_info[z].radioactive)
						tracer_module[z].conc_stack[lr][i][j]*=tracer_info[z].decfact;
				}
				//else
				//{
				//	if (zn[i][j]>nodata)
				//	{
				//		HEADoldold[lr][i][j]=HEADold[lr][i][j]=DHM[i][j];
				//	}
				//}
			}
		}
		// synchronize all grids that will be accessed in transmissivities() for not only the actual cell (i.e. some index like [i+1][j] or [i-1][j] or [i][j-1] or i[j+1]
		MPISynchronizeGridHalo(headindex[lr]);
		//MPISynchronizeGridHaloJustRecv(headindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustSend(headindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustWait(headindex[lr]); //LAURA

		MPISynchronizeGridHalo(headoldindex[lr]);
		//MPISynchronizeGridHaloJustRecv(headoldindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustSend(headoldindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustWait(headoldindex[lr]); //LAURA

		MPISynchronizeGridHalo(headoldoldindex[lr]);
		//MPISynchronizeGridHaloJustRecv(headoldoldindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustSend(headoldoldindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustWait(headoldoldindex[lr]); //LAURA
		// removed again 20160920 (was only in here to debug stripes) MPISynchronizeGridHalo(uplimindex[lr]);
	}
	
	// synchronize all grids that will be accessed in transmissivities() for not only the actual cell (i.e. some index like [i+1][j] or [i-1][j] or [i][j-1] or i[j+1]
	if (thawdepthindex > 0)
	{
		MPISynchronizeGridHalo(thawdepthindex);
		//MPISynchronizeGridHaloJustRecv(thawdepthindex); //LAURA
		//MPISynchronizeGridHaloJustSend(thawdepthindex); //LAURA
		//MPISynchronizeGridHaloJustWait(thawdepthindex); //LAURA
	}
	if (frozendepthindex > 0)
	{
		MPISynchronizeGridHalo(frozendepthindex);
		//MPISynchronizeGridHaloJustRecv(frozendepthindex); //LAURA
		//MPISynchronizeGridHaloJustSend(frozendepthindex); //LAURA
		//MPISynchronizeGridHaloJustWait(frozendepthindex); //LAURA
	}
	if (iPondAvailable)
	{
		MPISynchronizeGridHalo(pondindex);
		//MPISynchronizeGridHaloJustRecv(pondindex); //LAURA
		//MPISynchronizeGridHaloJustSend(pondindex); //LAURA
		//MPISynchronizeGridHaloJustWait(pondindex); //LAURA
	}
	MPISynchronizeGridHalo(gwstandindex);
	//MPISynchronizeGridHaloJustRecv(gwstandindex); //LAURA
	//MPISynchronizeGridHaloJustSend(gwstandindex); //LAURA
	//MPISynchronizeGridHaloJustWait(gwstandindex); //LAURA

	/* calling the function to calculate transmissivities */
	transmissivities();

	/* 2. do while loop -> iteration of crank nicholson approximation */
	for (int lr = 1; lr <= layers; lr++)
	{
		MPISynchronizeGridHalo(flowxindex[lr]); // other row or column used for mass flux calculation later
		//MPISynchronizeGridHaloJustRecv(flowxindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustSend(flowxindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustWait(flowxindex[lr]); //LAURA

		MPISynchronizeGridHalo(flowyindex[lr]); // other row or column used for mass flux calculation later
		//MPISynchronizeGridHaloJustRecv(flowyindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustSend(flowyindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustWait(flowyindex[lr]); //LAURA

		//removed again 20160920 (was only in here to debug stripes)		MPISynchronizeGridHalo(uplimindex[lr]);
		//removed again 20160920 (was only in here to debug stripes)		if (iPondAvailable)	MPISynchronizeGridHalo(pondindex);
		MPISynchronizeGridHalo(transxindex[lr]);
		//MPISynchronizeGridHaloJustRecv(transxindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustSend(transxindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustWait(transxindex[lr]); //LAURA

		MPISynchronizeGridHalo(transyindex[lr]);
		//MPISynchronizeGridHaloJustRecv(transyindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustSend(transyindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustWait(transyindex[lr]); //LAURA

		iter = 0;
		sor_sum = 0;
		sor_anz = 0;
#ifdef _OPENMP
		if (!omp_in_parallel())
		{
			if ((((60 * stunde + (int)((minute + 60.0*hour_offset))) % (60 * nachricht) == 0) && (stunde != 0)) || (nachricht == 1))
			{
				(*WSstream) << "GWlevel " << lr << ", MPI-rank " << mpi_rank << ", Steps ";
			}
		}
#else
		if ((((60 * stunde + (int)((minute + 60.0*hour_offset))) % (60 * nachricht) == 0) && (stunde != 0)) || (nachricht == 1))
		{
			(*WSstream) << "GWlevel " << lr << ", MPI-rank " << mpi_rank << ", Steps ";
		}
#endif
		double global_maxerror = 0;
		long global_num_oscill = 0;
		double global_sum_diff = 0;
		double global_maxsum_diff = 0;

		do
		{
			iter++;
			// now gather/scatter the Halo cells for the following entities in order to update all Halo cells for all MPI regions
			// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(flowxindex[lr]);
			// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(flowyindex[lr]);
			// if (methode == 1)
			// {
			// 	MPISynchronizeGridHalo(transxindex[lr]);
			// 	MPISynchronizeGridHalo(transyindex[lr]);
			// }
			
			MPISynchronizeGridHalo(headindex[lr]);
			//MPISynchronizeGridHaloJustRecv(headindex[lr]); //LAURA
			//MPISynchronizeGridHaloJustSend(headindex[lr]); //LAURA
			//MPISynchronizeGridHaloJustWait(headindex[lr]); //LAURA

			MPISynchronizeGridHalo(headoldindex[lr]);
			//MPISynchronizeGridHaloJustRecv(headoldindex[lr]); //LAURA
			//MPISynchronizeGridHaloJustSend(headoldindex[lr]); //LAURA
			//MPISynchronizeGridHaloJustWait(headoldindex[lr]); //LAURA

			MPISynchronizeGridHalo(headoldoldindex[lr]);
			//MPISynchronizeGridHaloJustRecv(headoldoldindex[lr]); //LAURA
			//MPISynchronizeGridHaloJustSend(headoldoldindex[lr]); //LAURA
			//MPISynchronizeGridHaloJustWait(headoldoldindex[lr]); //LAURA

			// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(uplimindex[lr]);
			// other row or column used for mass flux calculation later			if (iPondAvailable)	MPISynchronizeGridHalo(pondindex);

			if (methode == 1)
			{
				maxerror = gauss_seidel_iteration(lr);
			}
			else if (methode == 2)
			{
				maxerror = HorizontalEquationSystem(lr);
			}

#ifdef _MPI
			// and also mpi-reduce maximum of maxerror
			MPI_Reduce(&maxerror, &global_maxerror, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Bcast(&global_maxerror, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
			global_maxerror=maxerror;
#endif
		} while ((iter<maxiter) && (global_maxerror>maxerrorf) && (methode == 1));
		num_oscill = get_oscillations(HEAD[lr], HEADold[lr], HEADoldold[lr]);
		// mpi-reduce num_oscill here
#ifdef _MPI
		MPI_Reduce(&num_oscill, &global_num_oscill, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&sum_diff, &global_sum_diff, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&maxsum_diff, &global_maxsum_diff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Bcast(&global_num_oscill, 1, MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&global_sum_diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&global_maxsum_diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
		global_num_oscill=num_oscill;
		global_maxsum_diff = maxsum_diff;
		global_sum_diff = sum_diff;
#endif

		if (methode == 1)
		{ // print statistics on oscillating cells for method 1 only
#ifdef _OPENMP
			if (!omp_in_parallel())
			{
				if ((((60 * stunde + (int)((minute + 60.0*hour_offset))) % (60 * nachricht) == 0) && (stunde != 0) && (minute != 0)) || (nachricht == 1))
				//if ((stunde + (int)hour_offset) % nachricht == 0)
				{
					(*WSstream) << iter << " err " << global_maxerror << ", SOR " << (auto_sor ? sor_sum / sor_anz : sor)
						<< ", osc. at " << global_num_oscill << " points (av: " <<
						(global_num_oscill > 0 ? global_sum_diff / global_num_oscill : 0.0) << ", max " << global_maxsum_diff << ")" << "[end-entry]";
				}
			}
#else
			if ((((60 * stunde + (int)((minute + 60.0*hour_offset))) % (60 * nachricht) == 0) && (stunde != 0) && (minute != 0))||(nachricht==1))
			//if ((stunde+(int)hour_offset)%nachricht==0)
			{
				(*WSstream) <<iter<<" err "<<global_maxerror<<", SOR "<<(auto_sor?sor_sum/sor_anz:sor)
					<<", osc. at "<<global_num_oscill<<" points (av: "<<
					(global_num_oscill > 0 ? global_sum_diff / global_num_oscill : 0.0) << ", max " << global_maxsum_diff << ")" << "[end-entry]";
			}
#endif
		}

	}
	for (int lr = 1; lr <= layers; lr++)
	{
		// other row or column used for mass flux calculation later					MPISynchronizeGridHalo(flowxindex[lr]);
		// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(flowyindex[lr]);
		// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(transxindex[lr]);
		// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(transyindex[lr]);
		MPISynchronizeGridHalo(headindex[lr]);
		//MPISynchronizeGridHaloJustRecv(headindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustSend(headindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustWait(headindex[lr]); //LAURA
		// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(headoldindex[lr]);
		// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(headoldoldindex[lr]);
		// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(uplimindex[lr]);
		// other row or column used for mass flux calculation later			if (iPondAvailable)	MPISynchronizeGridHalo(pondindex);
	}

	/* calculation of flux fields */
	fluxfield();

	for (int lr=1;lr<=layers;lr++)
	{
		// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(flowxindex[lr]);
		// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(flowyindex[lr]);
		// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(transxindex[lr]);
		// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(transyindex[lr]);
		MPISynchronizeGridHalo(headindex[lr]);
		//MPISynchronizeGridHaloJustRecv(headindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustSend(headindex[lr]); //LAURA
		//MPISynchronizeGridHaloJustWait(headindex[lr]); //LAURA

		// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(headoldindex[lr]);
		// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(headoldoldindex[lr]);
		// other row or column used for mass flux calculation later			MPISynchronizeGridHalo(uplimindex[lr]);
		// other row or column used for mass flux calculation later			if (iPondAvailable)	MPISynchronizeGridHalo(pondindex);
	}

	/* resulting change in groundwater table perpendicular to the gw-surface in m/s as described by Darcy's law */
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
		{
			double icevalue = ICE ? ICE[i][j] : 0; // valid range of icevalue now between 0...2
			if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
			else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank)&&(icevalue<1.0))
			{
//				TG=zarray[(int)zn[i][j]];
				/* DELTAGW is the balance of inflows and outflows and is given as a velocity with respect to vertical change of water content, not lateral! */
				/* balance for lateral inflows and outflows as needed in the unsaturated zone model */
				DELTAGW[i][j]=(HEADold[1][i][j]-HEAD[1][i][j])*S0[1][i][j]/(60.0*InV/substepcnt);
			}
			else
			{
				if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
				{
					DELTAGW[i][j]=0; // no Groundwater/soil exchange in fully glacierized cells
				}
			}
		}
	}


	// balance checks and calculations of concentrations for lateral flows in the groundwater layers and leakage between
	// gw-layers (exchange with unsaturated zone was modelled in the unsaturated zone module already) 
	
	double balance=0; // internal use for debugging: total average error of storage change vs. flux balance (should be 0 except when there are boundary conditions)
	long count=0;
	
	for (int lr=1;lr<=layers;lr++)
	{
#pragma omp parallel for schedule(runtime) \
			shared(count,lr,DT)
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			double dx, dy, massflux1, massflux2, massflux3, massflux4, massflux5, massflux6, thick0, weight,
				n_eff,thick_old,water_old,water_new,flux1,flux2,flux3,flux4,flux5,flux6;

			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				double icevalue=ICE?ICE[i][j]:0; // valid range of icevalue now between 0...2
				if (icevalue>1.0) icevalue=1.0; // this is for compatibility with the old ice grid containing 1 for ice and 2 for firn. This sequence should still work for both versions. valid range now 0...1
				else if (icevalue<0.0) icevalue=0.0; // also for compatibility: the old ice grid contianed -9999 for non-iced cells
				if ((zn[i][j]>nodata)&&(MPIREG[i][j]==mpi_rank)&&(icevalue<1.0))
				{
					if (RelCellArea!=NULL)
					{
						dx=CellSizeX[i][j];
						dy=CellSizeY[i][j];
					}
					else
					{
						dx=dy=cellsize;
					}
					n_eff = S0[1][i][j];
					/* effective thicknesses in m if aquifer is not saturated (e.g. unconfined in lower layers) */

					/* for getting the old water content, the old head is needed  */
					if (HEADold[lr][i][j]>=UPLIM[lr][i][j])
					{
						thick_old = (double)THICK[lr][i][j];  /* confined layer */
					}
					else /* unconfined layer */
					{
						thick_old = (double)HEADold[lr][i][j]-(double)UPLIM[lr][i][j]+(double)THICK[lr][i][j];
					}
					if (thick_old<0) thick_old=0;

					/* actual cell: now actual water body thicknesses are calculated */
					if (HEAD[lr][i][j]>=UPLIM[lr][i][j])
					{
						thick0    = (double)THICK[lr][i][j];  /* confined layer */
					}
					else
					{
						thick0    = (double)HEAD[lr][i][j]-(double)UPLIM[lr][i][j]+(double)THICK[lr][i][j];
					}
					if (thick0<0) thick0=0;

					water_old=dx*dy*n_eff*thick_old;

					/* massflux between local index 1 and actual cell (from [i-1,j] to [i,j]) */
					if ((i>1)&&(zn[i-1][j]!=nodata)) 
					{
						if (FLOWY[lr][i-1][j]!=nodata) 
						{
							flux1=(HEAD[lr][i-1][j]-HEAD[lr][i][j])/dy*TRANSY[lr][i-1][j]*dx*DT;
						}
						else
						{
							flux1=0.0;
						}
					}
					else flux1=0;

					/* massflux between local index 2 and actual cell (from [i,j+1] to [i,j]) */
					if ((j<ncols)&&(zn[i][j+1]!=nodata)) 
					{
						if (FLOWX[lr][i][j]!=nodata) 
						{
							flux2=-(HEAD[lr][i][j]-HEAD[lr][i][j+1])/dx*TRANSX[lr][i][j]*dy*DT;
						}
						else
						{
							flux2=0.0;
						}
					}
					else flux2=0.0;

					/* massflux between local index 3 and actual cell (between [i+1,j] and [i,j]) */
					if ((i<nrows)&&(zn[i+1][j]!=nodata)) 
					{
						if (FLOWY[lr][i][j]!=nodata) 
						{
							flux3=-(HEAD[lr][i][j]-HEAD[lr][i+1][j])/dy*TRANSY[lr][i][j]*dx*DT;
						}
						else
						{
							flux3=0.0;
						}
					}
					else flux3=0.0;

					/* massflux between local index 4 and actual cell (between [i,j-1] and [i,j] */
					if ((j>1)&&(zn[i][j-1]!=nodata)) 
					{
						if (FLOWX[lr][i][j-1]!=nodata) 
						{
							flux4=(HEAD[lr][i][j-1]-HEAD[lr][i][j])/dx*TRANSX[lr][i][j-1]*dy*DT;
						}
						else
						{
							flux4=0.0;
						}
					}
					else flux4=0;

					/* massflux from leakage between actual layer and upper layer */
					if (lr>1) 
					{
						if (LEAKAGE[lr-1][i][j]<0) 
						{ 
							flux5= -LEAKAGE[lr-1][i][j]*dx*dy*DT; 
						}
						else 
						{ 
							flux5= -LEAKAGE[lr-1][i][j]*dx*dy*DT; 
						}   
					}
					else flux5=0;

					/* massflux between actual and lower layer */
					if (lr<layers) 
					{
						if (LEAKAGE[lr][i][j]<0) 
						{ 
							flux6=LEAKAGE[lr][i][j]*dx*dy*DT;  
						}
						else 
						{ 
							flux6=LEAKAGE[lr][i][j]*dx*dy*DT; 
						}
					}
					else flux6=0;

					water_new=water_old+flux1+flux2+flux3+flux4+flux5+flux6;
					// this may have been an error during the last 9 years. I do not know why a weight was introduced at all.
					// here. It lead to rising concentrations in case the thick0-value was capped when water rises above the ground.
					// so we will try now with a weight of 1.0.
					//			weight= water_new>0?dx*dy*n_eff*thick0/water_new:1.0;
					weight=1.0;
					water_new*=weight;

					double dCheckBalance = (water_new-water_old)/(dx*dy); // balance from fluxes

					dCheckBalance = (HEAD[lr][i][j]-HEADold[lr][i][j])*S0[lr][i][j] - dCheckBalance; // subtracting the massbalance from storage change as given by heads
					dCheckBalance*=1000; // conversion into mm

					BALANCE[lr][i][j]=(double)dCheckBalance;

					for (int z=1;z<=tracer_count;z++)
					{
						mass0=tracer_module[z].conc_stack[lr][i][j]*water_old;
						/* massflux between local index 1 and actual cell (between [i-1,j] and [i,j] */
						if ((i>1)&&(zn[i-1][j]!=nodata))
						{
							if (FLOWY[lr][i-1][j]>0)
							{
								massflux1=flux1*tracer_module[z].conc_stack[lr][i-1][j];   /* flux into the actual cell, so the concentration is the one from the northern cell */
							}
							else 
							{
								if (FLOWY[lr][i-1][j]!=nodata)
								{
									massflux1=flux1*tracer_module[z].conc_stack[lr][i][j];  
								}
								else
								{
									massflux1=0.0; 
								}
							}
						} 
						else
						{
								massflux1=0.0;
						}

						/* massflux between local index 2 and actual cell (from [i,j+1] to [i,j]) */
						if ((j<ncols)&&(zn[i][j+1]!=nodata))
						{
							if (FLOWX[lr][i][j]>0)
							{
								massflux2=flux2*tracer_module[z].conc_stack[lr][i][j]; 
							}
							else
							{
								if (FLOWX[lr][i][j]>nodata)
								{
									massflux2=flux2*tracer_module[z].conc_stack[lr][i][j+1]; 
								}
								else
								{
									massflux2=0.0; 
								}
							}
						} 
						else
						{
							massflux2=0.0;
						}
						/* massflux between local index 3 and actual cell (between [i+1,j] and [i,j]) */
						if ((i<nrows)&&(zn[i+1][j]!=nodata))
						{
							if (FLOWY[lr][i][j]>0)
							{
								massflux3=flux3*tracer_module[z].conc_stack[lr][i][j]; 
							}
							else
							{
								if (FLOWY[lr][i][j]>nodata)
								{
									massflux3=flux3*tracer_module[z].conc_stack[lr][i+1][j]; 
								}
								else
								{
									massflux3=0.0; 
								}
							}
						} 
						else 
						{
							massflux3=0.0;
						}

						/* massflux between local index 4 and actual cell (between [i,j-1] and [i,j] */
						if ((j>1)&&(zn[i][j-1]!=nodata))
						{
							if (FLOWX[lr][i][j-1]>0) /* if the flow from the neigbouring cell is positive, it flows from there to the actual cell */
							{
								massflux4=flux4*tracer_module[z].conc_stack[lr][i][j-1];
							}
							else
							{
								if (FLOWX[lr][i][j-1]>nodata) /* it flows to the left */
								{
									massflux4=flux4*tracer_module[z].conc_stack[lr][i][j]; /* but the concentration is taken from the actual cell */
								}
								else
								{
									massflux4=0.0;
								}
							}
						}
						else 
						{
							massflux4=0;
						}

						/* massflux from leakage between actual layer and upper layer */
						if (lr>1)
						{
							if (LEAKAGE[lr-1][i][j]<0)
							{
								massflux5=flux5*tracer_module[z].conc_stack[lr-1][i][j];  /* inflow counts positive! */
							}
							else
							{
								massflux5=flux5*tracer_module[z].conc_stack[lr][i][j];    /* if leakage >0, this flux is a loss from the layer above, using also the concentrations of this layer */
							}
						}
						else massflux5=0;

						/* massflux between actual and lower layer */
						if (lr<layers)
						{
							if (LEAKAGE[lr][i][j]<0)
							{
								massflux6=flux6*tracer_module[z].conc_stack[lr][i][j]; 
							}
							else
							{
								massflux6=flux6*tracer_module[z].conc_stack[lr+1][i][j]; 
							}
						}
						else massflux6=0;

						/* tracer-mass-balance */
						mass0=mass0+(massflux1+massflux2+massflux3+massflux4+massflux5+massflux6)*weight;
						/* konverting into concentrations */
						tracer_module[z].concflux[i][j]= water_new>0?mass0/water_new:0.0;

						if (RelCellArea!=NULL)
						{
							tracer_module[z].stat_stack[lr]->zone[zarray[(int)zn[i][j]]].sx +=
							tracer_module[z].concflux[i][j]* RelCellArea[i][j];
							//				tracer_module[z].conc_stack[lr][i][j] * water_new * RelCellArea[i][j];
							tracer_module[z].stat_stack[lr]->zone[zarray[(int)zn[i][j]]].n += water_new*RelCellArea[i][j];
							tracer_module[z].stat_stack[lr]->zone[zarray[(int)zn[i][j]]].RelCellSum += RelCellArea[i][j];
						}
						else
						{
							tracer_module[z].stat_stack[lr]->zone[zarray[(int)zn[i][j]]].sx +=
							tracer_module[z].concflux[i][j];
							//                tracer_module[z].conc_stack[lr][i][j] * water_new ;
							tracer_module[z].stat_stack[lr]->zone[zarray[(int)zn[i][j]]].n += water_new;
							tracer_module[z].stat_stack[lr]->zone[zarray[(int)zn[i][j]]].RelCellSum += water_new;
						}
						//						fpe_check(i,j,z,lr,massflux1,massflux2,massflux3,massflux4,massflux5);
					}
				}
				else
				{
				if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
					{
						BALANCE[lr][i][j]=0.0;
					}
				}
			}
		}
		
		// for internal debugging use
		//		balance=balance/count; // total average balance for the entire basin in mm

		/* and now replacing the old concentrations with the temporal values */
		for (int i = mpistartrow; i <= mpiendrow; i++)
		{
			for (int j = FirstLastCol[i].First; j <= FirstLastCol[i].Last; j++)
			{
				if (MPIREG[i][j]==mpi_rank)
				{
					if (zn[i][j]>nodata)
					{
						HEAD[lr][i][j]=double(HEAD[lr][i][j]);
						for (int z=1;z<=tracer_count;z++)
						{
							tracer_module[z].conc_stack[lr][i][j]=tracer_module[z].concflux[i][j];
						}
					}
					else
					{
						HEAD[lr][i][j]=nodata;
						for (int z=1;z<=tracer_count;z++)
						{
							tracer_module[z].conc_stack[lr][i][j]=nodata;
						}
					}
				}
			}
		}
	}

	if ((control_row>0) && (MPIREG[control_row] != NULL) && (MPIREG[control_row][control_col] == mpi_rank))
	{
		glogbufcnt++;
		glogbuffer[glogbufcnt].j=jahr;
		//    glogbuffer[glogbufcnt].j=jahr%100;
		glogbuffer[glogbufcnt].m=monat;
		glogbuffer[glogbufcnt].t=tag;
		glogbuffer[glogbufcnt].s = stunde;
		glogbuffer[glogbufcnt].min = minute;


		for (int i=1;i<=layers;i++)
		{
			glogbuffer[glogbufcnt].wert[6*(i-1)+1] = double(HEAD[i][control_row][control_col]);
			glogbuffer[glogbufcnt].wert[6*(i-1)+2] = control_row>1?HEAD[i][control_row-1][control_col]:-8888;
			glogbuffer[glogbufcnt].wert[6*(i-1)+3] = control_col<ncols?HEAD[i][control_row][control_col+1]:-8888;
			glogbuffer[glogbufcnt].wert[6*(i-1)+4] = control_row<nrows?HEAD[i][control_row+1][control_col]:-8888;
			glogbuffer[glogbufcnt].wert[6*(i-1)+5] = control_col>1?HEAD[i][control_row][control_col-1]:-8888;
			glogbuffer[glogbufcnt].wert[6*(i-1)+6] = DELTAGW[control_row][control_col-1];
		}
		
		if ( glogbufcnt>=   glogbufmax) 
		flushbuffer(glogfile,glogfilename,&glogbufcnt,glogbuffer);

	}
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void gwflowclass::post_run(long me,long se,long tg,long mt,long jr, long actsubstep)
{
	if (do_it == 0) return;
	// deltagw was already saved in unsatzonmodel, but for coupling two model runs, the value from the groundwater model is required when starting the model again. So it has to be overwritten here
	writegridtofile(deltagwindex);
	for (int lr=1;lr<=layers;lr++)
	{
		// do not gather grids from other instances since there will be zeros in it only
		writegridtofile(headindex[lr]);
		writegridtofile(headoldindex[lr]);
//		writegridtofile(headoldoldindex[lr]);
//		writegridtofile(headolditerindex[lr]);
		writegridtofile(flowxindex[lr]);
		writegridtofile(flowyindex[lr]);
		writegridtofile(balanceindex[lr]);
		for (int z=1;z<=tracer_count;z++)
		{
			writegridtofile(tracer_module[z].storindex_stack[lr]);
			tracer_module[z].stat_stack[lr]->calc_and_done(jahr,monat,tag,stunde, minute,0);
		}
	}
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void gwflowclass::flushbuffers()
{
	int z,lr;
	if (do_it==0) return;
	else
	{
		for (int z=1;z<=tracer_count;z++)
		for (int lr=1;lr<=layers;lr++)
		tracer_module[z].stat_stack[lr]->flushbuffer();

		if (control_row > 0)
		flushbuffer(   glogfile,   glogfilename,   &glogbufcnt,   glogbuffer);
	}
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void gwflowclass::ClearMemory()
{
	int i,lr;
	if (do_it>0)
	{

		if (control_row>0)
		{
			if (initialized) flushbuffer(   glogfile,   glogfilename,   &glogbufcnt,   glogbuffer);
			for (i=1;i<=glogbufmax;i++) 
			{
				if (glogbuffer[i].wert) 
				free(glogbuffer[i].wert);
				glogbuffer[i].wert=NULL;
			}
			if (glogbuffer) free(glogbuffer);
			glogbuffer=NULL;
		}
		if (glogfilename) 
		free(glogfilename);
		glogfilename=NULL;

		int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow; // sr = start row for mpi with halo
		int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow; // er = end row for mpi with halo

		for (int lr=1;lr<=layers;lr++)
		{
			free_matrix(HEAD[lr],        gridinfo[headindex[lr]].nrl,       gridinfo[headindex[lr]].nrh,        1, ncols, nrows);  HEAD[lr]        = NULL;
			free_matrix(HEADold[lr],     gridinfo[headoldindex[lr]].nrl,    gridinfo[headoldindex[lr]].nrh,     1, ncols, nrows);  HEADold[lr]     = NULL;
			free_matrix(HEADoldold[lr],  gridinfo[headoldoldindex[lr]].nrl, gridinfo[headoldoldindex[lr]].nrh,  1, ncols, nrows);  HEADoldold[lr]  = NULL;
			free_matrix(HEADolditer[lr], gridinfo[headolditerindex[lr]].nrl,gridinfo[headolditerindex[lr]].nrh, 1, ncols, nrows);  HEADolditer[lr] = NULL;
			free_matrix(BALANCE[lr],     gridinfo[balanceindex[lr]].nrl,    gridinfo[balanceindex[lr]].nrh,     1, ncols, nrows);  BALANCE[lr]     = NULL;
			free_matrix(TRANSX[lr],      gridinfo[transxindex[lr]].nrl,     gridinfo[transxindex[lr]].nrh,      1, ncols, nrows);  TRANSX[lr]      = NULL;
			free_matrix(TRANSY[lr],      gridinfo[transyindex[lr]].nrl,     gridinfo[transyindex[lr]].nrh,      1, ncols, nrows);  TRANSY[lr]      = NULL;
			free_matrix(FLOWX[lr],       gridinfo[flowxindex[lr]].nrl,      gridinfo[flowxindex[lr]].nrh,       1, ncols, nrows);  FLOWX[lr]       = NULL;
			free_matrix(FLOWY[lr],       gridinfo[flowyindex[lr]].nrl,      gridinfo[flowyindex[lr]].nrh,       1, ncols, nrows);  FLOWY[lr]       = NULL;
			free_matrix(LEAKAGE[lr],     gridinfo[leakageindex[lr]].nrl,    gridinfo[leakageindex[lr]].nrh,     1, ncols, nrows);  LEAKAGE[lr]     = NULL;
			free_matrix(UPLIM[lr],       gridinfo[uplimindex[lr]].nrl,      gridinfo[uplimindex[lr]].nrh,       1, ncols, nrows);  UPLIM[lr]       = NULL;

			delete_grid(headindex[lr]);
			delete_grid(headoldindex[lr]);
			delete_grid(headoldoldindex[lr]);
			delete_grid(headolditerindex[lr]);
			delete_grid(balanceindex[lr]);
			delete_grid(transxindex[lr]);
			delete_grid(transyindex[lr]);
			delete_grid(flowxindex[lr]);
			delete_grid(flowyindex[lr]);
			delete_grid(leakageindex[lr]);
			delete_grid(uplimindex[lr]);

			if (HEADgridname[lr])        free(HEADgridname[lr]);         HEADgridname[lr]=NULL;
			if (HEADoldgridname[lr])     free(HEADoldgridname[lr]);      HEADoldgridname[lr]=NULL;
			if (HEADoldoldgridname[lr])  free(HEADoldoldgridname[lr]);   HEADoldoldgridname[lr]=NULL;
			if (HEADolditergridname[lr]) free(HEADolditergridname[lr]);  HEADolditergridname[lr]=NULL;
			if (FLOWXgridname[lr]) free(FLOWXgridname[lr]); FLOWXgridname[lr]=NULL;
			if (FLOWYgridname[lr]) free(FLOWYgridname[lr]); FLOWYgridname[lr]=NULL;
			if (BALANCEgridname[lr]) free(BALANCEgridname[lr]); BALANCEgridname[lr]=NULL;
		}

		if (headindex) free(headindex); headindex=NULL;
		if (headoldindex) free(headoldindex); headoldindex=NULL;
		if (headoldoldindex) free(headoldoldindex); headoldoldindex=NULL;
		if (headolditerindex) free(headolditerindex); headolditerindex=NULL;
		if (balanceindex) free(balanceindex); balanceindex=NULL;
		if (transxindex) free(transxindex); transxindex=NULL;
		if (transyindex) free(transyindex); transyindex=NULL;
		if (flowxindex) free(flowxindex); flowxindex=NULL;
		if (flowyindex) free(flowyindex); flowyindex=NULL;
		if (leakageindex) free(leakageindex); leakageindex=NULL;
		if (uplimindex) free(uplimindex); uplimindex=NULL;

		if (HEADgridname)        free(HEADgridname);        HEADgridname=NULL;
		if (HEADoldgridname)     free(HEADoldgridname);     HEADoldgridname=NULL;
		if (HEADoldoldgridname)  free(HEADoldoldgridname);  HEADoldoldgridname=NULL;
		if (HEADolditergridname) free(HEADolditergridname); HEADolditergridname=NULL;
		if (FLOWXgridname) free(FLOWXgridname); FLOWXgridname=NULL;
		if (FLOWYgridname) free(FLOWYgridname); FLOWYgridname=NULL;
		if (BALANCEgridname) free(BALANCEgridname); BALANCEgridname=NULL;

		if (FLOWXwritegrid) free(FLOWXwritegrid); FLOWXwritegrid=NULL;
		if (FLOWYwritegrid) free(FLOWYwritegrid); FLOWYwritegrid=NULL;
		if (HEADwritegrid) free(HEADwritegrid); HEADwritegrid=NULL;
		if (HEADoldwritegrid) free(HEADoldwritegrid); HEADoldwritegrid=NULL;
		if (HEADoldoldwritegrid) free(HEADoldoldwritegrid); HEADoldoldwritegrid=NULL;
		if (HEADolditerwritegrid) free(HEADolditerwritegrid); HEADolditerwritegrid=NULL;
		if (BALANCEwritegrid) free(BALANCEwritegrid); BALANCEwritegrid=NULL;

		if (TRANSX ) free(TRANSX ); TRANSX =NULL;
		if (TRANSY ) free(TRANSY ); TRANSY =NULL;
		if (BOUNDH ) free(BOUNDH ); BOUNDH =NULL;
		if (BOUNDQ ) free(BOUNDQ ); BOUNDQ =NULL;
		if (HEAD   ) free(HEAD   ); HEAD   =NULL;
		if (HEADold) free(HEADold); HEADold=NULL;
		if (HEADoldold) free(HEADoldold); HEADoldold=NULL;
		if (HEADolditer) free(HEADolditer); HEADolditer=NULL;


		if (BALANCE) free(BALANCE); BALANCE=NULL;
		if (THICK  ) free(THICK  ); THICK  =NULL;
		if (KX     ) free(KX     ); KX     =NULL;
		if (KY     ) free(KY     ); KY     =NULL;
		if (FLOWX  ) free(FLOWX  ); FLOWX  =NULL;
		if (FLOWY  ) free(FLOWY  ); FLOWY  =NULL;
		if (S0     ) free(S0     ); S0     =NULL;
		if (LEAKAGE) free(LEAKAGE); LEAKAGE=NULL;
		if (UPLIM  ) free(UPLIM  ); UPLIM  =NULL;
		if (KOLM   ) free(KOLM   ); KOLM   =NULL;

	}
	modelobject::ClearMemory();
}
//#pragma optimize ("p",off)

//#pragma optimize ("p",on)
void gwflowclass::substance_transport_ini(int &retvalue,int p_tracer_count,tracer_info_type *p_tracer_info,TLastResults *LastResults,int maxLastResults) /* initializing substance_transport algorithms (getting tracer_info table adress) */
{
//	int j,l,k,z,lr;
	int gridindex;
	char *identify;
	if (do_it==0) return;
	else
	{
		int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow; // sr = start row for mpi with halo
		int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow; // er = end row for mpi with halo

		identify = newstring(300, 1);
		tracer_count = p_tracer_count;
		if (tracer_count>0)
		{
			tracer_info = p_tracer_info;
			tracer_module = (tracer_module_type*)calloc(tracer_count+1,sizeof(tracer_module_type));
			for (int z=1;z<=tracer_count;z++)
			{
				tracer_module[z].zone_stack         = (statrecord**)calloc(layers+1,sizeof(statrecord*));
				tracer_module[z].stat_stack         = (zonalstatistic**)calloc(layers+1,sizeof(zonalstatistic*));
				tracer_module[z].statfilename_stack = (char**)calloc(layers+1,sizeof(char*));
				tracer_module[z].gridname_stack     = (char**)calloc(layers+1,sizeof(char*));
				tracer_module[z].storindex_stack    = (int*)calloc(layers+1,sizeof(int));
				tracer_module[z].fluxindex_stack    = (int*)calloc(layers+1,sizeof(int));
				tracer_module[z].conc_stack         = (double***)calloc(layers+1,sizeof(double**));
				for (int lr=1;lr<=layers;lr++)
				{
					if (tracer_info[z].radioactive==1)
					/* this is the factor, with which the actual concentration of a radioactive
			tracer has to be multiplied to get the radioactive decay for one time step */
					tracer_info[z].decfact = exp(-(InV/substepcnt/60/24)/(tracer_info[z].half_time*_1_ln2_));
					/* initializing a statistc register set for this tracer */
					tracer_module[z].zone_stack[lr] = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
					if (!tracer_module[z].zone_stack[lr]) speichererror(1);
					for (int j=0;j<=zonmax;j++)
					{ 
						tracer_module[z].zone_stack[lr][j].n         = tracer_module[z].zone_stack[lr][j].sx        = 0;
						tracer_module[z].zone_stack[lr][j].zonenwert = zone[j].zonenwert;
						tracer_module[z].zone_stack[lr][j].rel       = zone[j].rel;
						tracer_module[z].zone_stack[lr][j].RelCellSum = 0;
					}
					/* initializing the output stream for the tracer statistics */
					tracer_module[z].statfilename_stack[lr]=newstring(_stringsize_+1,1);
					sprintf(tracer_module[z].statfilename_stack[lr],"%s%.3s%.4s%.1d.%s",
					tracer_info[z].path,tracer_info[z].name,instanzname,lr,tracer_info[z].statext);
					tracer_module[z].stat_stack[lr] = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,&tracer_module[z].zone_stack[lr],(int)zonmax,
					(int)tracer_info[z].statcode,(int)99,(char*)tracer_module[z].statfilename_stack[lr],DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
					if (retvalue!=0) throw retvalue;
					sprintf(identify,"tracer concentration for tracer %d in layer %d (%s) ",z,lr,tracer_info[z].name);
					tracer_module[z].stat_stack[lr]->schreibzonenwerte(identify,instanzname,LastResults,maxLastResults);
					/* allocating memory for the concentration in a storage */
					gridindex=1;
					while (gridinfo[gridindex].identification!=NULL)
					{
						gridindex++;
						if (gridindex==gridmaxanzahl)
						{
							(*WSstream) << "FATAL: no free grid left "<<gridmaxanzahl << "[end-entry]";
							throw(12);
						}
					}
					sprintf(identify,"concstor_tracer_%d_layer_%d_%s",z,lr,instanzname);
					tracer_module[z].gridname_stack[lr]=newstring(_stringsize_+1,1);
					sprintf(tracer_module[z].gridname_stack[lr],"%s%.3s%.1d%.4s.%s",
					tracer_info[z].path,tracer_info[z].name,lr,instanzname,tracer_info[z].gridext);
					if (gridslesen!=1)
					/* Grid wird nicht von Platte eingelesen */
					{
						tracer_module[z].conc_stack[lr]=matrix(sr,er,1,ncols,nrows);
						initialize_grid(gridindex, tracer_module[z].gridname_stack[lr], tracer_module[z].conc_stack[lr], NULL, 0, identify, tracer_info[z].writegrid, sr,er);
						for (int k = mpistartrow; k <= mpiendrow; k++)
						{
							for (int l=FirstLastCol[k].First;l<=FirstLastCol[k].Last;l++) 
							{
								if ((zn[k][l] != nodata) && (MPIREG[k][l]==mpi_rank))
								{
									tracer_module[z].conc_stack[lr][k][l] = tracer_info[z].evaporating ?
									tracer_info[z].gwini :
									THICK[lr][k][l] * S0[lr][k][l] * 1000.0*tracer_info[z].gwini; /* all concentrations(masses are set to soilini */
								}
								else
								{
									tracer_module[z].conc_stack[lr][k][l] = nodata;
								}
							}
						}
						tracer_module[z].storindex_stack[lr] = gridindex;
						MPISynchronizeGridHalo(tracer_module[z].storindex_stack[lr], tracer_module[z].conc_stack[lr]);
						//MPISynchronizeGridHaloJustRecv(tracer_module[z].storindex_stack[lr], tracer_module[z].conc_stack[lr]); //LAURA
						//MPISynchronizeGridHaloJustSend(tracer_module[z].storindex_stack[lr], tracer_module[z].conc_stack[lr]); //LAURA
						//MPISynchronizeGridHaloJustWait(tracer_module[z].storindex_stack[lr], tracer_module[z].conc_stack[lr]); //LAURA

					}
					else
					{
						/* alte Dateien einlesen */
						readgridfromfile(tracer_module[z].gridname_stack[lr], &(tracer_module[z].conc_stack[lr]), 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOn_);
						initialize_grid(gridindex, tracer_module[z].gridname_stack[lr], tracer_module[z].conc_stack[lr], NULL, 0, identify, tracer_info[z].writegrid, sr,er);
						tracer_module[z].storindex_stack[lr] = gridindex;
					}
				}
				/* allocating a grid for the concentration in the outflow */
				gridindex++;
				sprintf(identify,"concflux_tracer_%d_%s",z,instanzname);
				tracer_module[z].concflux=matrix(sr,er,1,ncols,nrows);
				initialize_grid(gridindex, identify, tracer_module[z].concflux, NULL, 0, identify, 0, sr, er);
				tracer_module[z].fluxindex = gridindex;
			}
		}
		if (identify) free(identify); identify=NULL;
	}
}
//#pragma optimize ("p",off)

/* deleting all the grids for substance concentrations;
this overwrites the inherited function of teilmodellklasse */
//#pragma optimize ("p",on)
void gwflowclass::substance_transport_delete()
{
	int z,lr;
	if (do_it==0) return;
	else
	{
		int sr = (mpistartrow > 1) ? mpistartrow - 1 : mpistartrow; // sr = start row for mpi with halo
		int er = (mpiendrow < nrows) ? mpiendrow + 1 : mpiendrow; // er = end row for mpi with halo

		if (tracer_count>0)
		{
			for (int z=1;z<=tracer_count;z++)
			{
				for (int lr=1;lr<=layers;lr++)
				{
					free_matrix(tracer_module[z].conc_stack[lr], gridinfo[tracer_module[z].storindex_stack[lr]].nrl, gridinfo[tracer_module[z].storindex_stack[lr]].nrh, 1, ncols, nrows);
					delete_grid(tracer_module[z].storindex_stack[lr]);
					if (tracer_module[z].gridname_stack[lr]) 
					{
						free(tracer_module[z].gridname_stack[lr]);
						tracer_module[z].gridname_stack[lr]=NULL;
					}
					tracer_module[z].stat_stack[lr]->ClearMemory();
					if (tracer_module[z].stat_stack[lr]) delete(tracer_module[z].stat_stack[lr]);
					if (tracer_module[z].statfilename_stack[lr])
					{
						free(tracer_module[z].statfilename_stack[lr]);
						tracer_module[z].statfilename_stack[lr]=NULL;
					}
					if (tracer_module[z].zone_stack[lr]) 
					{
						free(tracer_module[z].zone_stack[lr]);
						tracer_module[z].zone_stack[lr]=NULL;
					}
				}
				free_matrix(tracer_module[z].concflux, gridinfo[tracer_module[z].fluxindex].nrl, gridinfo[tracer_module[z].fluxindex].nrh, 1, ncols, nrows);
				delete_grid(tracer_module[z].fluxindex);
				if (tracer_module[z].zone_stack        ) 
				{
					free(tracer_module[z].zone_stack        );
					tracer_module[z].zone_stack        =NULL;
				}
				if (tracer_module[z].stat_stack        ) 
				{
					free(tracer_module[z].stat_stack        );
					tracer_module[z].stat_stack        =NULL;
				}
				if (tracer_module[z].statfilename_stack) 
				{
					free(tracer_module[z].statfilename_stack);
					tracer_module[z].statfilename_stack=NULL;
				}
				if (tracer_module[z].gridname_stack    ) 
				{
					free(tracer_module[z].gridname_stack    );
					tracer_module[z].gridname_stack    =NULL;
				}
				if (tracer_module[z].storindex_stack   ) 
				{
					free(tracer_module[z].storindex_stack   );
					tracer_module[z].storindex_stack   =NULL;
				}
				if (tracer_module[z].fluxindex_stack   ) 
				{
					free(tracer_module[z].fluxindex_stack   );
					tracer_module[z].fluxindex_stack   =NULL;
				}
				if (tracer_module[z].conc_stack        ) 
				{
					free(tracer_module[z].conc_stack        );
					tracer_module[z].conc_stack        =NULL;
				}
			}
		}
	}
}
//#pragma optimize ("p",off)


