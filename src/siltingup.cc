/*   $Revision: 23 $ 
	$Date: 7.03.11 0:36 $  */

/*
(c) 2007 by Joerg Schulla

*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "stringp.h"   /* String-funktionen */
#include "modconst.h"
#include "modef.h"     /* Definition des Ausgangsobjectes der Interpolationen und Teilmodelle, Code in modef.C */
#include "siltingup.h"
#include "wasim.h"

using namespace std;


siltingupclass::siltingupclass(int &retvalue,wasim::mystream *Stream,TLastResults *LastResults, int maxLastResults,
	double breitearg,double laengearg,double bezugsmeridianarg,char *bezeichnung,
	char *steuername,int LogLevel,int progaus,
	double p_min_to_int_end,int p_numsubb,int *p_listsubb,
	gridinfotype *p_gridinfo,
	int *p_zarray,statrecord *p_zone,firstlast *p_FirstLastCol ,
	int p_zonmax,int p_zonmaxwert,
	int p_ncols, int p_nrows, int p_mpistartrow, int p_mpiendrow,
	double p_xll,double p_yll,
	double p_cellsize,double p_nodata,
	long mine,long ste , long tge , long mte , long jre, double p_timekorr,
	LUSEtabletype *p_LUtable,
	SOILtabletype *p_SOILtable,
	MultiLayerLUSEType *p_MultiLayerLUTable,
	bool WriteASCIIGrids, bool WriteDoubleGrids, bool WriteNetcdfGrids, char *DefaultInputDirectory, char *DefaultOutputDirectory,
	long *v_mpiregioncnt,std::vector<mpihalo> *p_sendHalo,std::vector<mpihalo> *p_recvHalo,
	int p_doSpinUp, bool p_routingmode, std::vector<HRUMapElement*> *p_vHRU):
teilmodellklasse(Stream,breitearg,laengearg,bezugsmeridianarg,bezeichnung,steuername,LogLevel,progaus,p_min_to_int_end,p_numsubb,p_listsubb,
	p_gridinfo,p_zarray,p_zone,p_FirstLastCol,p_zonmax,p_zonmaxwert,
	p_ncols,p_nrows,p_mpistartrow,p_mpiendrow,p_xll,p_yll,p_cellsize,p_nodata,mine,ste,tge,mte,jre,p_timekorr,
	p_LUtable,p_SOILtable,p_MultiLayerLUTable,
	WriteASCIIGrids, WriteDoubleGrids, WriteNetcdfGrids, DefaultInputDirectory, DefaultOutputDirectory,
	v_mpiregioncnt,p_sendHalo,p_recvHalo,
	p_doSpinUp,p_routingmode,p_vHRU),
// initializing member variables here to 0/NULL or whatever...
method(0),
DHM(0),P(0),LUSE(0),SOIL(0),SDISP(0),RAINPAUSE(0),EKIN(0),TSBB(0),QDSU(0),
RainPause(0),code(0),
SDISPwritegrid(0),RAINPAUSEwritegrid(0),EKINwritegrid(0),TSBBwritegrid(0),QDSUwritegrid(0),
SDISPgridname(0),RAINPAUSEgridname(0),EKINgridname(0),TSBBgridname(0),QDSUgridname(0),
sdispindex(0),rainpauseindex(0),ekinindex(0),tsbbindex(0),qdsuindex(0),
SDISPstat(0),QDSUstat(0),zoneSDISP(0),zoneQDSU(0),SDIPSstatcode(0),QDSUstatcode(0),
SDISPstatfilename(0),QDSUstatfilename(0),iNumVegetationLayers(0),SurfaceRouting(false),
SurfaceRoutingChecked(false),
parA(65),
parB(12.2),
parC(0.52),
parD(-0.64),
parE(0.013),
parF(-1.03),
parG(0.7),
parH(-0.19),
interceptionBeforeSnow(false)
{
	try
	{
		retvalue=0;
		int i,j,k;

		/* ACHTUNG !!!!!: die Steuerdatei wird im Constructor der MutterKlasse (teilmodell)
	geoeffnet und bis zum Erreichen des

	ignoranzparameters 
	+ Intervallaenge

	gelesen.
	ist der Ignoranzparameter = 0, dann wird die Datei nach diesem Parameter (InV) nicht
	weitergelesen und gleich geschlossen.
	So kann der gesamte Ablauf der ini-Prozedur abgebrochen werden, da der
	Ignoranzparameter schon bekannt ist!
	*/

		if (routingmode)
		{
			do_it = 0; // wether the submodels code for do_it is set to 0 or 1: when routingmode is active, this submodel should never run
		}

		if (do_it == 0) return;

#ifdef __TIME__
		(*WSstream) << "file "<<__FILE__<<", compiled at "<<__DATE__<<", "<<__TIME__<<""<<"[end-entry]";
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

			if (ist_vorhanden(_MPIHelper_))
			{
				MPIHLP = gridinfo[getindex(_MPIHelper_)].pointer;
			}
			else
			{
				MPIHLP = NULL;
			}
		}

		if (ist_vorhanden(_GlacierizedCells_) || ist_vorhanden(_GlacierCodes_))
		{
			(*WSstream) << "FATAL: error: Silting-Up module cannot be used when glacier model is active. Remove grids for glacier codes and glaciated cell fractions from the grid list"<<"[end-entry]";
			throw(2);
		}


		// getting addresses of existing grids (standardgrids or interpolation results etc.)

		// elevation model
		DHM = gridinfo[2].pointer; /* Hoehenmodell von jetzt an unter DHM ansprechbar */


		// Multilayer-Land use or landuse grid
		if (MultiLayerLUTable)
		{
			iNumVegetationLayers=MultiLayerLUTable[1].getNumLayers();
		}
		else
		{
			iNumVegetationLayers=1;
		}
		LUSE   = new double**[iNumVegetationLayers+1];
		char name[500];
		int lay;
		for (lay=1;lay<=iNumVegetationLayers;lay++)
		{
			sprintf(name,"%s%d",_Landnutzung_,lay);
			LUSE[lay]   = gridinfo[getindex(name)].pointer;
		}

		// VCF
		VCF = new double**[iNumVegetationLayers+1];
		for (lay=1;lay<=iNumVegetationLayers;lay++)
		{
			sprintf(name,"%s%d",_Vegetationsbedeckungsgrad_,lay);
			VCF[lay]   = gridinfo[getindex(name)].pointer;
			//      VCF[lay]    = gridinfo[getindex(_Vegetationsbedeckungsgrad_)].pointer; // vegetation covered fraction of soil
		}

		// soil types
		SOIL  =  gridinfo[getindex(_Bodenarten_)].pointer; // soil types (texture)

		// searching surface storage of silting up module which is used in surface routing, too
		SSSU  =  gridinfo[getindex(_SURFSTORSILTINGUP_)].pointer;

		// searching precipitation grid (input into soil model from other model components
		//   -> if interception module was run: throughfall)
		//   -> if not: snowmelt
		//   -> if even no snow model was run: precipitation

// is done in method setInterceptionBeforeSnow now because sequence of snow model and interception mode may be reversed
//		sprintf(name,"%s%d",_Bestandesniederschlag_,0);
//		P     =  gridinfo[getindex_first_of_3(name,_Schneeausfluss_,_Niederschlag_)].pointer;

		// reserving memeory for the filenames
		SDISPgridname		= newstring(_stringsize_+1,1);
		RAINPAUSEgridname	= newstring(_stringsize_+1,1);
		EKINgridname		= newstring(_stringsize_+1,1);
		TSBBgridname		= newstring(_stringsize_+1,1);
		QDSUgridname		= newstring(_stringsize_+1,1);

		SDISPstatfilename	= newstring(_stringsize_+1,1);
		QDSUstatfilename	= newstring(_stringsize_+1,1);

		/* lesen der Werte aus der Steuerdatei */
		// todo: check if the next string in the control file is a number. If yes, this is the method, otherwise method is set to 1 by default
		fgets(tstr,_stringsize_,steuerdatei);
		if (sscanf(tstr,"%d",&method)!=1)
		{
			this->method = 1; // default
			if (sscanf(tstr,"%s%d",SDISPstatfilename,&SDIPSstatcode)!=2)
			{
				steuererror("name of the statistics for silting disposition SDISP ",this->instanzname);
			}
		}
		else 
		{   
			fscanfln(steuerdatei,"%s%d",SDISPstatfilename,&SDIPSstatcode ,"name of the statistics for silting disposition SDISP ",tstr);
		}

		fscanfln(steuerdatei,"%s%d",QDSUstatfilename,&QDSUstatcode ,"name of the statistics for direct runoff from silting up model (QDSU) ",tstr);

		fscanfln(steuerdatei,"%s",SDISPgridname,"name of the silting-up-disposition-grid, "    ,tstr);
		fscanfln(steuerdatei,"%d",&SDISPwritegrid,"Writegrid for silting-up-disposition-grid, ",tstr);

		fscanfln(steuerdatei,"%s",RAINPAUSEgridname,"name of the RAINPAUSE-grid, "    ,tstr);
		fscanfln(steuerdatei,"%d",&RAINPAUSEwritegrid,"Writegrid for RAINPAUSE-grid, ",tstr);

		fscanfln(steuerdatei,"%s",EKINgridname,"name of the EKIN-grid, "    ,tstr);
		fscanfln(steuerdatei,"%d",&EKINwritegrid,"Writegrid for EKIN-grid, ",tstr);

		fscanfln(steuerdatei,"%s",TSBBgridname,"name of the TSBB-grid, "    ,tstr);
		fscanfln(steuerdatei,"%d",&TSBBwritegrid,"Writegrid for TSBB-grid, ",tstr);

		fscanfln(steuerdatei,"%s",QDSUgridname,"name of the QDSU-grid, "    ,tstr);
		fscanfln(steuerdatei,"%d",&QDSUwritegrid,"Writegrid for QDSU-grid, ",tstr);

		// reading the zone-specific parameters. Step 1: reading zone codes, Step 2: reading each parameter
		code      = (int*)calloc(zonmax+1,sizeof(int)); // codes of the zones corresponding to the zone grid
		RainPause      = (double*)calloc(zonmax+1,sizeof(double)); // sample parameter
		//    ParB      = (double*)calloc(zonmax+1,sizeof(double)); // sample parameter
		k=zone[1].zonenwert; // first code is minmal for the lowest basin number
		do // reading ranges for subbasin codes and parameters 
		{
			if (fscanf(steuerdatei,"%d",&i)!=1)
			{
				(*WSstream) << "FATAL: \n-----------------error------------\n"
				"error reading code["<<(int)i<<"] for "<<instanzname<<"\n"
				"possible error by using the same parameters\n"
				"for groups of subbasins:\n"
				"at least the first AND the last subbasin code have\n"
				"to be specified with all their parameters like this:\n"
				"(if the smallest code is 1 and the largest code is 8)\n"
				"1\t8\n3\t3\n8\t8\n14\t14\n0.4\t0.4\n0.01\t0.01\n0.1\t0.1\n..\t..\n"
				"Then, to all subcatchments with zone-codes ranging from 1 to 7 the \n"
				"parameters of subbasin 1 will be assigned!"<<"[end-entry]";
				throw(10);
			}

			if ((i>zonmaxwert)||(i<1)||(zarray[i]==0))
			{
				(*WSstream) << "FATAL: \n-----------------error------------\n"
				"error! there is no code "<<i<<" (from "<<steuerdateiname
				<<") for a subbasin in "<<instanzname<<"\n"
				"possible error by using the same parameters\n"
				"for groups of subbasins:\n"
				"at least the first AND the last subbasin code have\n"
				"to be specified with all their parameters like this:\n"
				"(if the smallest code is 1 and the largest code is 8)\n"
				"1\t8\n3\t3\n8\t8\n14\t14\n0.4\t0.4\n0.01\t0.01\n0.1\t0.1\n..\t..\n"
				"Then, to all subcatchments with zone-codes ranging from 1 to 7 the \n"
				"parameters of subbasin 1 will be assigned!"<<"[end-entry]";
				throw(10);
			}
			else
			code[zarray[i]]=i;

			if (zarray[i]>zarray[k]+1) // filling codes inbetween with -1 
			if ((k>zone[1].zonenwert)||((k==zone[1].zonenwert)&&(code[1]>0)))
			for (j=zarray[k]+1;j<zarray[i];j++) code[j]=-1;
			else
			{
				code[zarray[k]]=zone[1].zonenwert;
				for (j=zarray[k]+1;j<=zarray[i];j++) code[j]=-1; // moving the first valid subbasin code to the first subbasin, if first code was higer
			}
			k=i;
		}
		while (i<zonmaxwert);
		tstr=fgets(tstr,_stringsize_,steuerdatei);

		// reading RainPause-set 
		for (i=1;i<=zonmax;i++)  
		{  
			if (code[i]>0)
			{
				if (fscanf(steuerdatei,"%lf",&(RainPause[i]))!=1)
				{ 
					(*WSstream) << "FATAL: \n-----------------error------------\n"
					"error reading RainPause["<<(int)code[i]<<"] for "<<instanzname<<" (minimum rainpause)\n"
					"possible error by using the same parameters\n"
					"for groups of subbasins:\n"
					"at least the first AND the last subbasin code have\n"
					"to be specified with all their parameters like this:\n"
					"(if the smallest code is 1 and the largest code is 8)\n"
					"1\t8\n3\t3\n8\t8\n14\t14\n0.4\t0.4\n0.01\t0.01\n0.1\t0.1\n..\t..\n"
					"Then, to all subcatchments with zone-codes ranging from 1 to 7 the \n"
					"parameters of subbasin 1 will be assigned!"<<"[end-entry]";
					throw(10);
				}
			} else RainPause[i]=RainPause[i-1];
		}  
		tstr=fgets(tstr,_stringsize_,steuerdatei);

		// reading ParB-set 
		/*    for (i=1;i<=zonmax;i++)  
	{  
	if (code[i]>0)
	{
		if (fscanf(steuerdatei,"%lf",&(ParB[i]))!=1)
		{ 
		(*WSstream) << "FATAL: \n-----------------error------------\n"
				"error reading ParB["<<(int)code[i]<<"] for "<<instanzname<<" (minimum Water depth for crop roughness)\n"
				"possible error by using the same parameters\n"
				"for groups of subbasins:\n"
				"at least the first AND the last subbasin code have\n"
				"to be specified with all their parameters like this:\n"
				"(if the smallest code is 1 and the largest code is 8)\n"
				"1\t8\n3\t3\n8\t8\n14\t14\n0.4\t0.4\n0.01\t0.01\n0.1\t0.1\n..\t..\n"
				"Then, to all subcatchments with zone-codes ranging from 1 to 7 the \n"
				"parameters of subbasin 1 will be assigned!"<<"[end-entry]";
		throw(10);
		}
	} else ParB[i]=ParB[i-1];
	}  
	tstr=fgets(tstr,_stringsize_,steuerdatei);
*/
		// reading the value of readgrids from control file
		if (fscanfln(steuerdatei,"%d",&gridslesen   ,"read control parameter for storage grids for ",tstr)!=1)
		{
			(*WSstream) << "FATAL: \n-----------------error------------\n"
			"error reading readgrid-Parameter for "<<instanzname<<"\n"
			"possible reason: parameter missing or a linefeed \nbefore the end of the parameter block"<<"[end-entry]";
			throw(10);
		}

		// reading parameters A to H for method 2
		if (method==2)
		{	
			fscanfln(steuerdatei,"%lf",&parA,"parameter A for method 2 (I_0) ",tstr);
			fscanfln(steuerdatei,"%lf",&parB,"parameter B for method 2 (I_inf) ",tstr);
			fscanfln(steuerdatei,"%lf",&parC,"parameter C for method 2 (I_inf) ",tstr);
			fscanfln(steuerdatei,"%lf",&parD,"parameter D for method 2 (I_inf) ",tstr);
			fscanfln(steuerdatei,"%lf",&parE,"parameter E for method 2 (C_v) ",tstr);
			fscanfln(steuerdatei,"%lf",&parF,"parameter F for method 2 (C_v) ",tstr);
			fscanfln(steuerdatei,"%lf",&parG,"parameter G for method 2 (C_v) ",tstr);
			fscanfln(steuerdatei,"%lf",&parH,"parameter H for method 2 (C_v) ",tstr);
		} else if (method==3)	{
			// reading the expressions
			// by first looking for a keyword "exprsssion" and then looking for matching curled brackets {}
			bool expressionfound=false;
			expressionlist.clear();
			do 
			{
				tstr=fgets(tstr,_stringsize_,steuerdatei);
				std::string record(tstr);
				if (record.find("SiltingUpExpressions")==0) // a line starting with the searched key word eas found
				{
					bool bracketfound=false;
					std::string subrecord=record.substr(20);
					if (record.find("{")!=std::string::npos)
					{
						bracketfound=true;
					}
					while ((record.find("{")==std::string::npos)&& (!feof(steuerdatei)) && (subrecord.find_first_not_of(" \n\t\r")==std::string::npos))
					{
						tstr=fgets(tstr,_stringsize_,steuerdatei);
						record.append(std::string(tstr));
						subrecord=std::string(tstr);
						if (record.find("{")!=std::string::npos)
						{
							bracketfound=true;
						}
					}
					if (feof(steuerdatei)||(!bracketfound))
					{
						(*WSstream) << "FATAL: error looking for opening bracket \"{\" after \"SiltingUpExpressions\" in control file \""<<steuerdateiname<<"\""<<"[end-entry]";
						throw(400);
					}
					// no looking for the closing bracket
					bracketfound=false;
					if (record.find("}")!=std::string::npos)
					{
						bracketfound=true;
					}
					while ((record.find("}")==std::string::npos) && (!feof(steuerdatei)) && (record.find_first_of(" [{")!=std::string::npos))
					{
						tstr=fgets(tstr,_stringsize_,steuerdatei);
						record.append(std::string(tstr));
						if (record.find("}")!=std::string::npos)
						{
							bracketfound=true;
						}
					}
					if (feof(steuerdatei))
					{
						(*WSstream) << "FATAL: error looking for closing bracket \"}\" after \"SiltingUpExpressions\" in control file \""<<steuerdateiname<<"\""<<"[end-entry]";
						throw(401);
					}
					record = record.substr(record.find("{")+1);
					record = record.substr(0,record.find("}"));
					if (record.find(";")!=std::string::npos)
						expressionfound=true;
					else
						expressionfound=false;
					if (expressionfound==true)
					{
						// now parsing the expressionlist
						record=record+std::string(" "); // make sure that at least one character follows the last semicolon
						while (record.find(";")!=std::string::npos)
						{
							std::string expression;
							int pos1,pos2;
							pos1=record.find_first_not_of(" \n\r\t");
							pos2=record.find(";");
							expression=record.substr(pos1,pos2-pos1);
							pos1 = expression.find_first_of(" \t\n\r");
							while (pos1!=std::string::npos)
							{
								expression=expression.substr(0,pos1)+expression.substr(pos1+1);
								pos1 = expression.find_first_of(" \t\n\r");
							}
							this->expressionlist.push_back(expression);
							// delete the actual expression from the original record
							record=record.substr(record.find(";")+1);
							// now, the expressions are stored in the expressionlist, 
							// where they will be taken from each time the expressions are processed				
						}
					}
				}
			}
			while ((!feof(steuerdatei))&&(index(tstr,"\n")>0)&&(!expressionfound));
			if (!expressionfound)
			{
				// clear the list if expressions or so whatever is required to ensure that the model will run properly, e.g. throwing an exception 
				// saying that at least X, Y and Z must be set by the expressions
				(*WSstream) << "FATAL: no SiltingUpExpressios found (at least X, Y and Z must be initialized) in control file \""<<steuerdateiname<<"\""<<"[end-entry]";
				throw(402);
			}
		}

		// now we can close the control file, since the last parameter was read successfully
		if ((fclose(steuerdatei))==EOF)
		{
			(*WSstream) << "FATAL: error closing control file \""<<steuerdateiname<<"\""<<"[end-entry]";
			throw(11);
		}

		// initializing the statistics registers for all statistics	tables we use in this module

		// statistics registers for siting up disposition
		zoneSDISP = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneSDISP) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{ 
			zoneSDISP[i].n           = zoneSDISP[i].sx        = 0;
			zoneSDISP[i].zonenwert   = zone[i].zonenwert;
			zoneSDISP[i].rel         = zone[i].rel;
			zoneSDISP[i].RelCellSum  = zone[i].RelCellSum;
		}

		// statistics registers for direct discharge from silting up module
		zoneQDSU = (statrecord *)calloc(zonmax+1,sizeof(statrecord));
		if (!zoneQDSU) speichererror(1);
		for (i=0;i<=zonmax;i++)
		{ 
			zoneQDSU[i].n           = zoneQDSU[i].sx        = 0;
			zoneQDSU[i].zonenwert   = zone[i].zonenwert;
			zoneQDSU[i].rel        = zone[i].rel;
			zoneQDSU[i].RelCellSum  = zone[i].RelCellSum;
		}

		// initializing the statistics class for each register set we use in this module
		// silting up disposition
		SDISPstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
		&zoneSDISP,(int)zonmax,(int)SDIPSstatcode,(int)4,(char*)SDISPstatfilename,
		DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		SDISPstat->schreibzonenwerte("Silting_up_disposition_",instanzname,LastResults,maxLastResults);

		// direct discharge from silting up module
		QDSUstat = new zonalstatistic(retvalue,WSstream,LastResults,maxLastResults,RelCellArea,nrows,ncols,
		&zoneQDSU,(int)zonmax,(int)QDSUstatcode,(int)4,(char*)QDSUstatfilename,
		DefaultOutputDirectory,(int)numsubb,(int*)listsubb, (InV < 60));
		if (retvalue!=0) throw retvalue;
		QDSUstat->schreibzonenwerte("Direct_Discharge_Silting_Up",instanzname,LastResults,maxLastResults);


		SDISP		=  NULL;
		RAINPAUSE	=  NULL;
		EKIN		=  NULL;
		TSBB		=  NULL;
		QDSU		=  NULL;


		// looking for the first free entry in the gridinfo array
		int sti=1;
		while (gridinfo[sti].identification!=NULL)
		{
			sti ++;
			if (sti==gridmaxanzahl)
			{
				(*WSstream) << "FATAL: no free grid left, allowed number of grids is (contact the author)"<<gridmaxanzahl<<"!"<<"[end-entry]";
				throw(12);
			}
		}

		if (gridslesen!=1)
		{
			// only allocate memeory and initialize gridinfo-array and the grid itself to 0
			SDISP = matrix(mpistartrow,mpiendrow,1,ncols,nrows);
			initialize_grid(sti, SDISPgridname, SDISP, NULL, 0, _SILTING_UP_DISPOSITION_, SDISPwritegrid, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					SDISP[k][j] = zn[k][j] == nodata ? nodata : 0;
				}
			}
			sdispindex = sti;
			sti ++;

			RAINPAUSE = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(sti, RAINPAUSEgridname, RAINPAUSE, NULL, 0, _RAINPAUSE_, RAINPAUSEwritegrid, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					RAINPAUSE[k][j] = zn[k][j] == nodata ? nodata : 0;
				}
			}
			rainpauseindex = sti;
			sti ++;

			EKIN = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(sti, EKINgridname, EKIN, NULL, 0, _KINETIC_ENERGY_, EKINwritegrid, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					EKIN[k][j] = zn[k][j] == nodata ? nodata : 0;
				}
			}
			ekinindex = sti;
			sti ++;

			TSBB = matrix(mpistartrow, mpiendrow, 1, ncols, nrows);
			initialize_grid(sti, TSBBgridname, TSBB, NULL, 0, _TIMESINCELASTSOILTILLAGE_, TSBBwritegrid, mpistartrow, mpiendrow);
			for (int k = mpistartrow; k <= mpiendrow; k++)
			{
				for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
				{
					TSBB[k][j] = zn[k][j] == nodata ? nodata : 0;
				}
			}
			tsbbindex = sti;
			sti ++;
		}
		else
		{
			// read existing files and initialize the gridinfo array for this grid
			readgridfromfile(SDISPgridname, &SDISP, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
			initialize_grid(sti, SDISPgridname, SDISP, NULL, 0, _SILTING_UP_DISPOSITION_, SDISPwritegrid, mpistartrow, mpiendrow);
			sdispindex = sti;
			sti ++;

			readgridfromfile(RAINPAUSEgridname, &RAINPAUSE, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
			initialize_grid(sti, RAINPAUSEgridname, RAINPAUSE, NULL, 0, _RAINPAUSE_, RAINPAUSEwritegrid, mpistartrow, mpiendrow);
			rainpauseindex = sti;
			sti ++;

			readgridfromfile(EKINgridname, &EKIN, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
			initialize_grid(sti, EKINgridname, EKIN, NULL, 0, _KINETIC_ENERGY_, EKINwritegrid, mpistartrow, mpiendrow);
			ekinindex = sti;
			sti ++;

			readgridfromfile(TSBBgridname, &TSBB, 1, _AllocateGridMemoryOn_, _ShowMessagesOn_, 0, _LimitAutoFillOn_, _MPIFullGridOff_, _MPIRank0OnlyOff_, _MPISyncHaloOff_);
			initialize_grid(sti, TSBBgridname, TSBB, NULL, 0, _TIMESINCELASTSOILTILLAGE_, TSBBwritegrid, mpistartrow, mpiendrow);
			tsbbindex = sti;
			sti ++;

		}

		// reserving memeory for direct discharge (will not read from grid since it is not a storage which has a lifetim for multiple time steps
		QDSU = matrix(mpistartrow, mpiendrow, 1, ncols,nrows);
		initialize_grid(sti, QDSUgridname, QDSU, NULL, 0, _DIRECTRUNOFFSILTINGUP_, QDSUwritegrid, mpistartrow, mpiendrow);
		for (int k = mpistartrow; k <= mpiendrow; k++)
		{
			for (int j = FirstLastCol[k].First; j <= FirstLastCol[k].Last; j++)
			{
				QDSU[k][j] = zn[k][j] == nodata ? nodata : 0;
			}
		}
		qdsuindex = sti;

		if (SDISPgridname) free(SDISPgridname); SDISPgridname=NULL;
		if (RAINPAUSEgridname) free(RAINPAUSEgridname); RAINPAUSEgridname=NULL;
		if (EKINgridname) free(EKINgridname); EKINgridname=NULL;
		if (TSBBgridname) free(TSBBgridname); TSBBgridname=NULL;
		if (QDSUgridname) free(QDSUgridname); QDSUgridname=NULL;
		
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

void siltingupclass::setInterceptionBeforeSnow(bool b_interceptionBeforeSnow)
{
	if (do_it == 0) return;
	interceptionBeforeSnow = b_interceptionBeforeSnow;
	char name[500];
	if (!P)
	{
		if (interceptionBeforeSnow) // new approach: first look for snow melt, then for interception model output and then for interpolated precipitation
		{
			sprintf(name, "%s%d", _Bestandesniederschlag_, 0);
			P = gridinfo[getindex_first_of_3(_Schneeausfluss_,name, _Niederschlag_)].pointer;
		}
		else // classic approach: first look for interception model output, then for snow melt and then for interpolated precipitation
		{
			sprintf(name, "%s%d", _Bestandesniederschlag_, 0);
			P = gridinfo[getindex_first_of_3(name, _Schneeausfluss_, _Niederschlag_)].pointer;
		}
	}
}

void siltingupclass::run(long me,long se,long tg,long mt,long jr)
{
	int z;
	// print a message to the output stream and register the actual date and time
	teilmodellklasse::run(me,se,tg,mt,jr,NULL,NULL);

	if (do_it == 0) return;

	if (!SurfaceRoutingChecked)
	{
		SurfaceRouting=false;
		SurfaceRoutingChecked=false;
		if (ist_vorhanden(_SURFACEROUTING_WEIGHTS_))
		{
			SurfaceRouting=true;
		}
		else
		{
			SurfaceRouting=false;
		}
		SurfaceRoutingChecked=true;
	}


	const double log004    = log10(0.04);
	const double log2      = log10(2.0);
	const double log6_3    = log10(6.3);
	const double log20     = log10(20.0);
	const double log36     = log10(36.0);
	const double log63     = log10(63.0);
	const double log100    = log10(100.0);
	const double log200    = log10(200.0);
	const double log630    = log10(630.0);
	const double log2000   = log10(2000.0);
	const double log6300   = log10(6300.0);
	const double log20000  = log10(20000.0);
	const double log63000  = log10(63000.0);
	const double log200000 = log10(200000.0);

	const double FClay   = log004+log2;
	const double FSilt   = 0.3326 * (log2+log6_3) +  0.3348 *(log6_3+log20)  + 0.1704 * (log20+log36)  + 0.1622 * (log36+log63) ;
	const double FSand   = 0.1336 * (log63+log100) + 0.2005 *(log100+log200) + 0.3318 *(log200+log630) + 0.3341 *(log630+log2000);
	const double FStones1 = (log2000+log6300);
	const double FStones2 = (log6300+log20000);
	const double FStones3 = (log20000+log63000);
	const double FStones4 = (log63000+log200000);


	// initializing statistics register
	SDISPstat->init(); 
	QDSUstat->init();
	for (z=1;z<=tracer_count;z++) tracer_module[z].stat->init();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//																											   
	// now processing the central loop
	//																											   
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma omp parallel for schedule(runtime) shared(tg,mt,jr)  
	// do not use any common variable in the openmp-loop for writing! multiple threads are accessing the same variable at the same time!
	// all variables with write access must be listed in the private-part of the openmp-clause. If the variables are even of interest 
	// for the results outside the loop, then the accumulation must be done using the suitable openmp-directives (or do it manually as 
	// we do it with the zonal statistics by doing such accumulation in a seperate loop after the central processing.
	for (int i=mpistartrow;i<=mpiendrow;i++) /* Reihen des Grids */
	{
		Variablelist *varlist;
		char *parserresult;
		Parser *parser;
		varlist=NULL;
		parser=NULL;
		parserresult=NULL;

		if (method==3)
		{
			// create a parser for each thread in order to separate varlist entries
			varlist = new Variablelist();
			parser = new Parser(varlist);
		}
		for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) /* Spalten des Grids */
		{
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank) && !(doSpinUp && HRUMaster[i][j]==nodata))
			//		  if (zn[i][j]!=nodata)
			{
				int TG=zarray[(int)zn[i][j]];
				if (!SurfaceRouting){
					SSSU[i][j]=0; // when no surface routing is active, all direct runoff was routed into the single linear storage cascade in unsatzonmodel, so SSSU must be set to 0 here
					// otherwise direct runoff would add up here. When the model runs with surface routing, the runoff of the water in SSSU is determined by the surface routing algorithm
				}
				// step 1: balancing tracer concentrations in the surface runoff storage with incoming fluxes
				// and apply radioactive decay, if required for radiactive tracers
				for (int z=1;z<=tracer_count;z++)
				{
					if (tracer_info[z].radioactive)
					tracer_module[z].concstor[i][j]*=tracer_info[z].decfact;
					if (P[i][j]>0) // no change in concentrations without new input
					{
						// evaporating tracers will be handled as concentrations, all other tracers will be handled as absolute masses
						tracer_module[z].concstor[i][j] = tracer_info[z].evaporating ? 
						mixing(tracer_module[z].concstor[i][j],tracer_module[z].cinput[i][j],
						SSSU[i][j],P[i][j]) :
						tracer_module[z].concstor[i][j]+tracer_module[z].cinput[i][j]*P[i][j]; 
					}
				}
				// get the array, the soil tillage days are stored in
				int* iSoilTillageDays;
				int  iSoilTillageCount;
				LUtable[(int)LUSE[1][i][j]].GetSoilTillage(iSoilTillageCount,iSoilTillageDays);
				// iSoilTillageDays[0] will contain a negative number, which is the last entry in the array - 366, so it is the day of the last year, the last soil tillage took place
				// all entries may be accessed by iSoilTillageDays[i] with i ranging from 0..iSoilTillageCount.
				// to get the last soil tillage event, the actual day has to be converted into a Julian day first.
				// do the central processing here for water and tracer fluxes

				// get the fractions of sand, silt and clay
				double dFractionSand,
				dFractionSilt,
				dFractionClay,
				dFractionStones1,
				dFractionStones2,
				dFractionStones3,
				dFractionStones4;
				// get local landuse code luse, vegetation covered fraction of soil vcf, and soil type s
				int maxlayer = iNumVegetationLayers;
				while ((int)LUSE[maxlayer][i][j] == nodata && maxlayer>1) { maxlayer--; };
				int luse = (int)LUSE[maxlayer][i][j];
				double vcf = (double)VCF[maxlayer][i][j];
				double ResidCovrg = LUtable[(int)LUSE[1][i][j]].GetSR_ResidCovrg(WSstream,tg,mt,jr,DHM[i][j],1,i,j);
				//			  vcf = vcf + (1.0-vcf)*ResidCovrg;
				// new in 8.02.04
				vcf = vcf + ResidCovrg;
				if (vcf>1.0) // limiting vcf to 1.0
				{
					vcf=1.0;
				}

				int s=(int)SOIL[i][j];

				// distribution of sand, clay and silt for stationary infiltration (should be k_sat?)
				SOILtable[s].getGrainSizeDist(dFractionSand,dFractionSilt,dFractionClay,dFractionStones1,dFractionStones2,dFractionStones3,dFractionStones4);

				// old method double dg = ((ln004+ln2)*dFractionClay* + (ln2+ln63)*dFractionSilt + (ln63+ln2000)*dFractionSand)/2.0;
				//				dg = pow(double(10.0),dg);
				//				double inf_infinite = dg/6.0;
				// new 20081016
				double dg = (FClay*dFractionClay + 
				FSilt*dFractionSilt + 
				FSand*dFractionSand + 
				FStones1*dFractionStones1 + 
				FStones2*dFractionStones2 + 
				FStones3*dFractionStones3 + 
				FStones4*dFractionStones4)/2.0;

				dg = pow(double(10.0),dg);

				// time since last soil tillage TSBB
				int JD=LUtable[luse].GetJulDay(tag,monat,jahr);
				if (iSoilTillageDays[0]>-9998)
				{
					if (JD>=iSoilTillageDays[iSoilTillageCount])
					{ // after last day of soiltillage array, so simply 
						TSBB[i][j]=JD-iSoilTillageDays[iSoilTillageCount];
						if (TSBB[i][j]<=0) TSBB[i][j]=0.3; // with 0, the pow-function is undefined
					}
					else
					{
						for (int ijd=1;ijd<=iSoilTillageCount;ijd++)
						{
							if (iSoilTillageDays[ijd]>JD)
							{
								TSBB[i][j]=JD-iSoilTillageDays[ijd-1];
								if (TSBB[i][j]<=0) TSBB[i][j]=0.3; // with 0, the pow-function is undefined
								break;
							}
						}
					}
				}
				else
				{
					TSBB[i][j]=0;
				}

				// convert precip into an intensity value
				double rain_intensity = P[i][j]/(InV/60.0); // precipitation as intensity in mm/h

				// rain pause
				RAINPAUSE[i][j]=(rain_intensity>0.05)?0:RAINPAUSE[i][j]+InV/1440.0;

				if (method==3)
				{
					varlist->add("A",SOILtable[s].get_SU_PAR(1));
					varlist->add("B",SOILtable[s].get_SU_PAR(2));
					varlist->add("C",SOILtable[s].get_SU_PAR(3));
					varlist->add("D",SOILtable[s].get_SU_PAR(4));
					varlist->add("E",SOILtable[s].get_SU_PAR(5));
					varlist->add("F",SOILtable[s].get_SU_PAR(6));
					varlist->add("G",SOILtable[s].get_SU_PAR(7));
					varlist->add("H",SOILtable[s].get_SU_PAR(8));
					varlist->add("I",SOILtable[s].get_SU_PAR(9));
					varlist->add("J",SOILtable[s].get_SU_PAR(10));
					varlist->add("K",dg);
					varlist->add("L",dFractionSand);
					varlist->add("M",dFractionClay);
					varlist->add("N",dFractionSilt);
					varlist->add("O",TSBB[i][j]);
					varlist->add("P",rain_intensity);
					varlist->add("Q",EKIN[i][j]);
					varlist->add("R",InV);
				}

				double inf_infinite=0; 
				double energy=0;
				double c1=0;
				double c2=0;
				double c3=0;
				double inf_start=0;
				double inf_pot=0;
				double inf_act=0;

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// start of dynamic expressions section /////////////////////////
////////////////////////// each block may also be calculated by /////////////////////////
////////////////////////// an expression defined in the control /////////////////////////
////////////////////////// file (when using method 3)           /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

				int actExpression=0;
				if (method==3)
				{
					// sample expressions:
					//
					// W = ((P>0.05)&(P<76.2))*(11.89+8.73*log10(Abs(P+0.001))) + (P>=76.2)*28.33; output, energy (used internally later for cumulated EKIN calculation
					// X = 65.2; output, start infiltration rate (same as parameter SU_PAR01 may be, but WaSiM needs it internally in X)
					// Y = B*K^C*(L*100)^D; output, inf_infinite
					// Z = ; output, Silting-up disposition
					// C1 = (100*L)^F; # F = SU_PAR06
					// C2 = K^G; # G = SU_PAR07
					// C3 = (O+0.001)^H; # H = SU_PAR08
					// Z = (O<=0) + (O>0)*(E*C1*C2*C3); # E = SU_PAR05
					while ((actExpression < expressionlist.size()) && (expressionlist[actExpression].find("V=")==std::string::npos))
					{
						parserresult = parser->parse(expressionlist[actExpression].c_str());
						if (std::string(parserresult).find("Error")!=std::string::npos)
						{
							(*WSstream) << "FATAL: \n-----------------error------------\n"
							"in instance "<<instanzname<<"\n"
							"expression "<<expressionlist[actExpression].c_str()<<" could not be parsed correctly. "
							<<"Parser message: "<<parserresult<<"[end-entry]";
							throw(403);
						}
						actExpression++;
					}
					// values should now be set; WaSiM will use the custom values to initialize internal variables
					double retval;
					varlist->get_value("Z",&retval);
					SDISP[i][j]=(double)retval;
					varlist->get_value("Y",&inf_infinite);
					varlist->get_value("X",&inf_start);
					varlist->get_value("W",&energy);
				}
				else
				{
					// end infiltration rate
					// Y = B*K^C*(L*100)^D
					inf_infinite = parB*pow(double(dg),double(parC))*pow(double(dFractionSand*100.0),double(parD));

					// actual cinetic energy
					// W = ((P>0.05)&(P<76.2))*(11.89+8.73*log10(Abs(P+0.001))) + (P>=76.2)*28.33
					if ((rain_intensity >0.05)&&(rain_intensity <76.2))
					{
						energy=(11.89+8.73*log10(rain_intensity));
					}
					else
					{
						if (rain_intensity>=76.2)
						{
							// old: energy=28.33*P[i][j];
							energy=28.33;
						}
						else
						{
							energy=0;
						}
					}
					
					// silting up disposition
					// C1 = (100*L)^F; # F = SU_PAR06
					// C2 = K^G; # G = SU_PAR07
					// C3 = (O+0.001)^H; # H = SU_PAR08
					// Z = (O<=0) + (O>0)*(E*C1*C2*C3); # E = SU_PAR05
					if (TSBB[i][j]>0)
					{
						c1=pow(double(100.0*dFractionSand),double(parF));
						c2=pow(double(dg),double(parG));
						c3=pow(double(TSBB[i][j]),double(parH));
						SDISP[i][j]= parE*c1*c2*c3;
					}
					else
					{
						SDISP[i][j]= 0.0;
					}
					
					// starting infiltration rate (quite high!)
					// X = 65.2
					inf_start=parA;
				}


				// cumulative energy for an event
				if (TSBB[i][j]>0) // means: iSoilTillageDays[0]>-9998
				{
					if (RAINPAUSE[i][j]<RainPause[zarray[(int)zn[i][j]]])
					{
						EKIN[i][j]=energy*rain_intensity*(1.0-vcf)*InV/60.0+EKIN[i][j];
					}
					else
					{
						EKIN[i][j]=energy*rain_intensity*(1.0-vcf)*InV/60.0; // will be 0 if P=0 since energy was set to 0
					}
				}
				else
				{
					EKIN[i][j]=0;
				}

				if (method==3)
				{
					varlist->set_value("Q",EKIN[i][j]);
					// sample expressions:
					// V -> output, potential infiltration (may be limited by actual precipitation
					// uses the previously calculated values for 
					// X = inf_start
					// Y = inf_infinite
					// Z = SDISP
					// Q = EKIN (using W = energy)
					// V = ((X-Y)*exp(-Z*Q)+Y)*R/60
					while (actExpression < expressionlist.size())
					{
						parserresult = parser->parse(expressionlist[actExpression].c_str());
						if (std::string(parserresult).find("Error")!=std::string::npos)
						{
							(*WSstream) << "FATAL: \n-----------------error------------\n"
							"in instance "<<instanzname<<"\n"
							"expression "<<expressionlist[actExpression].c_str()<<" could not be parsed correctly. "
							<<"Parser message: "<<parserresult<<"[end-entry]";
							throw(403);
						}
						actExpression++;
					}
					varlist->get_value("V",&inf_pot);
				}
				else
				{
					inf_pot = ((inf_start-inf_infinite)*exp(-SDISP[i][j]*EKIN[i][j])+inf_infinite)*InV/60.0;
				}

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  end of dynamic expressions section  /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


				// the existing water in the surface storage for silting up must be added to the new precipitation input
				P[i][j]+=SSSU[i][j];
				SSSU[i][j]=P[i][j]; // the total content in storage for silting up is now equal to the precipitation (especially for later tracer balancing)

				if (inf_pot>P[i][j])
				{
					inf_act=P[i][j];
				}
				else
				{
					inf_act=inf_pot;
				}

				// direct runoff from silting up modul -> is stored and transfered to unsatzon model and later on to surface routing and coming back from there in the next time step
				if (iSoilTillageDays[0]>-9998)
				{
					QDSU[i][j]=P[i][j]-inf_act; 
				}
				else
				{
					QDSU[i][j]=0;
				}

				P[i][j]-=QDSU[i][j]; // for input into soil model the remaining precipitation must be calculated
				// balancing tracer concentrations for output-flux
				for (int z=1;z<=tracer_count;z++)
				{
					if (QDSU[i][j]>0)
					{
						if (tracer_info[z].evaporating)
						tracer_module[z].concflux[i][j]=tracer_module[z].concstor[i][j];
						else
						{   // subtracting massflux from storage and add it to the outflow 
							double trms    = (QDSU[i][j]+SSSU[i][j])>0.0 ? tracer_module[z].concstor[i][j]*SSSU[i][j]/(QDSU[i][j]+SSSU[i][j]) : 0.0;
							double trmsmax = tracer_info[z].max_conc*SSSU[i][j];
							if (trmsmax<trms) trms=trmsmax;
							tracer_module[z].concstor[i][j]-=trms;
							tracer_module[z].concflux[i][j]=(QDSU[i][j]>0)?trms/QDSU[i][j]:0;
						}
					}
				}
			}
			if (method==3)
			{
				varlist->clear();
			}
		}
		if (varlist) { delete varlist; varlist=NULL;}
		if (parser) { delete parser; parser=NULL; }
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//																											   
	// end of processing the central loop
	//																											   
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// statistics must be done outside the parallel section because common variables will be accessed by wertzufuegen
	// as well as by tracer_module[z].stat->zone[].sx and other operators
	for (int i = mpistartrow; i <= mpiendrow; i++)
	{
		for (int j=FirstLastCol[i].First;j<=FirstLastCol[i].Last;j++) 
		{
			if ((zn[i][j]!=nodata)&&(MPIREG[i][j]==mpi_rank))
			//    if (zn[i][j]!=nodata)
			{
				if (doSpinUp!=0 && HRUMaster[i][j]==nodata)
				{ 
					// getting the coordinates of the master cell
					int imaster = vHRU->at(int(HRU[i][j]))->MasterRow;
					int jmaster = vHRU->at(int(HRU[i][j]))->MasterCol;
					// filling the grid cells
					SSSU[i][j]		= SSSU[imaster][jmaster];
					TSBB[i][j]		= TSBB[imaster][jmaster];
					RAINPAUSE[i][j]	= RAINPAUSE[imaster][jmaster];
					EKIN[i][j]		= EKIN[imaster][jmaster];
					SDISP[i][j]		= SDISP[imaster][jmaster];
					P[i][j]			= P[imaster][jmaster];
					QDSU[i][j]		= QDSU[imaster][jmaster];
				}

				int TG=zarray[(int)zn[i][j]];
				SDISPstat->wertzufuegen(i,j,zarray[(int)zn[i][j]],SDISP[i][j]);
				QDSUstat->wertzufuegen(i,j,zarray[(int)zn[i][j]],QDSU[i][j]);

				// since the tracer concentration is not changed in this module (no mixing takes place with any other flux)
				// the resulting concentrations will be equivalent to the input
				for (z=1;z<=tracer_count;z++)
				{
					if (RelCellArea!=NULL)
					{
						tracer_module[z].stat->zone[TG].sx+=tracer_module[z].concflux[i][j]*RelCellArea[i][j]*(QDSU[i][j]);
						tracer_module[z].stat->zone[TG].n+=(QDSU[i][j]);
						tracer_module[z].stat->zone[TG].RelCellSum+=(QDSU[i][j])*RelCellArea[i][j];
					}
					else
					{
						tracer_module[z].stat->zone[TG].sx+=tracer_module[z].concflux[i][j]*(QDSU[i][j]);
						tracer_module[z].stat->zone[TG].n+=(QDSU[i][j]);
						tracer_module[z].stat->zone[TG].RelCellSum+=(QDSU[i][j]);
					}
				}
			}
		}
	}

	SDISPstat->calc_and_done(jahr,monat,tag,stunde, minute);
	QDSUstat->calc_and_done(jahr,monat,tag,stunde, minute);

	// writing result grids
	writegridtofile(sdispindex);
	writegridtofile(rainpauseindex);
	writegridtofile(ekinindex);
	writegridtofile(tsbbindex);
	writegridtofile(qdsuindex);

	// writing tracer result grids
	for (z=1;z<=tracer_count;z++)
	{
		tracer_module[z].stat->calc_and_done(jahr,monat,tag,stunde, minute);
		writegridtofile(tracer_module[z].fluxindex);
	}
}


/* 3.3.97: neu anstelle des Destructors wird der Speicher extra freigegeben.
C++ soll Standarddestructor erzeugen, da bei Borland C++ 4.5 oft beim
Freigeben (Aufruf des Destructors) unter Win32 Fehler in CW3215.dll erzeugt wurden. */
void siltingupclass::ClearMemory()
{

	if (do_it > 0)
	{
		// free allocated memeory

		if (RainPause) delete[] RainPause; RainPause=NULL;
		//    if (ParAB delete[] ParB; ParB=NULL;

		SDISPstat->ClearMemory();
		QDSUstat->ClearMemory();

		if (SDISPstat) delete(SDISPstat);
		if (QDSUstat) delete(QDSUstat);

		if (zoneSDISP) free(zoneSDISP); zoneSDISP=NULL;
		if (zoneQDSU)  free(zoneQDSU);  zoneQDSU=NULL;

		free_matrix(gridinfo[sdispindex].pointer, gridinfo[sdispindex].nrl, gridinfo[sdispindex].nrh, 1, ncols,nrows);
		free_matrix(gridinfo[rainpauseindex].pointer, gridinfo[rainpauseindex].nrl, gridinfo[rainpauseindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[qdsuindex].pointer, gridinfo[qdsuindex].nrl, gridinfo[qdsuindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[ekinindex].pointer, gridinfo[ekinindex].nrl, gridinfo[ekinindex].nrh, 1, ncols, nrows);
		free_matrix(gridinfo[tsbbindex].pointer, gridinfo[tsbbindex].nrl, gridinfo[tsbbindex].nrh, 1, ncols, nrows);

		delete_grid(sdispindex);
		delete_grid(rainpauseindex);
		delete_grid(qdsuindex);
		delete_grid(ekinindex);
		delete_grid(tsbbindex);

		if (SDISPstatfilename) free(SDISPstatfilename); SDISPstatfilename=NULL;
		if (QDSUstatfilename) free(QDSUstatfilename); QDSUstatfilename=NULL;

		if (LUSE) delete[] LUSE; LUSE=NULL;

		this->expressionlist.clear();

	}
	modelobject::ClearMemory();

}

void siltingupclass::substance_transport_ini(int &retvalue,int p_tracer_count, tracer_info_type *p_tracer_info_table,TLastResults *LastResults,int maxLastResults)
{
	int z;
	char *tempstr1,*tempstr2,*tempstr3,*tempstr4;
	if (do_it==0) return;
	else
	{
		teilmodellklasse::substance_transport_ini(retvalue,p_tracer_count,p_tracer_info_table,LastResults,maxLastResults);
		if (retvalue!=0) throw retvalue;
		if (tracer_count<=0) return;
		else
		{
			tempstr1=newstring(200,sizeof(char));
			tempstr2=newstring(200,sizeof(char));
			tempstr3=newstring(200,sizeof(char));
			tempstr4=newstring(200,sizeof(char));

			for (z=1;z<=tracer_count;z++)
			{
				/* if the tracer is evaporating, the input concentration is always the same as in the last storage,
		but if it is not evaporating (salts), the input-concentration may not rise above the maximum
		concentration, therefore another grid containing the concentration in the flux is allocated and
		consequently searched for here in such cases */

				sprintf(tempstr4,"concflux_tracer_%1d_%s",z,"input"); 
				sprintf(tempstr3,"concflux_tracer_%1d_%s",z,_irrigation_parameters_); 
				sprintf(tempstr2,"concflux_tracer_%1d_%s",z,_Schneemodell_);
				sprintf(tempstr1,"concflux_tracer_%1d_%s",z,_Interzeptionsmodell_);
				tracer_module[z].cinput=gridinfo[getindex_first_of_4(tempstr1,tempstr2,tempstr3,tempstr4)].pointer;
			}
			if (tempstr1) free(tempstr1); tempstr1=NULL;
			if (tempstr2) free(tempstr2); tempstr2=NULL;
			if (tempstr3) free(tempstr3); tempstr3=NULL;
			if (tempstr4) free(tempstr4); tempstr4=NULL;
			return;
		}
	}
}

void siltingupclass::substance_transport_delete()
{
	if (do_it==0) return;
	else
	teilmodellklasse::substance_transport_delete();
}


void siltingupclass::flushbuffers()
{
	int z;
	if (do_it==0) return;
	else
	{
		if (SDISPstat) SDISPstat->flushbuffer();
		if (QDSUstat) QDSUstat->flushbuffer();
		for (z = 1; z <= tracer_count; z++)
		{
			tracer_module[z].stat->flushbuffer();
		}
	}
}

/* ++++++++++++++++++++ end of silting up routines ++++++++++++++++++++++++++++++ */

