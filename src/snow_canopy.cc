/* Canopy_Interception.f90
   Modification of meteorological parameters for a forest canopy
   Second try: according to: Liston & Elder, JHM 2005, und diese v.a. nach
   Pomeroy et al., HP 1998, und Hellstroem, HP 2000.

!Meaning of the variables in G.L.'s paper and coded hereafter in this hack (see Canopy.fi):
! Qcs  =  sublimation of snow held within the canopy (kg/m**2)                                     = SUBLIMSNOW
! PSIs =  sublimation-loss rate coefficient for an ice sphere (1/s)                                = LOSSRATECOEFF
! dm/dt=  rate of mass loss from an ice sphere                                                     = MASSLOSSRATE
! I    =  intercepted canopy load (kg/m**2)                                                        = CANOPINTLOAD
! Iold =  intercepted canopy load in the previous time step (kg/m**2)                              = CANOPINTLOADOLD
! Ce   =  canopy exposure coefficient                                                              = CANEXPCOEFF
! Imax =  maximum interception storage (kg/m**2)                                                   = MAXINTERCEPSTOR
! P    =  snow precipitation in the current time step (kg/m**2)                                    = SNOWPRECIP
! LAI  =  effective LAI (includes stems, leaves, branches) (m**2/m**2)                             = LAI
! kc   =  dimensionless coefficient = 0.010 related to the shape of the intercepted snow deposits  = DIMLESSCOEFF
! m    =  particle mass (kg)                                                                       = PARTICLEMASS
! r    =  radius of a sperical ice particle (500 micrometer)                                       = ICERADIUS
! rhoi =  ice density (kg/m**2)                                                                    = ICEDENSITY
!                                                                                                    OMEGA
! M    =  molecular weight of water = 18.01 (kg/mole)                                              = MOLWEIGHTWAT
! R    =  universal gas constant = 8313 (J/(mole K))                                               = UNIVGASCONST
! lt   =  thermal conductivity of the atmosphere = 0.024 (J/(m s K))                               = THERMCONDATMO
! D    =  diffusivity of water vapour in the atmosphere (m**2/s)                                   = WATVAPDIFF
! rhov =  saturation density of water vapour (kg/m**3)                                             = WATVAPSATDENS
! Rd   =  gas constant of dry air = 287 (J/(deg kg))                                               = GASCONSTDRYAIR
! es   =  saturation vapor pressure over ice (Pa)                                                  = SATVAPPRESSICE
! Nu   =  Nusselt number                                                                           = NUSSELT
! Sh   =  Sherwood number                                                                          = SHERWOOD
! Re   =  Reynolds number                                                                          = REYNOLDS
! v    =  kinematic viscosity of air = 1.3*10-5 (m**2/s)                                           = KINVISCOSAIR
! uc   =  ventilation velocity = wind speed within the canopy                                      = WIND
! h    =  canopy height (m)                                                                        = CANHEIGHT
! z    =  canopy wind speed reference level = 0.6*h (m)                                            = CANWINDREFLEV
! a    =  canopy flow index                                                                        = CANFLOWIND
! beta =  scaling factor = 0.9                                                                     = SCALEFAC
! Sp   =  solar radiation absorbed by the snow particle (W)                                        = ABSORBRAD
! Lm   =  melt-unloading rate per time step (kg/m**2)                                              = MELTUNLOAD
! Qsif =  solar radiation reaching the snow surface underneath the canopy (W/m=)                   = modified DIR & DIF
! taun =  fraction of incoming solar radiation transmitted through the canopy                      = CANOPTRANSM
! k    =  vegetation-dependent extinction coefficient                                              = EXTINCTCOEFF
*/

#include "snow_canopy.h"
#include <math.h>

using namespace std;
using namespace snow_canopy;


//!Berechnung der Albedo
//!Uli Strasser, Maerz 2004, LMU
bool snow_canopy::SIM_ALBEDO(double ALBMIN, double ALBADD, double ABKLINGPOS, 
	double ABKLINGNEG, double SNOWPRECIP,	double &ALBEDO, 
	double &SNOWAGE,double AIRTEMP, double INTERVALL)
{

	if((SNOWPRECIP/INTERVALL) > 0.0001389)        //Bei mehr als 0.5 mm/h Schnee
	{
		ALBEDO  = ALBMIN + ALBADD;                //Albedo auf Startwert setzen
		SNOWAGE = 0.0;                             //Alter re-initialisieren!
	}


	if(ALBEDO > ALBMIN && SNOWAGE > 0)
	{
		if(AIRTEMP >= 273.16)
			ALBEDO = ALBMIN + (ALBEDO - ALBMIN) * exp(-ABKLINGPOS * INTERVALL/86400.);
		else
			ALBEDO = ALBMIN + (ALBEDO - ALBMIN) * exp(-ABKLINGNEG * INTERVALL/86400.);
	}

	// Update snowage
	SNOWAGE  = SNOWAGE + (INTERVALL/86400.);
	return true;
}

//!Berechnung der Schneeschmelze mit der T-Index Methode
//!Uli Strasser, Maerz 2006, LMU
double snow_canopy::T_INDEX(double AIRTEMP, double ALBEDO, double DIR,
		double DIF, int INTERVALL, double TI_temp, double TI_rad)
{
	double temperatureFactor, albedoFactor, TINDEXMELT;

	//********************************************************************
	//!Nach G. Liston
	//TINDEXMELT = (5.0/86400.0) * (AIRTEMP - 273.16) * INTERVALL;
	//********************************************************************

	//********************************************************************
	//Nach F.Pellicciotti
	temperatureFactor = TI_temp; // 0.05;   // !mm/(h C)
	albedoFactor      = TI_rad;  // 0.0094; // !m**2 mm/(Wh)
	TINDEXMELT = (temperatureFactor/(60.*60.) * (AIRTEMP - 273.16) +
				   albedoFactor/(60.*60.) * (1. - ALBEDO) * (DIR + DIF)) * INTERVALL;
	//********************************************************************
	return TINDEXMELT;
}

double snow_canopy::max2(double a, double b)
{
	if(a>b) return a;
	else    return b;
}

double snow_canopy::min2(double a, double b)
{
	if(a>b) return b;
	else    return a;
}

//!Berechnung des Dampfdruckes mit der Magnus-Formel
//!Uli Strasser, Juli 2004, LMU
double snow_canopy::VAP_PRESSURE(double Temperature, double Humidity)
{
	// MAGNUS-Formel @todo !!!!!!!!!!!
	//***********************************************************************************
	// double PRESSURE = 6.1078 * exp((17.08085*(Temperature-273.16))/(234.175+(Temperature-273.16)))*(Humidity/100.);
	//***********************************************************************************
	//return PRESSURE;
	
	// Use the formulae provided by
	// Alduchov, O. and Eskridge, R.: Improved Magnus� form approximation of saturation vapor pressure, Tech. rep., 
	// Department of Commerce, Asheville, NC (United States), doi:10.2172/548871, http://dx.doi.org/10.2172/548871, 1997.
	double vap_press;
	if (Temperature >= 273.15)
		vap_press = 6.112 * exp((17.62 * (Temperature - 273.15)) / (243.12 + (Temperature - 273.15))) * Humidity / 100.;
	else
		vap_press = 6.112 * exp((22.46 * (Temperature - 273.15)) / (272.62  + (Temperature - 273.15))) * Humidity / 100.;
	return vap_press;
}

//!Berechnung der Feuchttemperatur nach Wasti Weber
//!Uli Strasser, Mai 2007, LMU
// SUBROUTINE WET_TEMPERATURE
double snow_canopy::WETBULB_TEMPERATURE(double AIRTEMP, double HUM, double ALTITUDE)
{
	double WATVAPPRESSDIFF, PSYCHKONST, AIRPRESSURE;
	int x;
	//!********************************************************************
	double WETTEMP  = AIRTEMP;
	AIRPRESSURE     = 1013.25 * exp(-0.0001184 * ALTITUDE);
	PSYCHKONST      = 0.67 * (AIRPRESSURE/1007.0);
	WATVAPPRESSDIFF = VAP_PRESSURE(AIRTEMP, HUM) -
		(VAP_PRESSURE(WETTEMP, 100.) - PSYCHKONST * (AIRTEMP - WETTEMP));
	//!********************************************************************
	x = 0;
    do{
		x++;
		if(x==150) break;
		WETTEMP         = WETTEMP - 0.1;
		WATVAPPRESSDIFF = VAP_PRESSURE(AIRTEMP, HUM) -
			(VAP_PRESSURE(WETTEMP, 100.) - PSYCHKONST * (AIRTEMP - WETTEMP));
	}while(WATVAPPRESSDIFF < 0.0);
	//********************************************************************
	if(x==150) WETTEMP = AIRTEMP;
	//********************************************************************
	return WETTEMP;
}

bool snow_canopy::PhasePartitioning(double AIRTEMP, double HUM, double PRECIP,
	double ALTITUDE, double &SNOWPRECIP, double &RAINPRECIP)
{
	//****************************************************************************
	// Phasenuebergang mit der Feuchttemperatur
	// CALL WET_TEMPERATURE
	double WETTEMP = WETBULB_TEMPERATURE(AIRTEMP, HUM, ALTITUDE);
	// double SNOWPRECIP, RAINPRECIP;
	// REAL :: PHASE //!0=Regen, 1=Schnee
	double PHASE = 0;
	if(PRECIP > 0.0)
		if(WETTEMP < 275.16 && AIRTEMP < 276.16)	// Steinacker (1980) ?
			PHASE = 1.;
	if(PHASE==1.)
	{
		SNOWPRECIP = PRECIP;
		RAINPRECIP = 0.0;
	}
	else
	{
		RAINPRECIP = PRECIP;
		SNOWPRECIP = 0.0;
	}
	return true;
}

// $kf 2018-06-15: from WaSiM snow modell
// this is basically a copy from snowmod.cc
// @todo: move to separate function that we can call from both snow and interception model
void snow_canopy::PhasePartitioningWaSiM(double temperatur, double T0R, double transzone,
	double precip, double &SNOWPRECIP, double &RAINPRECIP)
{
	double schneeanteil;
	if (temperatur <= (T0R+transzone))
	{
		// Lufttemperatur unter der Grenztemperatur T0R 
		// Ermitteln des Schneeanteils am Niederschlag 
		if (temperatur>(T0R-transzone))
		{
			if (temperatur < 0 && precip > 0.)
			{
				int z=0;
			}
			schneeanteil=((T0R+transzone)-temperatur)/(transzone+transzone);
			
		}
		else
			schneeanteil = 1.0;
	}
	else
		schneeanteil = 0.0;
	SNOWPRECIP = schneeanteil*precip;
	RAINPRECIP = (1-schneeanteil)*precip;
}

//Berechnung der atmosphaerischen Gegenstrahlung
//Uli Strasser, Juli 2004, LMU
//Erweiterung durch Thomas Marke fuer Methodenunterscheidung und Berechnung der
//Bewoelkung der langwelligen Einstrahlung ueber die Feuchte), Januar 2012, KFU Graz
double snow_canopy::COUNTERRAD(double STATIONELEVATION, double AIRTEMP, double HUM, double PRECIP, double NB){
	int NB_PROVIDED;      // !Schalter fuer die Bereitstellung der Bewoelkung
	/*double TEMPLR;        // !Temperatur Laps Rate (�C/km)
	double DEWTEMPLR;     // !Taupunkt Temperatur Laps Rate (�C/km)
	double DEWTEMP;       // !Taupunkt Temperatur (K)
	double DEWTEMP3000;   // !Taupunkt Temperatur 3000m (K)
	double TEMP3000;      // !Lufttemperatur 3000m (K)
	double RH3000;        // !Relative Luftfeuchte 3000m (K)
	double PRESAIR3000;   // !Dampfdruck in 3000m (hPa)
	double SATPRESAIR3000;// !Saettigungsdampfdruck in 3000m (hPa)
	double ELEDIFF;       // !Differenz zwischen Stationshoehe und 3000m (m)*/

	//Aktuellen Dampfdruck berechnen
	double PRESAIR  = VAP_PRESSURE(AIRTEMP,HUM);
	//!Schalter fuer die Bereitstellung/Berechnung der Bewoelkung setzen (muss derzeit noch hier im Code gesezt werden)
	//NB_PROVIDED = 0; //(1 = bereitgestellte Bewoelkung verwenden, 0 = Bewoelkung berechnen nach Liston & Elder 2006 (Micromet))

/*    double NB = 0;
	if(NB_PROVIDED == 0)
	{
		TEMPLR = 0.00378;    // !Temperatur Lapse Rate �C/m aus Daten fuer die Obere Donau (1971-2000) f�r den Zeitraum Nov.-Apr.
		DEWTEMPLR = 0.00467; // !Taupunkt Temperatur Lapse Rate �C/m aus Daten f�r die Obere Donau (1971-2000) f�r den Zeitraum Nov.-Apr.

		//! Taupunkttemperatur berechnen
		// $kf 2012-08-08 ALOG durch log aus <math.h> ersetzt
		if(AIRTEMP > 273.15)
			DEWTEMP = 240.97 * log(PRESAIR / 6.1078) / (17.08085 - log(PRESAIR / 6.1078)) +272.15;
		else
			DEWTEMP = 272.55 * log(PRESAIR / 6.1071) / (22.4429 - log(PRESAIR / 6.1078)) +272.15;

		//!Unterschied Stationshoehe-3000m in m berechnen
		ELEDIFF = (3000.0 - STATIONELEVATION);

		//!Taupunkttemperatur in 3000m berechnen
		DEWTEMP3000 = DEWTEMP - (DEWTEMPLR * ELEDIFF);

		// !Temperatur in 3000m berechnen
		TEMP3000 =  AIRTEMP - (TEMPLR * ELEDIFF);

		//Aktuellen Dampfdruck und Saettigungsdampfdruck in 3000m berechnen
		SATPRESAIR3000 = VAP_PRESSURE(TEMP3000, 100.);
		PRESAIR3000 = VAP_PRESSURE(DEWTEMP3000, 100.);

		//!Relative Feuchte berechnen
		RH3000 = 100.0 * PRESAIR3000/SATPRESAIR3000;

		//!Bewoeklung berechnen
		NB = 0.832 * exp((RH3000 - 100.0) / 41.6);
	}*/
	//*************************************************************************************************
	//Nach Brunt (1932), wolkenloser Fall..
	double at = AIRTEMP*AIRTEMP;
	at = at*at; // airtemp^4
	//js
	double IR = STEFBOLTZ * at * (0.610 + 0.050*sqrt(PRESAIR));
	// original version uses two pow-calls, which is usually much slower for exponents > 2 and also slower than sqrt()
	//	double IR = STEFBOLTZ*pow(AIRTEMP,4.)*(0.610+0.050*pow(PRESAIR,0.5));
	
	//Nach Bolz (1949), Modifikation durch Wolken
	double BOLZ;
	if(PRECIP == 0.0)
	{
		BOLZ = 0.17;				  //Koeffizienten nach Escher-Vetter (2001)!
		if(NB <= 0.2) BOLZ = 0.08;
		if(NB >= 0.6) BOLZ = 0.20;
	}
	else
		BOLZ = 0.24;
	//js
	IR = IR * (1 + BOLZ* NB * NB);
//	IR = IR * (1 + BOLZ*pow(NB, 2.));
	return IR;
}

bool snow_canopy::CANOPY_INTERCEPTION(SnowCanopyParam* p)
{
//******************************************************************************
//!1.) Solar radiation - same for snowfall, rainfall or no precipitation
//Solar radiation absorbed by a snow particle in the canopy
	// CALL SIM_ALBEDO
	// SIM_ALBEDO via call-by-reference
//	SIM_ALBEDO(p, index, h->SNOWPRECIP, h->ALBEDO, h->SNOWAGE, h->AIRTEMP);

	// workaround (as we only have the global radiation):
	double DIF = 0.;
	double DIR = p->SHORTWAVE;

    // Albedo computation for the snow layer on the trees
	SIM_ALBEDO(p->ALBMIN, p->ALBADD, p->ABKLINGPOS, p->ABKLINGNEG, p->SNOWPRECIP, 
		p->ALBEDO, p->SNOWAGE, p->AIRTEMP, p->INTERVALL);
	//js
	double ABSORBRAD = 3.14159 * ICERADIUS * ICERADIUS * (1. - p->ALBEDO) * (DIR + DIF);
//	double ABSORBRAD = 3.14159 * pow(ICERADIUS, 2.) * (1. - p->ALBEDO) * (DIR + DIF);

	//Solar radiation reaching the ground beneath the canopy
	double CANOPTRANSM = exp(-EXTINCTCOEFF * p->LAI);
	DIR = CANOPTRANSM * DIR;	// modifizierte Strahlung speichern
	DIF = CANOPTRANSM * DIF;
	p->SHORTWAVE = DIR + DIF; // return shortwave radiation

//!2.) Thermal radiation
	double CANOPYFRACTION = 0.55 + 0.29 * log(p->LAI);
	if(CANOPYFRACTION < 0.0) CANOPYFRACTION = 0.0;
	if(CANOPYFRACTION > 1.0) CANOPYFRACTION = 1.0;
	// call counter radiation function
	p->IR = COUNTERRAD(p->ALTITUDE, p->AIRTEMP, p->HUM, p->PRECIP, p->CL);
	//mk: modify counter radiation in the canopy -> overwritten in CEB
        //assumption: Canopy temp = Air temp !!
	//js
	double at = p->AIRTEMP*p->AIRTEMP;
	at = at*at;
	p->IR = (1. - CANOPYFRACTION) * p->IR + CANOPYFRACTION * STEFBOLTZ * at;
//	p->IR = (1. - CANOPYFRACTION) * p->IR + CANOPYFRACTION * STEFBOLTZ * pow(p->AIRTEMP, 4.);

//******************************************************************************

//******************************************************************************
//!2.) Wind speed in the canopy - same for snowfall, rainfall or no precipitation
	//Ventilation velocity = canopy wind speed
    p->WIND_OLD = p-> WIND;
	double CANFLOWIND = SCALEFAC * p->LAI;               // Canopy flow index
	double CANWINDREFLEV = 0.6 * p->CANHEIGHT;           // Reference level
	p->WIND = p->WIND * exp(-CANFLOWIND * (1. - CANWINDREFLEV/p->CANHEIGHT));
//******************************************************************************
	//Reynolds, Nusselt and Sherwood numbers for the mass loss rate description
	double REYNOLDS = 2. * ICERADIUS * p->WIND / KINVISCOSAIR;
	if(REYNOLDS < 0.7) REYNOLDS = 0.7;
	if(REYNOLDS > 10.) REYNOLDS = 10.;
	//js
	double NUSSELT = 1.79 + 0.606 * sqrt(REYNOLDS);
	double SHERWOOD = 1.79 + 0.606 * sqrt(REYNOLDS);
//	double NUSSELT = 1.79 + 0.606 * pow(REYNOLDS, 0.5);
//	double SHERWOOD = 1.79 + 0.606 * pow(REYNOLDS, 0.5);
	//******************************************************************************

//******************************************************************************
//!3.) Humidity in the canopy - to be included
        p->HUM_OLD = p->HUM;
	p->HUM = p->HUM * (1 + 0.1 * CANOPYFRACTION);
	if(p->HUM > 100.) p->HUM = 100.;
//******************************************************************************

//******************************************************************************
//!4.) Temperature, modified by the canopy
        p->AIRTEMP_OLD = p->AIRTEMP; //mk: store original air temperature
        //p->AIRTEMP = p->CANTEMP_OLD; //mk: big leaf temperature
	    p->AIRTEMP = p->AIRTEMP - CANOPYFRACTION * (p->AIRTEMP -(0.8 * (p->AIRTEMP - p->MEANTEMP) + p->MEANTEMP - p->DELTATEMP));
//******************************************************************************

//******************************************************************************
/*!5.) Precipitation: snow or rain
!** Already done in SnoMod.f90
!IF (PHASE .EQ. 1) THEN
!  SNOWPRECIP = PRECIP
!  RAINPRECIP = 0.0
!ELSE
!  RAINPRECIP = PRECIP
!  SNOWPRECIP = 0.0
!ENDIF
*/
// Phase partitioning has already been computed.
//******************************************************************************

//******************************************************************************
//!6.) Now computation of interception, sublimation and melt unload
	//Saturation vapour pressure over ice
	double SATVAPPRESSICE = 611.15 * exp(22.452 * (p->AIRTEMP - 273.16) / (p->AIRTEMP - 0.61));

	//Saturation density of water vapour
	double WATVAPSATDENS = 0.622 * (SATVAPPRESSICE / (GASCONSTDRYAIR * p->AIRTEMP));

	//Diffusivity of water vapour in the atmosphere
	double WATVAPDIFF = 0.0000206 * pow((p->AIRTEMP / 273.),1.75);

	//Omega
	double OMEGA = ((CONSUBSNOW * MOLWEIGHTWAT) / (UNIVGASCONST * p->AIRTEMP) - 1.)
			/ (THERMCONDATMO * p->AIRTEMP * NUSSELT);

	//Mass loss rate
	double MASSLOSSRATE = (2. * 3.14159 * ICERADIUS * (p->HUM/100. - 1.) - ABSORBRAD * OMEGA) /
			(CONSUBSNOW * OMEGA + (1./(WATVAPDIFF * WATVAPSATDENS * SHERWOOD)));

	//Particle mass
	//js
	double PARTICLEMASS = 4. / 3. * 3.14159 * ICEDENSITY * ICERADIUS* ICERADIUS* ICERADIUS;
//	double PARTICLEMASS = 4. / 3. * 3.14159 * ICEDENSITY * pow(ICERADIUS, 3.);

	//Sublimation loss rate coefficient
	double LOSSRATECOEFF = MASSLOSSRATE / PARTICLEMASS;

	//Maximum interception storage --> REFEREMCE??
	double MAXINTERCEPSTOR = 4.4 * p->LAI; //orig: 4.4

	//Canopy-intercepted load
	p->CANOPINTLOAD = p->CANOPINTLOADOLD + 0.7 * (MAXINTERCEPSTOR - p->CANOPINTLOADOLD) *
		(1. - exp((-p->SNOWPRECIP)/MAXINTERCEPSTOR));

	//Snowfall which has intercepted in this time step
	double INTERCEPSNOWPRECIP = p->CANOPINTLOAD - p->CANOPINTLOADOLD;

	//Sublimation rate - only if intercepted snow is available - is a negative rate
	if(p->CANOPINTLOAD > 0.0)
	{
		double CANEXPCOEFF  = DIMLESSCOEFF * pow((p->CANOPINTLOAD/MAXINTERCEPSTOR),(-0.4));            //Canopy exposure coefficient
		p->SUBLIMSNOW   = min2(-(CANEXPCOEFF * p->CANOPINTLOAD * LOSSRATECOEFF * p->INTERVALL),p->CANOPINTLOAD);
		p->CANOPINTLOAD = p->CANOPINTLOAD - p->SUBLIMSNOW;
	// !IF (CANOPINTLOAD < 0.0) STOP
		// workaround: unbedingt noch pruefen (**)
		if(p->CANOPINTLOAD < 0.0)
		{
			p->SUBLIMSNOW += p->CANOPINTLOAD;
			p->CANOPINTLOAD = 0.;
        }
	}

	double POTMELTUNLOAD = 0.;
	double MELTUNLOAD = 0.;
	// p->MELTUNLOAD = 0.;	// $kf 2012-08-11: zuruecksetzen...!
	//Unloading by melt and fall down
	if(p->AIRTEMP > 273.16)
	{
		// CALL T_INDEX
		double TINDEXMELT = T_INDEX(p->AIRTEMP, p->ALBEDO, DIR,
				DIF, p->INTERVALL, p->tIndexTemp, p->tIndexRad);

		// Add components melt + falldown (calibrated with Liston & Elder 2005)
		POTMELTUNLOAD = 1.0 * TINDEXMELT + 2.3 * TINDEXMELT;
		// IF (POTMELTUNLOAD < 0.0) STOP ! // (**)
		MELTUNLOAD    = min2(POTMELTUNLOAD,p->CANOPINTLOAD);
		// Set humidity to saturation if melt occurs
		if(MELTUNLOAD > 0.0) p->HUM = 100.;
  		p->CANOPINTLOAD = p->CANOPINTLOAD - MELTUNLOAD;
	}
	
	// $kf 2016-06-20: separate liquid and solid
	p->MELTUNLOAD_L = 1.0 / 3.3 * MELTUNLOAD;
	p->MELTUNLOAD_S = 2.3 / 3.3 * MELTUNLOAD;

	// Precipitation reaching the ground ($kf: separate solid and liquid water)
	p->SNOWPRECIP = (p->SNOWPRECIP - INTERCEPSNOWPRECIP) + p->MELTUNLOAD_S; // p->MELTUNLOAD;
	p->RAINPRECIP += p->MELTUNLOAD_L;
	
	// $kf 2016-07-18: reset state variables of intercepted snow when it's gone
	if(p->CANOPINTLOAD == 0.)
	{
		p->ALBEDO = p->ALBMIN;
		p->SNOWAGE = 0.;
	}
        
	return true;
}
/*
 * CANENBAL EXTENSIONS mostly after Gouttevin et al., 2015 with further literature
 *  
 */

//////////////////////////////////////////////////////////////////////////////////////////
//CALC SHORTWAVE NET RADIATION IN CANOPY LAYER
//@SW_NET_CAN
/////////////////////////////////////////////////////////////////////////////////////////
double snow_canopy::CALC_SW_NET_CAN(double LAI, double SW_IN, double ALBEDO_CAN, double ALBEDO_SURF)
{
    double CANOPTRANSM = exp(-EXTINCTCOEFF * LAI);
    double ABSORPFACT = (1 - CANOPTRANSM);
    
    double SW_NET_CAN = SW_IN * (1 - ALBEDO_CAN) * ABSORPFACT * (1 + (ALBEDO_SURF * (1 - ABSORPFACT)/(1 - ABSORPFACT * ALBEDO_CAN * ALBEDO_SURF)));
    
    return SW_NET_CAN;  
}
/////////////////////////////////////////////////////////////////////////////////////////
// LONGWAVE NET RADIATION IN CANOPY LAYER
//@LW_NET_CAN
////////////////////////////////////////////////////////////////////////////////////////
double snow_canopy::CALC_LW_NET_CAN(double LAI, double SURFTEMP, double CANTEMP, double LW_IN)
{
    if (SURFTEMP == -9999){
        SURFTEMP = 273.16;
    }
    double CANOPTRANSM = exp(-EXTINCTCOEFF * LAI);
    double ABSORPFACT = (1 - CANOPTRANSM);
    double T_SURF = SURFTEMP * SURFTEMP * SURFTEMP * SURFTEMP;
    double T_CAN = CANTEMP * CANTEMP * CANTEMP * CANTEMP;
    
    double LW_NET_CAN = ABSORPFACT * (LW_IN + STEFBOLTZ * T_SURF - 2*STEFBOLTZ * T_CAN);
     
    return LW_NET_CAN;  
}
//////////////////////////////////////////////////////////////////////////////////////
//WIND CONDITIONS
//@WIND
////////////////////////////////////////////////////////////////////////////////////////

//CALCULATE FRICTION VELOCITY FRICTION_VELOCITY 
//NO MEASURMENTS OVER CANOPY AVAILABLE!!!! --> Calculate friction velocity from measurements, apply on canopy!!!
//@FRICTION_VELOCITY
double snow_canopy::CALC_FRICTION_VELOCITY(double WIND)
{
    double ROUGHNESS_SNOW = 0.005; //SNOW CONDITIONS! -> 0.005! WMO Guide to Meteorological Instruments and Methods of Observation WMO-No. 8 page I.5-13
    double MEASUREMENT_HEIGHT = 2.5;

    double FRICTION_VELOCITY = (WIND * KARMAN)/(log(MEASUREMENT_HEIGHT/ROUGHNESS_SNOW));

    return FRICTION_VELOCITY;
}

// ZERO DISPLACEMENT PLANE 
//@ZERO_DISPLACEMENT Between 2/3 and 3/4 of CANHEIGHT
double snow_canopy::CALC_ZERO_DISPLACEMENT(double CANHEIGHT)
{
    double ZERO_DISPLACEMENT = 0.67 * CANHEIGHT;
    
    return ZERO_DISPLACEMENT;  
}

// CALCULATE ROUGHNESS LENGTH --> DEWALLE
//@MOMENTUM_ROUGHNESS
double snow_canopy::CALC_ROUGHNESS_LENGTH(double CANHEIGHT)
{
    double ROUGHNESS_LENGTH = 0.136 * CANHEIGHT;
    
    return ROUGHNESS_LENGTH;  
}
// HEAT_VAPOR_ROUGHNESS --> Irmak 2010
//@HEAT_VAPOR_ROUGHNESS
double snow_canopy::CALC_MOMENTUM_ROUGHNESS(double CANHEIGHT)
{
    double MOMENTUM_ROUGHNESS = 0.123 * CANHEIGHT; // 0.136 -> deWalle!
        
    return MOMENTUM_ROUGHNESS;  
}


// HEAT_VAPOR_ROUGHNESS --> Irmak 2010
//@HEAT_VAPOR_ROUGHNESS
double snow_canopy::CALC_HEAT_VAPOR_ROUGHNESS(double CANHEIGHT)
{
    double MOMENTUM_ROUGHNESS = 0.136 * CANHEIGHT; // 0.136 -> deWalle!
    
    double HEAT_VAPOR_ROUGHNESS = 0.1 * MOMENTUM_ROUGHNESS;
    
    return HEAT_VAPOR_ROUGHNESS;  
}

//CALCULATE WIND SPEED AT CANHEIGHT
//Assumption: Zero displacement height + roughness length
//@WINDSPEED_ABOVE
double snow_canopy::CALC_WINDSPEED_ABOVE(double WIND, double CANHEIGHT, double Z0)
{
    //Calculate windspeed in measurement height 12 m
    double ROUGHNESS_SNOW = 0.005; //SNOW CONDITIONS! -> 0.005! WMO Guide to Meteorological Instruments and Methods of Observation WMO-No. 8 page I.5-13
    double D0 = CALC_ZERO_DISPLACEMENT(CANHEIGHT);
	double FRICTION_VELOCITY = CALC_FRICTION_VELOCITY(WIND);
    double WINDSPEED_ABOVE = FRICTION_VELOCITY/KARMAN * log((CANHEIGHT - D0)/ROUGHNESS_SNOW);

    return WINDSPEED_ABOVE;
}



//////////////////////////////////////////////////////////////////////////////////////
//RESISTANCES
//@RESISTANCES
////////////////////////////////////////////////////////////////////////////////////////

//AERODYNAMIC RESISTANCE --> Musselmann 2017
//For neutral conditions is heat resistance = momentum resistance
//@RESIST_AERO
double snow_canopy::CALC_RESIST_AERO(double CANHEIGHT, double WIND, double Z0)
{
    const double KARMAN = 0.4;
    const double z_m = 15.0; //measurement height wind speed --> CALCULATION ACCORDING TO LOGARITHMIC WIND PROFILE!
    const double z_h = 15.0; //measurement height humidity ??
    
    if (WIND == 0.0){
        WIND = 0.01;
    }
    double WIND_AC = CALC_WINDSPEED_ABOVE(WIND, CANHEIGHT, Z0);
    
    double ZERO_DISPLACEMENT = CALC_ZERO_DISPLACEMENT(CANHEIGHT);
    double MOMENTUM_ROUGHNESS = CALC_MOMENTUM_ROUGHNESS(CANHEIGHT);
    double HEAT_VAPOR_ROUGHNESS = CALC_HEAT_VAPOR_ROUGHNESS(CANHEIGHT);
    
    //OLD --> REF??
    //double RESIST_AERO = log((z_m - ZERO_DISPLACEMENT)/(MOMENTUM_ROUGHNESS)) * log((z_h - ZERO_DISPLACEMENT)/(HEAT_VAPOR_ROUGHNESS))/(KARMAN * KARMAN * WIND_AC);

    //Musselmann 2017
	double RESIST_AERO = (1./(KARMAN * KARMAN * WIND))*(log((z_m - ZERO_DISPLACEMENT)/(MOMENTUM_ROUGHNESS)))*(log((z_m - ZERO_DISPLACEMENT)/(MOMENTUM_ROUGHNESS))); //WIND_AC
    
    return RESIST_AERO;  
}

// WET CANOPY FRACTION --> Gouttevin
// @WET_CAN_FRAC
double snow_canopy::CALC_WET_CAN_FRAC(double WATINTERCEP, double MAXWATINTERCEP)
{
    
    double WET_CAN_FRAC = pow((WATINTERCEP/MAXWATINTERCEP),2/3);
    
    return WET_CAN_FRAC;  
}
// VAPOR RESISTANCE --> Gouttevin, from LU table
// RSI, RSC FROM LU Table!
//@RESIST_LE
double snow_canopy::CALC_RESIST_LE(double RSC, double RSI, double WATINTERCEP, double MAXWATINTERCEP)
{   
    //double RESIST_AERO = CALC_RESIST_AERO(CANHEIGHT, WIND, Z0);
    //double RESIST_SURF = RSC; //-> FROM LU TABLE!!!
    //double RESIST_LE = RESIST_AERO + RESIST_SURF; //LHOMME 1994
    
    double WET_FRAC = CALC_WET_CAN_FRAC(WATINTERCEP, MAXWATINTERCEP);
	//js
	double RESIST_LE = 1.0/((1 / RSI) * WET_FRAC + (1 / RSC)*(1 - WET_FRAC));
//	double RESIST_LE = pow(((1 / RSI) * WET_FRAC + (1 / RSC)*(1 - WET_FRAC)), -1);

    return RESIST_LE;  
}

//////////////////////////////////////////////////////////////////////////////////////
//BULK TRANSFER APPROACH
//@BULK TRANSFER COEFFICIENT
////////////////////////////////////////////////////////////////////////////////////////

// BULK RICHARDSON NUMBER --> DeWalle
// @RICHARDSON
double snow_canopy::CALC_RICHARDSON(double CANTEMP, double AIRTEMP, double WIND, double CANHEIGHT, double Z0)
{   
    double WIND_AC = CALC_WINDSPEED_ABOVE(WIND, CANHEIGHT, Z0);
	//wind speed over canopy = measured wind speed!!
    double RICHARDSON = (2 * 9.81)/(CANTEMP + AIRTEMP) * CANHEIGHT*((AIRTEMP - CANTEMP)/(WIND*WIND));
    
    return RICHARDSON;  
}

// BULK TRANSFER COEFFICIENT, NEUTRAL CONDITIONS --> DeWalle
// @C_HN
double snow_canopy::CALC_C_HN(double CANHEIGHT, double Z0)
{   
    double KARMAN = 0.4;
    double MEASUREMENT_HEIGHT = 15;
    double MOMENTUM_ROUGHNESS = CALC_MOMENTUM_ROUGHNESS(CANHEIGHT); //Z0
    double C_HN = KARMAN * KARMAN * pow((log(MEASUREMENT_HEIGHT/MOMENTUM_ROUGHNESS)),-2);
    
    return C_HN;  
}

// DETERMINE STABILITY CONDITIONS --> CALC BULK TRANSFER COEFFICIENTS
// @BULK_TRANS_COEFF
double snow_canopy::CALC_BULK_TRANS_COEFF(double CANTEMP, double AIRTEMP, double WINDSPEED, double CANHEIGHT, double Z0)
{   
    double RICHARDSON = CALC_RICHARDSON(CANTEMP, AIRTEMP, WINDSPEED, CANHEIGHT, Z0 );
    double TEMPDIFF = AIRTEMP - CANTEMP;
    
    double C_HN = CALC_C_HN(CANHEIGHT, Z0);
    double BULK_TRANS_COEFF = 1;
    
    
    //CONDITIONS?? --> REF!!
    if (fabs(TEMPDIFF) < 0.1){
        //neutral
        BULK_TRANS_COEFF = CALC_C_HN(CANHEIGHT, Z0);
    }
    else{
        
        if(TEMPDIFF < 0){
            //CANTEMP > AIRTEMP --> UNSTABLE
            BULK_TRANS_COEFF = C_HN * pow((1 - 16* RICHARDSON),0.75);
        }
        else{
            //CANTEMP < AIRTEMP --> STABLE
            BULK_TRANS_COEFF = C_HN * pow((1 - 5* RICHARDSON),2);
        }
    }
      
    return BULK_TRANS_COEFF;  
}


//////////////////////////////////////////////////////////////////////////////////////
//SENSIBLE HEAT FLUX
//@H_CAN
////////////////////////////////////////////////////////////////////////////////////////
double snow_canopy::CALC_SENS_AIR_CAN(double CANHEIGHT, double RHO_AIR, double HEATCAP_AIR, double CANTEMP, double AIRTEMP, double WIND, double Z0)
{
    double SENS_AIR_CAN = 0;
    double C_H = 1;
    
    if (WIND == 0.0){
        WIND = 0.01;
    }
	//CHECK WIND!!!
	double WIND_AC = CALC_WINDSPEED_ABOVE(WIND,CANHEIGHT,Z0);
    C_H = CALC_BULK_TRANS_COEFF(CANTEMP, AIRTEMP, WIND,  CANHEIGHT, Z0);

    double RESIST_AERO = CALC_RESIST_AERO(CANHEIGHT, WIND, Z0);

    //SENS_AIR_CAN = ((C_H * RHO_AIR * HEATCAP_AIR)/(RESIST_AERO_HEAT)) * (CANTEMP - AIRTEMP); //old
	//SENS_AIR_CAN = ((C_H * RHO_AIR * HEATCAP_AIR) * WIND) * (CANTEMP - AIRTEMP); // deWalle --> Kustas 1994 WIND_AC
    SENS_AIR_CAN = ((RHO_AIR * HEATCAP_AIR) /(RESIST_AERO)) * (CANTEMP - AIRTEMP ); //Gouttevin --> RESIST_AERO has no stability correction!!!

    return SENS_AIR_CAN;  
}
/////////////////////////////////////////////////////////////////////////////////////////////
//LATENT HEAT FLUX
//@LATHEAT
///////////////////////////////////////////////////////////////////////////////////////////

// CALC AIR PRESSURE WITH "Internationale Höhenformel" using a standard atmosphere (Temperatur 15 °C = 288,15 K, Luftdruck 1013,25 hPa, Temperaturgradient 0,65 K pro 100 m)
//@AIR_PRESSURE for LE formulation of deWalle
double snow_canopy::CALC_PRESSURE(double ALTITUDE, double AIRTEMP)
{
    double AIR_PRESSURE = 0;

    //Internationale Höhenformel 
    //AIR_PRESSURE = 1013.25 * pow((1 - (0.0065 * ALTITUDE)/(278.15)),5.255) * 100; //Pa
    //PRESSURE @ 900m : 905 hPa
    AIR_PRESSURE = 905 * pow((1 - (0.0065 * ALTITUDE)/(AIRTEMP)),5.255) * 100; //Pa
    //ALTERNATIVE?? INCLUDE PRESSURE IN WASIM???
    
    return AIR_PRESSURE;  
}
// LATENT HEAT TRANSFER IN CANOPY LAYER --> dependent on air temp within canopy --> Foerster!
//@LATHEAT
double snow_canopy::CALC_LATHEAT(double CANTEMP)
{
    double LATHEAT = 0;
    if (CANTEMP < 273.15){
        LATHEAT = 3.335e5; //Sublimation, Resublimation
    }
    else{
        LATHEAT = 2.257e6; //Vaporization, Condensation
    }    
    return LATHEAT;  
}
double snow_canopy::CALC_LE_CAN(double AIRTEMP, double GASCONSTDRYAIR, double CANTEMP, double HUMIDITY, double RSC, double RSI, double WATINTERCEP, double MAXWATINTERCEP, double WIND, double CANHEIGHT, double Z0, double ALTITUDE, double RHO_AIR)
{
    double SAT_PRES_AIR_CAN = 0;
    double LATHEAT = CALC_LATHEAT(CANTEMP);
    double C_e = 1;
    
    C_e = CALC_BULK_TRANS_COEFF(CANTEMP, AIRTEMP, WIND,  CANHEIGHT, Z0); //???
    
    //!!Einheit Vorfaktor!!! --> So Pa!! 6.112 -> hPa
    if (CANTEMP >= 273.15){
        SAT_PRES_AIR_CAN = 611.2 * exp((17.62 * (CANTEMP -273.15))/ (243.12 + (CANTEMP -273.15))); //Magnus
        //SAT_PRES_AIR_CAN = 6.107 * pow(10,((7.5*(CANTEMP -273.15))/(235 + (CANTEMP -273.15)))); //see Sonja
    }
    else {
        SAT_PRES_AIR_CAN  = 611.2 * exp((22.46 * (CANTEMP -273.15))/ (272.62 + (CANTEMP -273.15))); //Magnus
        //SAT_PRES_AIR_CAN = 6.107 * pow(10,((9.5*(CANTEMP -273.15))/(265.5 + (CANTEMP -273.15)))); //see Sonja
    }
    
    double RESIST_LE = CALC_RESIST_LE(RSC, RSI,WATINTERCEP, MAXWATINTERCEP ); //maybe include stability here
    
    double WIND_AC = CALC_WINDSPEED_ABOVE(WIND,CANHEIGHT,Z0);
    double AIR_PRESSURE = CALC_PRESSURE(ALTITUDE,AIRTEMP);
    
    double VAP_PRES_AIR_ATMO = VAP_PRESSURE(AIRTEMP, HUMIDITY);
    double KARMAN = 0.41;
    double FACTOR = 0.1;
    
    double LE_CAN = ((0.622*LATHEAT)/(GASCONSTDRYAIR*AIRTEMP*RESIST_LE))*(SAT_PRES_AIR_CAN - VAP_PRES_AIR_ATMO ); //old --> Gouttevin
    //double LE_CAN = RHO_AIR * LATHEAT * C_e * KARMAN * KARMAN * pow(log(15/Z0),-2) * WIND_AC * (1) ; // KUSTAS 1994 --> specific humidity??
    //double LE_CAN = ((0.622*LATHEAT)/(GASCONSTDRYAIR*AIRTEMP*RESIST_LE))* C_e *(SAT_PRES_AIR_CAN - VAP_PRES_AIR_ATMO); //stability!   -> läuft mit sensflux old
    //double LE_CAN =  ((RHO_AIR * 0.622*LATHEAT)/(AIR_PRESSURE))* C_e * WIND *(SAT_PRES_AIR_CAN - VAP_PRES_AIR_ATMO ); //dewalle!! --> for snowpack!  No resistances!!!
    
    return LE_CAN;  
}

//////////////////////////////////////////////////////////////////////////////////////
//ADVECTIVE FLUX RAIN SNOW
//@SENS_PREC_CAN
///////////////////////////////////////////////////////////////////////////////////////
double snow_canopy::CALC_PREC_HEAT_CAN(double AIRTEMP, double CANTEMP, double RAIN, double SNOW, double INTERVALL, double RSC, double RSI, double WATINTERCEP, double MAXWATINTERCEP )
{
    double RESIST_LE = CALC_RESIST_LE(RSC, RSI,WATINTERCEP, MAXWATINTERCEP ); //maybe include stability here
    //like Warscher + SNOW + RESISTANCES
    double SENS_PREC_CAN = (((RAIN/INTERVALL * 4180.0 * (AIRTEMP-CANTEMP)))+((SNOW/INTERVALL * 2100.0 * (AIRTEMP-CANTEMP))))/RESIST_LE;
    
    return SENS_PREC_CAN;  
}
/////////////////////////////////////////////////////////////////////////////////////////////
// ENERGY STORAGE IN BIOMASS OF CANOPY
// @BM_CAN
/////////////////////////////////////////////////////////////////////////////////////////////
double snow_canopy::CALC_BM_CAN(double CANTEMP, double CANTEMP_OLD, double LAI, double CANHEIGHT,double INTERVALL)
{   
//	return 0; //js: temporary disable this term
    double LEAF_THICKNESS = 1e-3; // [m]
    double BIOMASS_DENSITY = 900; // [kg/m�]
    double BIOMASS_HEAT_MASS = 2800; // [J/(kg*K)]
    double STAND_BASAL_AREA = 3 * 1e-3; // [m�/m�] --> implement in landuse_table!!! Calibration?? --> REFERENZ!!!
    
    //
//    double DELTA_t = 3600;
    
    double HM_LEAVES = LAI * LEAF_THICKNESS * BIOMASS_DENSITY * BIOMASS_HEAT_MASS;
    double HM_TRUNK = 0.5 * STAND_BASAL_AREA * CANHEIGHT * BIOMASS_DENSITY * BIOMASS_HEAT_MASS;
    
//	double BM_CAN = (HM_LEAVES + HM_TRUNK) * (CANTEMP - CANTEMP_OLD) / DELTA_t; //--> HOURLY!!!
	double BM_CAN = (HM_LEAVES + HM_TRUNK) * (CANTEMP - CANTEMP_OLD) / INTERVALL; //--> HOURLY!!!

    return BM_CAN;  
}
////////////////////////////////////////////////////////////////////////////////////////
//@CANENBAL
////////////////////////////////////////////////////////////////////////////////////////
double snow_canopy::CALC_CANENBAL(SnowCanopyParam* p)
{   
    //recalculate original COUNTERRADIATION with unaltered values of Airtemp
    p->IR = COUNTERRAD(p->ALTITUDE, p->AIRTEMP_OLD, p->HUM_OLD, p->PRECIP, p->CL);
    
    double CANENBAL = 0;
    double SW_NET_CAN       = CALC_SW_NET_CAN( p->LAI, p->GLOBRAD, p->ALBEDO_CAN, p->ALBEDO_SURF); //Albedo from SCN and wasim works!!! --> Albedo = gridinfo[getindex(_Albedo_)].pointer; in interzep.cc
	//js using dSurftemp as airtemp - 0.5 when surftemp is -9999 (which, by the way, should never happen since surftemp will always be available for SSTO>0 and CALC_ENBAL will be called for SSTO>0 only
	double dSurftemp = (p->SURFTEMP > -9999) ? p->SURFTEMP : p->AIRTEMP - 0.5;
	double LW_NET_CAN = CALC_LW_NET_CAN(p->LAI, dSurftemp, p->CANTEMP, p->IR);
//	double LW_NET_CAN = CALC_LW_NET_CAN(p->LAI, p->SURFTEMP, p->CANTEMP, p->IR);
	double LE_CAN = CALC_LE_CAN(p->AIRTEMP_OLD, GASCONSTDRYAIR, p->CANTEMP, p->HUM_OLD, p->RSC, p-> RSI, p->SI, p->SI_MAX, p->WIND_OLD, p->CANHEIGHT, p->Z0, p->ALTITUDE, p->RHO_AIR);
    double SENS_AIR_CAN     = CALC_SENS_AIR_CAN(p->CANHEIGHT, p->RHO_AIR, p->HEATCAP_AIR, p->CANTEMP, p->AIRTEMP_OLD, p->WIND_OLD, p-> Z0);
    double PREC_HEAT_CAN    = CALC_PREC_HEAT_CAN(p->AIRTEMP_OLD,  p->CANTEMP, p->RAINPRECIP, p->SNOWPRECIP, p->INTERVALL, p->RSC, p-> RSI, p->SI, p->SI_MAX);
    double BM_CAN           = CALC_BM_CAN(p->CANTEMP,p->CANTEMP_OLD,p->LAI, p->CANHEIGHT, p->INTERVALL); //working!!! --> cantemp < airtemp!!!

    // Calculation of energy balance
	CANENBAL = SW_NET_CAN + LW_NET_CAN - SENS_AIR_CAN - BM_CAN; //+ PREC_HEAT_CAN - LE_CAN - BM_CAN

//	PREC_HEAT_CAN = 0;
//	BM_CAN = 0;
//	CANENBAL = SW_NET_CAN + LW_NET_CAN + LE_CAN + SENS_AIR_CAN;

    return CANENBAL;
}

////////////////////////////////////////////////////////////////////////////////////////
//@SOLVE_CEB
////////////////////////////////////////////////////////////////////////////////////////

bool snow_canopy::SOLVE_CEB(SnowCanopyParam* p)
{   
    //loop parameters
    int itval = 0;
    int maxiter = 300;
    bool on = true;
    
    ////////////////////////////////////////
    //Bisection method
    ///////////////////////////////////////
    double epsilon = 1e-8;
    //1. Initial guess and calculation of CEB
	//double a = 200; //CANTEMP for positive CANENBAL
	//js
	double CEB_a=0;
	double a = 240; //CANTEMP for positive CANENBAL

	// for debugging
//	for (double te = -25.0; te <= 25.0; te += 1)
//	{
//		double a = p->AIRTEMP + te; //CANTEMP for positive CANENBAL
//		p->CANTEMP = a;
//		CEB_a = CALC_CANENBAL(p);
//	}
//	a = 200; 

	p->CANTEMP = a;
	CEB_a = CALC_CANENBAL(p);

	double b = 350; //CANTEMP for negative CANENBAL
	p->CANTEMP = b;
    double CEB_b = CALC_CANENBAL(p);
     
    double CEB_c = 666;
    //@todo CANHEIGHT IS CONSTANT!! --> < 3
    if(p->CANHEIGHT > 3. && on == true){
        
        if (CEB_a * CEB_b >= 0) 
        { 
			p->CANTEMP = a = p->AIRTEMP;
			double CEB_c = CALC_CANENBAL(p);
		    on = false;
        }
        
        double c = a;
        
        //ITERATIVE CALCULATION OF new CANTEMP
        while ((b - a) > epsilon && itval < maxiter && on == true){
            
            //middle point
            c = (a+b)/2;
            
            //calc CEB(c)
            p->CANTEMP = c;
            CEB_c = CALC_CANENBAL(p);
            
            // Check if middle point is solution
            if (fabs(CEB_c) < 1.0) 
                break; 
  
            // Decide the side to repeat the steps 
            else if (CEB_c * CEB_a < 0) 
                b = c; 
            else
                a = c;
            
            itval++;
        }
        
        //store Results    
        
        p->CANTEMP = c; //T_temp
        p->CANENBAL = CEB_c; 
        p->DUMMY = (p->CANTEMP - p->AIRTEMP_OLD)  ;//(p->CANTEMP - p->AIRTEMP_OLD) 
        
        /*
        //SET MAXDIFF CANTEMP <-> AIRTEMP
        double maxdiff = 30;
        if((p->CANTEMP - p->AIRTEMP_ABOVE) > maxdiff){
            p->CANTEMP = p->AIRTEMP_ABOVE + maxdiff;
            p->DUMMY = maxdiff;
            itval = maxiter;
        }
        if((p->CANTEMP - p->AIRTEMP_ABOVE) < - maxdiff){
            p->CANTEMP = p->AIRTEMP_ABOVE - maxdiff;
            p->DUMMY = - maxdiff;
            itval = maxiter;
            }
         */
    }
    //USE new CANTEMP for alteration of incoming IR --> 2 parts: canopy + atmosphere!!!
    
    double CANOPYFRACTION = 0.35 + 0.29 * log(p->LAI); //wie FOERSTER; VCF aus LU table?? 0.55 + 0.29 * log(p->LAI)
    if(CANOPYFRACTION < 0.0) CANOPYFRACTION = 0.0;
    if(CANOPYFRACTION > 1.0) CANOPYFRACTION = 1.0;
    
    
    double IR_ATMO = (1- CANOPYFRACTION) * p->IR;
    double emissivity = 0.98; //-> REFERENCE?
    double IR_CANOPY = CANOPYFRACTION * (emissivity * STEFBOLTZ * p->CANTEMP * p->CANTEMP * p->CANTEMP * p->CANTEMP); 
    p->IR = IR_ATMO + IR_CANOPY;
    
   	/*
    //alternative from Gouttevin  --> BETTER!!??
    double emissivity = 0.98;
    double CANOPTRANSM = exp(-EXTINCTCOEFF * p->LAI);
    double ABSORPFACT = (1 - CANOPTRANSM);
    p->IR = (1 - ABSORPFACT) * p->IR + ABSORPFACT * emissivity  * STEFBOLTZ * pow(p->CANTEMP,4); // - STEFBOLTZ * pow(p->SURFTEMP,4); 
    */   
    return true;
}