/* snow_canopy.h
 * 
 * This file includes definitions of constants and functions. It is mostly based on
 * the file canopy.fi which is part of the ESCIMO model written by U. Strasser
 * (Srasser et al., 2002, 2004, 2008, 2011)
 * 
 * The conversion to C++ was done for the implementation of the model into
 * the PANTA RHEI model in 2012 (Förster et al., 2014 HESS)
 * 
 */

#ifndef SNOW_CANOPY
#define SNOW_CANOPY

#include <vector>

namespace snow_canopy {

	// aus Escimo.fi
	static const double STEFBOLTZ = 5.67E-8; //!Stefan-Boltzmann-Konstante (W/(m**2*K**4))
	static const double CONSUBSNOW = 2835500.0;     //!Kondens./Sublim.-Waerme von Schnee (J/kg)

	// aus Canopy.fi
	static const double DIMLESSCOEFF = 0.01;    // ! kc   =  dimensionless coefficient = 0.010 related to the shape of the intercepted snow deposits
	static const double EXTINCTCOEFF = 0.41;    // ! k    =  vegetation-dependent extinction coefficient 0.41 --> orig:0.71
	static const double ICEDENSITY = 917.0;     // ! rhoi =  ice density (kg/m≥)
	static const double ICERADIUS = 0.0005;  // ! r    =  radius (m) of a sperical ice particle (500 µm)
	static const double MOLWEIGHTWAT = 0.01801; // ! M    =  kg/mole, molecular weight of water = 18.01 (g/mole)
	static const double UNIVGASCONST = 8.313;   // ! R    =  J/(mole K) universal gas constant
	static const double THERMCONDATMO = 0.024;   // ! lt   =  thermal conductivity of the atmosphere = 0.024 (J/(m s K))
	static const double GASCONSTDRYAIR = 287.0;     // ! Rd   =  gas constant of dry air = 287 (J/(deg kg))
	static const double KINVISCOSAIR = 0.000013;// ! v    =  kinematic viscosity of air = 1.3*10-5 (m≤/s)
	static const double SCALEFAC = 0.9;     // ! beta =  scaling factor = 0.9
	static const double GASCONSTAIR = 287.058;     // ! J/(kg*K)
	static const double KARMAN = 0.40;     // []

	// This class organizes data exchange between the interception model of WaSiM and the new snow canopy model
	class SnowCanopyParam
	{
	public:
		// METEO-INPUT
		double AIRTEMP;
		double RHO_AIR;
		double PRES_AIR;
		double HEATCAP_AIR;
		double MEANTEMP;
		double DELTATEMP;
		double HUM;
		double PRECIP;
		double SNOWPRECIP;
		double RAINPRECIP;
		double SHORTWAVE;
		double WIND;
		double IR;
		double CL;

		// PARAMETER
		double INTERVALL;
		double ALTITUDE;
		double ALBMIN;
		double ALBADD;
		double ABKLINGPOS;
		double ABKLINGNEG;
		double LAI;
		double SCALELAI;
		double CANHEIGHT;
		// $kf 2012-10-14: Parameter zum Testen
		double tIndexTemp;
		double tIndexRad;

		// State variables
		double ALBEDO;
		double CANOPINTLOAD;	// Canopy-SWE TFL
		double CANOPINTLOADOLD;
		double MELTUNLOAD_L;
		double MELTUNLOAD_S;
		double SUBLIMSNOW;	// (Re-)Sublimation Canopy
		double SNOWAGE;

		//Extension CANENBAL
		//albedo stuff
		double ALBEDO_CAN;
		double ALBEDO_CAN_DRY;
		double ALBEDO_CAN_WET;
		double FRACTION_WET;
		double ALBEDO_SURF;
		//Radiation
		double SURFTEMP;
		double CANTEMP;
		double CANTEMP_OLD;
		double LW_NET_CAN;
		double GLOBRAD;
		double SW_NET_CAN;
		//Turbulent fluxes
		double AIRTEMP_OLD;
		double HUM_OLD;
		double WIND_OLD;
		double RSC;
		double RSI;
		double Z0;
		double SI;
		double SI_MAX;
		double SAT_PRES_AIR_CAN;
		double VAP_PRES_AIR_ATMO;
		double H_CAN;
		double LE_CAN;
		//heat mass
		double BM_CAN; //
		//THE CANENBAL
		double CANENBAL;
		double DUMMY;


		SnowCanopyParam()
		{
			// Meteorological input
			this->AIRTEMP = 273.16; // KELVIN
			this->MEANTEMP = 273.16; // KELVIN
			this->RHO_AIR = 1.257; // kg/m^3
			this->PRES_AIR = 101300; // Pa
			this->HEATCAP_AIR = 1005; // J/(kg*K)
			this->DELTATEMP = 0.;
			this->HUM = 90.; // PERCENT
			this->PRECIP = 0.;
			this->SNOWPRECIP = 0.;
			this->RAINPRECIP = 0.;
			this->SHORTWAVE = 0.;
			this->WIND = 0.;
			this->IR = 316.; // downward longwave radiative flux (default value 273.15K)
			this->CL = 0.0;
			// Parameters
			this->INTERVALL = 3600.;
			this->ALTITUDE = 0.;
			this->ALBMIN = 0.1;
			this->ALBADD = 0.7;
			this->ABKLINGPOS = 0.3;
			this->ABKLINGNEG = 0.1;
			this->LAI = 4.;
			this->SCALELAI = 1.0;
			this->CANHEIGHT = 10.; // ABOVE CANOPY = 12 m
			this->tIndexTemp = 0.025;//0.05;
			this->tIndexRad = 0.005;//0.0094;
			// States variables
			this->ALBEDO = 0.2;
			this->CANOPINTLOADOLD = 0.;
			this->SNOWAGE = 0.;
			// Additiona results (no state vars)
			this->MELTUNLOAD_L = 0.;
			this->MELTUNLOAD_S = 0.;
			this->CANOPINTLOAD = 0.;
			this->SUBLIMSNOW = 0.;
			//PARAMETERS AND VARIABLES FOR CANENBAL
			this->SURFTEMP = 273.16; //Get from snowmod!!
			this->AIRTEMP_OLD = 272.16; //Get from snowmod!!
			this->ALBEDO_SURF = 0.75;
			this->ALBEDO_CAN = 0.1; //
			//this->FRACTION_WET   = 0.3;
			this->GLOBRAD = 0.;
			this->Z0 = 2.;
			//this->RSC             = 60;
			//this->RESIST_AERO_HEAT    = 50.;
			//this->RESIST_AERO_VAPOR  = 50.;
			this->CANTEMP = 273.16;
			this->CANTEMP_OLD = 273.15;

			//initialize SUMMANDS CANENBAL --> values written if SSTO = 0!!
			this->LW_NET_CAN = 0.;
			this->SW_NET_CAN = 0.;
			this->H_CAN = 0.;
			this->LE_CAN = 0.;
			this->BM_CAN = 0.;
			this->CANENBAL = -2000.; //> limit!
			this->DUMMY = 0.;
		}
	};



	// Functions
	bool CANOPY_INTERCEPTION(SnowCanopyParam* p);
	bool SIM_ALBEDO(double ALBMIN, double ALBADD, double ABKLINGPOS, double ABKLINGNEG, double SNOWPRECIP, double &ALBEDO, double &SNOWAGE, double AIRTEMP, double INTERVALL);
	double T_INDEX(double AIRTEMP, double ALBEDO, double DIR, double DIF, int INTERVALL, double TI_temp, double TI_rad);
	double min2(double a, double b);
	double max2(double a, double b);
	double VAP_PRESSURE(double Temperature, double Humidity);
	double WETBULB_TEMPERATURE(double AIRTEMP, double HUM, double ALTITUDE);
	bool PhasePartitioning(double AIRTEMP, double HUM, double PRECIP, double ALTITUDE, double &SNOWPRECIP, double &RAINPRECIP);
	void PhasePartitioningWaSiM(double temperatur, double T0R, double transzone, double precip, double &SNOWPRECIP, double &RAINPRECIP);
	double COUNTERRAD(double STATIONELEVATION, double AIRTEMP, double HUM, double PRECIP, double NB);
	//CEB
	double CALC_CANENBAL(SnowCanopyParam* p);
	bool SOLVE_CEB(SnowCanopyParam* p);

	double CALC_ZERO_DISPLACEMENT(double CANHEIGHT);
	double CALC_WINDSPEED_ABOVE(double WIND,double CANHEIGHT, double Z0);
	double CALC_WET_CAN_FRAC(double WATINTERCEP, double MAXWATINTERCEP);
	double CALC_SW_NET_CAN(double LAI, double SW_IN, double ALBEDO_CAN, double ALBEDO_SURF);
	double CALC_LW_NET_CAN(double LAI, double SURFTEMP, double CANTEMP, double LW_IN);
	double CALC_SNOW_DENSITY(double AIRTEMP);
	double CALC_SENS_AIR_CAN(double CANHEIGHT, double RHO_AIR, double HEATCAP_AIR, double CANTEMP, double AIRTEMP, double WIND, double Z0);
	double CALC_RICHARDSON(double CANTEMP, double AIRTEMP, double WIND, double CANHEIGHT, double Z0);
	double CALC_HEAT_VAPOR_ROUGHNESS(double CANHEIGHT);
	double CALC_MOMENTUM_ROUGHNESS(double CANHEIGHT);
	double CALC_RESIST_LE(double RSC, double RSI, double WATINTERCEP, double MAXWATINTERCEP);
	double CALC_LATHEAT(double CANTEMP);
	double CALC_LE_CAN(double AIRTEMP, double GASCONSTDRYAIR, double CANTEMP, double HUMIDITY, double RSC, double RSI, double WATINTERCEP, double MAXWATINTERCEP, double WIND, double CANHEIGHT, double Z0, double ALTITUDE, double RHO_AIR);
	double CALC_PREC_HEAT_CAN(double AIRTEMP, double CANTEMP, double RAIN, double SNOW, double INTERVALL, double RSC, double RSI, double WATINTERCEP, double MAXWATINTERCEP );
	double CALC_RESIST_AERO(double CANHEIGHT, double WIND, double Z0);
	double CALC_BULK_TRANS_COEFF(double CANTEMP, double AIRTEMP, double WINDSPEED, double CANHEIGHT,double Z0);
	double CALC_C_HN(double CANHEIGHT, double Z0);
	double CALC_BM_CAN(double CANTEMP, double CANTEMP_OLD, double LAI, double CANHEIGHT, double INTERVALL);
	double CALC_PRESSURE(double ALTITUDE, double AIRTEMP);
	double CALC_FRICTION_VELOCITY(double WIND);
	double CALC_ROUGHNESS_LENGTH(double WIND);

}

#endif
