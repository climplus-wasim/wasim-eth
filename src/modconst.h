/*   $Revision: 54 $ 
	$Date: 19.02.11 20:35 $  */
//   (c) 1995-2007 by Joerg Schulla
//   more version history: see source safe comments

#ifndef MODCONSTH
#define MODCONSTH

#undef  _ENGLISH_
#undef  _16_bit_
#define _ENGLISH_
//#define _16_bit_

#if defined _16_bit_
#define buffersize      (long)65535
#else
#define buffersize      (long)2097152 /* maximum size of the internal output buffer 2097152 = 2^21 = 2 MByte */
//#define buffersize      (long)1048576 /* maximum size of the internal output buffer 1048576 = 2^20 = 1 MByte */
#endif

#define _stringsize_      2400        /* alle Dateinamen und andere reservierte Strings haben diesen Speicherplatz zur Verfuegung */
#define _meteosatzlaenge_ 8000        /* maximale Laenge des Puffers, in den ein Satz aus meteorologischen Dateien im ganzen eingelesen werden kann. Sind die Saetze laenger, wird nachgelesen */
#define _Dicke_Schicht_1_old_ 0.05         /* Dicke der ersten Bodenschicht in m */
#define gridmaxanzahl     1000         /* maximale Anzahl von im SPeicher zu haltenden Grids */
#define maxzonenzahl      16000        /* soviele Zonen koennen im Zonengrid codiert sein */
#define maxzonenwert      (long)32767  /* maximale Groesse der Codes fuer die Zonen */
#define _max_regionen_zahl_ 30        /* anzahl an Regionen, die fuer die regression-interpolation genutzt werden koennen */

#define true              1
#define false             0
#define maxsymbolanzahl   2000         /* fuer Steuerdatei */

#ifdef WIN32
#define dirzeichen       '\\'         // Windows style directory seperator
#else
#ifdef WIN64
#define dirzeichen       '\\'         
#else
#define dirzeichen       '/'  // UNIX/LINUX style directory seperator
#endif
#endif
#define maxteilgebiete    10000        /* fuer Abflussbildungsmodell -> Bodenmodellklasse
									pro Teilgebiet werden im Bodenmodell 168 Byte belegt (zonale Statistik u.a.)
									bei 500 Teilgebieten also ca. 82 Kbyte */
#define Pi              (double)3.14159265359

#define _maxResults_            200  /* maximum number of entries in the lastresult-array */

#define NO_VALID_INSTANCE_INDEX 102
#define UNKNOWN_EXCEPTION       103

#define OTHERDATACOUNT 16 /* how many data columns will be writen to the control points "otherdata" file */

/* nun einige Konstanten, die fuer eine rationelle Speicherreservierung noetig sind,
besonders, wenn mit GNU-Compilern gearbeitet wird, die keine eigene Speicherverwaltung ala
Borland haben. Damit bei einer Matrix-reservierung nicht immer maximal eine halbe 2er-Potenz
Speicher verschenkt wird, wird der Speicher in optimalen Bloecken angefordert, die moeglichst
wenig Rest lassen, aber auch, moeglichst wenig Verwaltungsaufwand erfordern sollen, weshalb
ein maximaler overhead angegeben werden kann, der die Blance zwischen Speicherverschwendung und
minimalem Verwaltungsaufwand angibt. Der Wert sollte <= 0.05 sein.
Je nach Betriebssystem kann der maximal linear zu adressierende Speicher in _maxlinmem_ angegeben werden.
16-Bit hat 2^16-24 = 65520, 32-bit hat 2^32 Byte direkt zu adressieren.  */

#if defined _16_bit_
#define _maxlinmem_      (unsigned long)65520
#else
#define _maxlinmem_      (unsigned long)2147480000
#endif

#define _1_ln2_          (double)1.442695041111   /* 1/ln2 */
#define _maxmemoverhead_ (double)0.1    /* max. ca. 5% Speicherverschwendung sollen im Interesse kleinen Verwaltungsaufwandes moeglich sein */

/* hier werden die im Modell genutzten und z.T. auch in der Steuerdatei im Klartext angegebenen
Bezeichnungen einzelner Groessen, Teilmodelle und anderer Bezeichnungen als Makros definiert,
so kann spaeter das Umschreiben der Meldungen und der Steuerdatei in eine andere Sprache
schneller erfolgen. Ausserdem hat das den Vorteil, dass sich Aenderungen in der Benennung eines
Abschnittes schnell im Modell realisieren lassen, da nicht mehr jedes einzelne AUftreten
des Textes gesucht werden muss, sondern die Aenderung in dieser Datei ueberall wirksam ist. */


/* Abschnitte und Bezeichnungen aus der Steuerdatei */
#define _Programmausgabe_             "[output_interval]"
#define _output_list_                 "[output_list]"
#define _Koordinaten_                 "[coordinates]"
#define _region_transition_distance_  "[region_transition_distance]"
#define _Hoehenmodell_                "[elevation_model]"
#define _Zonengrid_                   "[zone_grid]"
#define _Standardgrids_               "[standard_grids]"
#define _variable_Grids_              "[variable_grids]"
#define _Modellzeit_                  "[model_time]"
#define _anzahl_                      "[meteo_data_count]"
#define _meteonamen_                  "[meteo_names]"
#define _Niederschlagskorrektur_      "precipitation_correction"
#define _Expositionskorrektur_        "radiation_correction"
#define _Verdunstungsmodell_          "evapotranspiration"
#define _Schneemodell_                "snow_model"
#define _Interzeptionsmodell_         "interception_model"
#define _SiltingUpModell_			  "SiltingUpModel"
#define _SurfaceRoutingModell_		  "SurfaceRoutingModel"
#define _Infiltrationsmodell_         "infiltration_model"
#define _Bodenmodell_                 "soil_model"
#define _Unsatzonmodell_              "unsatzon_model"
#define _gwflow_                      "groundwater_flow"
#define _irrigation_amounts_          "irrigation_amounts"
#define _irrigation_parameters_       "irrigation"
#define _Abflussmodell_               "routing_model"
#define _Defaults_                    "defaults"
#define _MultiLayerLUSETable_         "[multilayer_landuse]"
#define _Landnutzungstabelle_         "[landuse_table]"
#define _Bodenartentabelle_           "[soil_table]"
#define _irrigation_table_            "[irrigation_table]"
#define _lake_model_				  "[lake_model]"
#define _RegionalSuperposition_       "RegionalSuperposition"
#define _heat_transfer_				  "[heat_transfer]"
#define _special_output_			  "special_output"
#define _data_assimilation_			  "data_assimilation"
#define _unsatzon_tuning_			  "unsatzon_tuning"
#define _groundwater_tuning_		  "groundwater_tuning"
/* Meteogroessen */

#define _Globalstrahlung_             "global_radiation"
#define _Nettostrahlung_              "net_radiation"
#define _Sonnenscheindauer_           "sunshine_duration"
#define _SunshineDurationAbsolute_    "sunshine_duration_absolute"
#define _Temperatur_                  "temperature"
#define _Temperatur_14_               "temperature_14"
#define _Windgeschwindigkeit_         "wind_speed"
#define _Dampfdruck_                  "vapor_pressure"
#define _Dampfdruck_14_               "vapor_pressure_14"
#define _Luftfeuchte_                 "air_humidity"
#define _Luftfeuchte_14_              "air_humidity_14"
#define _Luftdruck_                   "air_pressure"
#define _Niederschlag_                "precipitation"

/* Standardgrids und veraenderliche Standardgrids */

#define _regression_regions_1_        "regression_regions"
#define _regression_regions_2_        "regression_regions2"
#define _regression_regions_3_        "regression_regions3"
#define _regression_weights_1_        "regression_weights1"
#define _regression_weights_2_        "regression_weights2"
#define _regression_weights_3_        "regression_weights3"
#define _RelCellArea_                 "RelCellArea"
#define _CellSizeX_                   "CellSizeX"
#define _CellSizeY_                   "CellSizeY"
#define _CoordY_                      "CoordY"
#define _CoordX_                      "CoordX"
#define _Bodenarten_                  "soil_types"
#define _Landnutzung_                 "landuse"
#define _ice_firn_                    "ice_firn"
#define _Nutzbare_Feldkapazitaet_     "available_soil_moisture"
#define _Gefaelle_                    "slope_angle"
#define _Exposition_                  "slope_aspect"
#define _Topographischer_Faktor_      "topographic_faktor"
#define _Zonengrid_Bodenmodell_       "zonegrid_soilmodel"
#define _Hydraulische_Leitfaehigkeit_ "hydraulic_conductivity"
#define _SkyViewFactor_               "sky_view_factor"
#define _Fliesszeiten_                "flow_times"
#define _Albedo_                      "albedo"
#define _Bodenspeicher_               "soil_storage"
#define _max_ponding_storage_         "max_ponding_storage"
#define _clay_depth_                  "clay_depth"
#define _irrigation_codes_            "irrigation_codes"
#define _ETR2ETP_					  "ETR2ETP"
#define _river_depth_                 "river_depth"
#define _river_width_                 "river_width"
#define _river_links_                 "river_links"
#define _kolmation_                   "kolmation"

/* Speichergrids fuer Zustaende/Fluesse */

#define _Auffuellbare_Porositaet_     "fillable_porosity"
#define _Saettigungsdefizit_          "saturation_deficite"
#define _Bestandesniederschlag_       "throughfall"
#define _Schneeausfluss_              "snowcover_outflow"
#define _Snow_out_                    "melt_from_snowcover"
#define _Schneetage_                  "days_snow"
#define _Schneealter_                 "snow_age"
#define _Schneeanteil_                "snow_fraction"
#define _Regenanteil_                 "rain_fraction"
#define _melt_from_firn_              "melt_from_firn"
#define _melt_from_ice_               "melt_from_ice"

/* #define _Wasseraequivalent_           "Wasseraequivalent" */

#define _Unsaturated_Zone_Storage_    "unsaturated_zone_storage"
#define _Interflow_Storage_           "interflow_storage"
#define _Saturated_Soil_Partition_    "saturated_soil_partition"
#define __sbalt_grid__                "$$sbold_grid$$"
#define _Bodenspeicher_Alt_           "soil_storage_old"
#define _Basisabfluss_                "baseflow"
#define _Direktabfluss_               "directflow"
#define _Interflow_                   "interflow"
#define _Peakflow_                    "peakflow"
#define _Nacht_Temperatur_            "temperature_night"
#define _Korrigierte_Temperatur_      "corr_temperature"
#define _Potentielle_Verdunstung_     "pot_transpiration"
#define _Reale_Verdunstung_           "real_transpiration"
#define _Potential_Evaporation_       "potential_evaporation"
#define _Reale_Evaporation_           "real_evaporation"
#define _Snow_Evaporation_            "snow_evaporation"
#define _Potential_InterceptionEvap_  "potential_InterceptionEvap"
#define _Reale_InterceptionEvap_      "real_InterceptionEvap"
#define _Nadelwald_                   "pineforest"
#define _Schatten_                    "Shadow"
#define _IEFF_MAX_                    "Ieff_max"
#define __Nacht_Temperatur__          "$$temperature_night$$"
#define _Infiltrationsexzess_         "infiltration_excess"
#define _Total_Infiltration_		  "total_infiltration"
#define _Infiltrationssaettigung_     "infiltration_saturation"
#define _ungueltiger_Name_            "ungueltiger_Name"
#define _Interzeptionsspeicher_       "interception_storage"
#define _Interzeptionsverdunstung_    "interception_evaporation"
#define _snowstorage_solid_			  "snowstorage_solid"
#define _snowstorage_liquid_		  "snowstorage_liquid"
#define _Wurzeltiefe_                 "root_depth"
#define _Vegetationsbedeckungsgrad_   "vegetation_coverage_degree"
#define _Blattflaechenindex_          "leaf_area_index"
#define _inaktive_Bodenfeuchte_       "inaktive_soil_moisture"
#define _inaktive_Feldkapazitaet_     "inaktive_field_capacity"

// $kf 2016-06-20: Add grid names for snow canopy model
#define _Snow_Canopy_Albedo_          "Snow_canopy_albedo"
#define _Snow_Canopy_Load_            "Snow_canopy_load"
#define _Snow_Canopy_Age_             "Snow_canopy_age"
#define _Snow_Canopy_Sublim_			"Snow_canopy_sublim"
#define _Snow_Canopy_Unload_solid_    "Snow_canopy_unload_solid"
#define _Snow_Canopy_Unload_liquid_   "Snow_canopy_unload_liquid"
#define _Snow_Canopy_Temperature_     "Snow_canopy_temperature"
#define _Snow_Canopy_RH_              "Snow_canopy_rh"
#define _Snow_Canopy_RAD_             "Snow_canopy_rad"
#define _Snow_Canopy_Wind_            "Snow_canopy_wind"
#define _Snow_Canopy_IR_	      "Snow_canopy_IR_radiation"
//CEB
#define _Canopy_Energy_Balance_     "Canopy_Energy_balance"
#define _Canopy_Temperature_              "CEB_Canopy_Temperature"
#define _Canopy_Temperature_Old_          "CEB_Canopy_Temperature_Old"
#define _Canopy_Dummy_              "CEB_DUMMY"



/* einige Definitionen fuer die Abflussoptimierung mit flow_opt */

#define _Abflussoptimierung_          "discharge_optimizing"
#define _Mh_                          "Mh"
#define _Bh_                          "Bh"
#define _kh_                          "kh"
#define _Th_                          "Th"
#define _Mv_                          "Mv"
#define _Bv_                          "Bv"
#define _kv_                          "kv"

/* fuer Routingclass */
#define _TG_                          "TG"
#define _OL_                          "OL"
#define _SUMTRIB_                     "SUMTRIB"
#define _ZL_                          "ZL"
#define _AL_                          "AL"
#define _SP_                          "SP"
#define _constant_                    "constant"
#define _Schwelle_                    "Threshold"
#define _extern_                      "extern"
#define _intern_                      "intern"
#define _extern_with_rule_			  "extern_with_rule"
#define _intern_with_rule_			  "intern_with_rule"
#define _intern_with_rule_from_reservoir_  "intern_with_rule_from_reservoir"
#define _intern_from_file_from_reservoir_  "intern_from_file_from_reservoir"
#define _from_file_					  "from_file"

#define _Theta_                         "theta"
#define _hydraulic_head_                "hydraulic_heads"
#define _geodetic_height_               "geodetic_heads"
#define _soil_bottom_fluxes_            "q_down"
#define _Concentrations_                "concentrations"
#define _Grundwasserstand_              "groundwater_distance"
#define _Theta_in_GWLEVEL_              "Theta_in_GWLEVEL"
#define _Grundwasserneubildung_         "groundwater_recharge"
#define _Grundwasserlevel_              "groundwater_level"
#define _Drainageabfluss_               "drainage_flow"
#define _Drainagedistance_              "drainage_distance"
#define _Drainagehoehen_                "drainage_depth"
#define _Drainagelevel_                 "drainage_level"
#define _Tangens_Beta_                  "tangens_beta"
#define _ponding_storage_               "ponding_storage"
#define _ponding_storage_top_           "ponding_storage_top"
#define _Max_Pumping_                   "Max_Pumping"
#define _SB_1_grid_                     "SB_1_grid"
#define _SB_1_REL_grid_					"SB_1_REL_grid"
#define _SB_2_grid_                     "SB_2_grid"
#define _ROOTgrid_                      "ROOTgrid"
#define _SBLOW_grid_                    "SBLOW_grid"
#define _SBHIGH_grid_                   "SBHIGH_grid"
#define _THETA_LIMIT_FOR_ETP_REDUCTION_ "theta_limit_for_etp_reduction"
#define _river_infiltration_            "river_infiltration"
#define _river_exfiltration_            "river_exfiltration"
#define _groundwater_balance_           "groundwater_balance"
#define _macropore_infiltration_        "macropore_infiltration"
#define _number_of_substeps_			"number_of_substeps"

#define _gw_k_x_                        "gw_k_x_"
#define _gw_k_y_                        "gw_k_y_"
#define _gw_k_x_1_ori_                  "gw_k_x_1_ori"
#define _gw_k_y_1_ori_                  "gw_k_y_1_ori"
#define _gw_head_old_                   "gw_head_old_"
#define _gw_head_                       "gw_head_"
#define _gw_storage_coeff_              "gw_storage_coeff_"
#define _gw_boundary_fix_h_             "gw_boundary_fix_h_" 
#define _gw_boundary_fix_q_             "gw_boundary_fix_q_"
#define _gw_transmissivity_x_           "gw_transmissivity_x_"
#define _gw_transmissivity_y_           "gw_transmissivity_y_"
#define _aquifer_thickness_             "aquifer_thickness_"
#define _gw_flow_x_                     "gw_flow_x_"
#define _gw_flow_y_                     "gw_flow_y_"
#define _gw_kolmation_                  "gw_kolmation_"
#define _gw_leakage_                    "gw_leakage_"
#define _aquifer_upper_limit_           "aquifer_upper_limit"
#define _gw_balance_                    "gw_balance"

#define _actual_parameters_           "actual parameters"
#define _Meteodatenzeit_              "(time of meteo data)"
#define _Schneeverdunstung_           "snow_evaporation"
#define _ETP_                         "pot. evaptranspir."
#define _ETR_                         "real evapotranspir."
#define _Scenarios_                   "using scenarios!"
#define _suchradius_                  "no station within the specified search radius!"

/* transport algorithms */

#define _tracer_ini_                  "tracer_ini"
#define _substance_transport_         "substance_transport"
#define _radioactive_                 "radioactive"
#define _non_radioactive_             "non_radioactive"
#define _half_time_                   "half_time"
#define _days_                        "days"
#define _mgl_                         "mg/l"
#define _kgkg_                        "kg/kg"
#define _gl_                          "g/l"
#define _evaporating_                 "evaporating"
#define _non_evaporating_             "non_evaporating"
#define _min_conc_                    "min_conc"
#define _max_conc_                    "max_conc"
#define _soilini_                     "soilini"
#define _gwini_                       "gwini"


// definitions for new WaSiM-ETH 7 and higher does not use german definitions any longer 

#define _ExternalCoupling_             "[ExternalCoupling]"
#define _GWTableExtern_                "GWTableExtern"
#define _RenameMode_Delete_            0
#define _RenameMode_Rename_            1
#define _BALANCE_                      "Balance"
#define _SurfaceCanopyResistance_      "SurfaceCanopyResistance"
#define _SurfaceInterzepResistance_    "SurfaceIntercepResistance"
#define _SurfaceEvaporationResistance_ "SurfaceEvaporationResistance"
#define _RoughnessLength_              "RoughnessLength"
#define _Capillar_Uprise_              "Capillar_Uprise"
#define _Percolation_                  "Percolation"
#define _Schneespeicher_Total_         "snow_storage_total"
#define _AspectCorrectedGlobalRadiation_ "AspectCorrectedGlobalRadiation"

// definitions for dynamic phenology
#define _SumOfForcingUnits_				"SumOfForcingUnits"
#define _Pheno_start_					"Pheno_start"
#define _SumOfChillingUnits_			"SumOfChillingUnits"
#define _FStar_ForcingThreshold_		"FStar_ForcingThreshold"
#define _MultipleStartTimes_	    	"MultipleStartTimes"
#define _THROOT_grid_					"MeanThetaInRootZone"

// definitions for silting up and surface routing model
#define _SILTING_UP_DISPOSITION_		"SiltingUpDisposition"
#define _RAINPAUSE_						"RainPause"
#define _KINETIC_ENERGY_				"KineticEnergy"
#define _TIMESINCELASTSOILTILLAGE_		"TimeSinceLastSoilTillage"
#define _DIRECTRUNOFFSILTINGUP_			"DirectRunoffSiltingUp"
#define _SURFSTORSILTINGUP_				"SurfStorSiltingUp"
#define _FLOWDIRECTION_					"FlowDirection"
#define _DIRECTRUNOFFSURFACEROUTING_    "DirectRunoffSurfaceRouting"
#define _SURFACEROUTING_VELOCITY_       "SurfaceRoutingVelocity"
#define _SURFACEROUTING_FLUXES_         "SurfaceRoutingFluxes"
#define _SURFACEROUTING_WEIGHTS_		"SurfaceRoutingWeights"
#define _SURFACEROUTING_HEAD_		    "SurfaceRoutingHead"
#define _SURFACEROUTING_HEAD_OLD_		"SurfaceRoutingHeadOld"
#define _SURFACEROUTING_ROUGHNESS_X_	"SurfaceRoutingRoughnessX"
#define _SURFACEROUTING_ROUGHNESS_Y_	"SurfaceRoutingRoughnessY"
#define _SURFACEROUTING_MANNING_COEFF_	"SurfaceRoutingManningCoeff"
#define _PLANCURVATURE_                 "Plancurvature"
#define _SURFACEROUTING_QIN_            "SurfaceRoutingQIn"
#define _SURFACEROUTING_QOUT_           "SurfaceRoutingQOut"
#define _CRITICAL_SHEAR_STRESS_         "CriticalShearStress"

// definitions for the lake model
#define _LAKE_CODES_					"lake_codes"
#define _mean_ponding_storage_          "mean_ponding_storage" // max_ponding_storage has been used already for normal pond calculation 
#define _VegetationStart_               "VegetationStart"
#define _VegetationStop_                "VegetationStop"
#define _VegetationDuration_            "VegetationDuration"
#define _TemperatureAVG_			    "TemperatureAVG"
#define _TemperatureMOVAVG7_		    "TemperatureMOVAVG7"
#define _StartDayOfLeaveFall_			"StartDayOfLeaveFall"

// MPI-definitions
#define _MPIRegions_					"MPIRegions"
#define _MPIHelper_						"MPIHelper"
#define _MPIHalo_						"MPIHalo"

// definitions for Permafrost model
#define _permafrost_					"permafrost"
#define _ThawCoeffPermaFrost_           "ThawCoeffPermaFrost"
#define _SnowFreeDaysGrid_              "SnowFreeDaysGrid"
#define _SnowCoverDaysGrid_             "SnowCoverDaysGrid"
#define _PermafrostThawDepth_           "PermafrostThawDepth"

// definitions for dynamic glacier model
#define _firn_WE_stack_					"firn_WE_stack"
#define _GlacierMassBalance_			"GlacierMassBalance"
#define _OldGlacierMassBalance_			"OldGlacierMassBalance"
#define _GlacierizedCells_				"GlacierizedCells"
#define _GlacierizedCellsOld_			"GlacierizedCellsOld"
#define _GlacierCodes_					"GlacierCodes"
#define _GlacierCodesOld_				"GlacierCodesOld"
#define _debris_on_glaciers_			"debris_on_glaciers"
#define _soil_surface_groundwater_substeps_ "[soil_surface_groundwater_substeps]"
#define _TotalSnowStorage_onGlacier_	"TotalSnowStorage_onGlacier"
#define _SolidSnowStorage_onGlacier_	"SolidSnowStorage_onGlacier"
#define _LiquidSnowStorage_onGlacier_	"LiquidSnowStorage_onGlacier"
#define _SnowMelt_onGlacier_			"SnowMelt_onGlacier"
#define _SnowEvaporation_onGlac_        "SnowEvaporation_onGlacier"

// definitions for soil stability index
#define _root_cohesion_             "root_cohesion"
#define _soil_cohesion_             "soil_cohesion"
#define _wet_soil_density_			"wet_soil_density"
#define _internal_friction_angle_   "internal_friction_angle"
#define _soil_depth_				"soil_depth"
#define _soil_stability_index_		"soil_stability_index"

// definition for temperature modelling
#define _T_Lower_Boundary_Condition_ "T_Lower_Boundary_Condition"
#define _HeatFlux_Lower_Boundary_Condition_ "HeatFlux_Lower_Boundary_Condition"
#define _Soil_Temperature_Stack_     "Soil_Temperature_Stack"
#define _Soil_SE_Stack_				 "Soil_SE_Stack"
#define _ThawDepthTMod_              "ThawDepthTMod"
#define _frozenSoilInGroundwater_  "frozenSoilInGroundwater"
// definitions for temperature dependent colmation modeling between groundwater layers (aquifers)
#define _colmation_frozen_			"gw_1_colmation_frozen"
#define _colmation_unfrozen_		"gw_1_colmation_unfrozen"

// some definitions for ArcVeg integration

#define _IceWedgeUpperDepth_		"IceWedgeUpperDepth"
#define _IceWedgeLowerDepth_		"IceWedgeLowerDepth"
#define _ExcessIceFactor_			"ExcessIceFactor"
#define _AVINHLP_					"AVINHLP"
#define _AVGDD_						"AVGDD"
#define _AVINMAX7DAVG_				"AVINMAX7DAVG"
#define _AVTDAVG_					"AVTDAVG"
#define _AVINDAVG_					"AVINDAVG"
#define _AVSNOWFREE_				"AVSNOWFREE"
#define _AVPGP_						"AVPGP"
#define _AVINSCAL_					"AVINSCAL"
#define _AVNMIN_                    "AVNMIN"
#define _AVNMINTOT_                 "AVNMINTOT"
#define _AVSON_                     "AVSON"
#define _AVBIOM_                    "AVBIOM"
#define _ArcVegIntegration_         "[ArcVegIntegration]"

// definitions for SpecialOutputClass

// flags used when checking Lusetable for flags regarding the external or internal calculation in int_tvar_params
#define ALBEDOFLAG 0x01 // = 1   or 00000001b
#define RSEFLAG    0x02 // = 2   or 00000010b
#define RSCFLAG    0x04 // = 4   or 00000100b
#define RSIFLAG    0x08 // = 8   or 00001000b
#define LAIFLAG    0x10 // = 16  or 00010000b
#define VCFFLAG    0x20 // = 32  or 00100000b
#define ROOTFLAG   0x40 // = 64  or 01000000b
#define Z0FLAG     0x80 // = 128 or 10000000b
// "or"-ed combinations of these flags should be used, the default is 255/11111111b

#define PTYPE 1
#define TTYPE 0

// additional definitions for snow model extensions (Warscher, 2011-013)
#define _snowtemp_                    "snowtemp"
#define _surftemp_                    "surftemp"
#define _SurfaceTemperature_          "surfacetemperature"
#define _snowtemp_onGlac_             "snowtemp_onGlacier"
#define _surftemp_onGlac_             "surftemp_onGlacier"

#define _radflux_                     "radflux"
#define _sensflux_                    "sensflux"
#define _latflux_                     "latflux"
#define _adflux_                      "adflux"
#define _soilflux_                    "soilflux"
#define _slidefraction1_              "slidefraction1"	
#define _slidefraction2_              "slidefraction2"
#define _slidefraction3_              "slidefraction3"
#define _slidefraction4_              "slidefraction4"
#define _deposition_max_              "depostion_max"
#define _elev_order_row_              "elev_order_row"
#define _elev_order_col_              "elev_order_col"
#define _input_mass_                  "input_mass"
#define _mobile_mass_                 "mobile_mass"
#define _deposition_                  "deposition"
#define _max_deposition_              "max_deposition"
#define _depoindex_                   "deposition_index"
#define _erosindex_                   "erosion_index"
#define _Cloudiness_                  "cloudiness"

// Spinup module
#define _SpinUp_					    "[SpinUp]"
#define _HydrologicResponseUnits_	    "HydrologicResponseUnits"
#define _MasterHydrologicResponseUnits_ "MasterHydrologicResponseUnits"
#define _FreezingDegreeDays_			"FreezingDegreeDays"
#define _ThawingDegreeDays_				"thawingDegreeDays"
#define _DailyTemperature_				"DailyTemperature"
#define _AnnualTemperature_				"AnnualTemperature"

#define _ScalingGrid_					"ScalingGrid"

#define _MPIFullGridOn_					true
#define _MPIFullGridOff_				false
#define _LimitAutoFillOn_				true
#define _LimitAutoFillOff_				false
#define _AllocateGridMemoryOn_			true
#define _AllocateGridMemoryOff_			false
#define _ShowMessagesOn_				true
#define _ShowMessagesOff_				false
#define _MPIRank0OnlyOn_				true
#define _MPIRank0OnlyOff_				false
#define _MPISyncHaloOn_					true
#define _MPISyncHaloOff_				false

#define _PrecOnLakes_					"PrecOnLakes"

// additional defines for layered snow model
#define _EnergyBalance_					"EnergyBalance"
#define _EnergyBalance_onGlac_			"EnergyBalance_onGlacier"
#define _LayeredSnowModel_				"[LayeredSnowModel]"
#define _LayeredSnowProperties_         "[LayeredSnowProperties]"
#define _Theta_Snow_                    "theta_Snow"
#define _hydraulic_head_Snow_           "hydraulic_heads_Snow"
#define _geodetic_height_Snow_          "geodetic_height_Snow"
#define _soil_bottom_fluxes_Snow_       "q_down_Snow"
#define _Soil_Temperature_Stack_Snow_   "Temperature_Stack_Snow"
#define _SnowWaterEquivalentSolid_		"SnowWaterEquivalentSolid"
#define _SnowWaterEquivalentLiquid_		"SnowWaterEquivalentLiquid"
#define _SnowDensity_					"SnowDensity"
#define _SnowLayerThickness_			"SnowLayerThickness"
#define _upperIndexSnow_				"upperIndexSnow"
#define _fromRefreezingTrue_			true
#define _fromRefreezingFalse_			false
#define _onGlacierTrue_					true
#define _onGlacierFalse_				false

#define _MAX_OPEN_NETCDF_FILES_			1000

// defeinitions for output of actual field capacity values
#define _usableFieldCapacity_mm_       "usableFieldCapacity_mm"  // mm water in root zone at actual soil moisture content; output valoue and also used for calculating percentage of usable water content 
#define _usableFieldCapacity_perc_	   "usableFieldCapacity_perc" // quotient of usableFieldCapacity_mm / maxFieldCapacity_mm * 100%
#define _maxFieldCapacity_mm_          "maxFieldCapacity_mm"  // mm water in root zone at field capacity (depending on parameter psi_FC per soil layer) used for calculating percentage of usable water content 
#define _minFieldCapacity_mm_	       "minFieldCapacity_mm"  // mm water in rrot zone at theta_res 

// definitions for computing soil types from grain size fraction grids for clay, sand and silt
#define _clayfraction_					"clay_fraction" 
#define _siltfraction_					"sand_fraction" 
#define _sandfraction_					"silt_fraction" 

#endif // enclosing ifdef for modconst.h

