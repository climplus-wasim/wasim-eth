#ifdef _NETCDF

#ifndef NETCDFHLPH
#define NETCDFHLPH

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <cmath>

#include <netcdf>
#include "stringp.h"
#include "modconst.h"
#include "wsstream.h"

class netcdfhlp {

public:
	
	// Given a file name for a grid file, figure out which variable
	// is to be written out. Do this by removing slashes, underscores 
	// and dots and the strings outside of these separators
	static void findVarNameFromGridName(char* varname, char* gridname);
	static void findVarNameFromGridNameExt(char* varname, char* gridname);
	
	// Given a file name for a grid file, figure out a file name for the
	// netcdf-file
	static void findFileNameFromGridName(char* filename, char* gridname, bool isSumFile);
	
};


// Manage knowledge about coordinates etc.
class netcdfGeometry {
	
public:
	netcdfGeometry();
	~netcdfGeometry();

	void setGeometry(double xll, double yll, double cellsize);

	double lat(double x, double y);
	double lon(double x, double y);

protected:
	double xll_, yll_, cellsize_;

private:
	int sgn(double val);

};


class netcdfWriter : public netcdfGeometry {
	
public:
	
#ifdef SG
	netcdfWriter();
#else
        netcdfWriter(wasim::mystream* WSstream);
#endif
        ~netcdfWriter();
	void setLatGrid(double* latgrid);
	void setLonGrid(double* longrid);
	void writeTimestep(bool doublemode, double* dataD, float* dataF,
		char* filename, char* varname, int nrows,
		int ncols, double nodata, double cellsize,
		int year = -1, int month = -1, int day =-1,int hour=-1,int minute=-1);
	
private:
	
#ifndef SG
	wasim::mystream* WSstream;
#endif
	bool has_latgrid_;
	bool has_longrid_;
	double *latgrid_;
	double *longrid_;
	
	// Maintain a list of files we have already written to so that we 
	// can append later time steps to the correct file
	std::vector<char*> openFiles_;
	bool fileIsOpen(char* filename);
	void setFileOpen(char* filename);

	// Maintain a list of variables we have already written to a given 
	// file. For every new variable in every new file, we push a new 
	// pair to the following vectors.
	std::vector<char*> listOfFilesForVars_;
	std::vector<char*> listOfVarsInFile_;
	bool varExistsInFile(char* filename, char* varname);
	void setVarInFile(char* filename, char* varname);

};


class netcdfReader : public netcdfGeometry {
	
public:
	
#ifdef SG
        netcdfReader();
#else
	netcdfReader(wasim::mystream* WSstream);
#endif
	~netcdfReader();
	
//	void readTimestep(bool doublemode, double *dataD, float *dataF,
//		char* filename, int nrows, int ncols, double nodata);
	void readTimestep(bool doublemode, double* &stackbufD, char* filename,char *pvarname,
		int &nc_rows, int &nc_cols, double &nc_xll, double &nc_yll, double &nc_cellsize,
		double &nc_nodata, int &nc_levels,int actLevel=0,
		int year = -1, int month = -1, int day = -1, int hour = -1, int minute = -1);

	std::string readTimeValue(char* filename,int actLevel=0);

private:
	netCDF::NcFile *dataFile;
	std::map<std::string, netCDF::NcFile*> mapDataFiles;
	bool isTimeVariable;
#ifndef SG
	wasim::mystream* WSstream;
#endif        

};

#endif

#endif
