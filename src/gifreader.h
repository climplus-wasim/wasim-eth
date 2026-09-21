// class for reading Swiss GIF Radar Images 

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <vector>

typedef struct {
	int prevCode;
	char code;
} codeType;

class GIFReader {

private:

	char *buffer;
	std::ifstream s;
	std::string infilename;
//	std::stringstream &output;

	std::vector<codeType> vCodes;

	char header[6];
	char logScreenDescr[7];
	char *colortable;
	char GraphContrlExt[8];
	char *comment;
	char imgDescr[10];
	char blockidentifyer;
	char extensionidentifyer;

	int nextCodeExtensionIndex;
	int actCodeLength;
	int mincodesize;
	int endOffset;

	unsigned long decoderDWord;
	unsigned long tempDWord;
	unsigned long CodeMask;
	unsigned long actCode;
	unsigned long lastCode;
	unsigned long clearCode;
	unsigned long EOIcode;

	unsigned int gifcols;
	unsigned int gifrows;
	unsigned int infofield;
	unsigned int colortablesize;

	double **grid;
	int rowcnt;
	int colcnt;

void initializeCodeVector(int numColors);
void writeOutputToStream(int index);
char getFirstIndexOfCode(int index);
void addNewCode(int prevCode,char code);

public:
GIFReader(std::string infilename);
//GIFReader(std::string infilename,std::stringstream &output);
~GIFReader();
void ReadGIF(double **outputgrid);
void CutModelGrid(double **GIFOrig,int mincolsgo,int maxcolsgo,int minrowsgo,int maxrowsgo,double xllgo,double yllgo,double csgo,
				  double **ModelGrid,int mincolsmod,int maxcolsmod,int minrowsmod,int maxrowsmod,double xllmod,double yllmod,double csmod,double **maskgrid=NULL);
void GetDimensions(int &width,int &height);

void MapValues(double **ModelGrid, int mincolsmod, int maxcolsmod, int minrowsmod, int maxrowsmod, double f1 = 8.916787e-5, double x1 = 0.0767528404, int cut = 127, double f2 = 6.00e-006, double x2 = 0.098);
//void MapValues(double **ModelGrid, int mincolsmod, int maxcolsmod, int minrowsmod, int maxrowsmod, double f1, double x1, int cut, double f2, double x2);
//void MapValues(double **ModelGrid, int mincolsmod, int maxcolsmod, int minrowsmod, int maxrowsmod);

};
