// class for reading Swiss GIF Radar Images 

#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <memory.h>
#include "gifreader.h"


void GIFReader::initializeCodeVector(int numColors)
{
	codeType ctCode;
	vCodes.clear();
	for (int i=0;i<=numColors+1;i++) // including clear code and EOI code it's two more than the colortable length
	{
		ctCode.code=i;
		ctCode.prevCode=-1;
		vCodes.push_back(ctCode);
	}

	CodeMask=0;
	for (int i=0;i<mincodesize;i++)
	{
		CodeMask = CodeMask | (0x0000001<<i);
	}
	nextCodeExtensionIndex = (1<<mincodesize) -1;
	actCodeLength=mincodesize;
}

// recurively getting the output sequence for the given code
void GIFReader::writeOutputToStream(int index)
{
	if (vCodes[index].prevCode>=0)
	{
		writeOutputToStream(vCodes[index].prevCode);
	}
	unsigned short int value = (unsigned short int)(vCodes[index].code & 0x0ff);
	if (value==255)
	{
		if (grid)
		{
			grid[rowcnt][colcnt]=-9999;
			if (++colcnt>gifcols) 
			{
				rowcnt++;
				colcnt=1;
			}
		}
//			output << -9999 << '\t';
	}
	else
	{
		if (grid)
		{
			grid[rowcnt][colcnt]=value;
			if (++colcnt>gifcols) 
			{
				rowcnt++;
				colcnt=1;
			}
		}
//			output << value << '\t';
	}
}

char GIFReader::getFirstIndexOfCode(int index)
{
	if (vCodes[index].prevCode>=0)
	{
		return getFirstIndexOfCode(vCodes[index].prevCode);
	}
	else
	{
		return vCodes[index].code;
	}
}

void GIFReader::addNewCode(int prevCode,char code)
{
	codeType ctCode;
	ctCode.prevCode=prevCode;
	ctCode.code=code;
	vCodes.push_back(ctCode);
	if ((vCodes.size()==nextCodeExtensionIndex+1)&&(actCodeLength<12))
	{
		CodeMask = (CodeMask<<1) + 1;
		nextCodeExtensionIndex = (nextCodeExtensionIndex << 1) + 1;
		actCodeLength++;
	}
}

GIFReader::~GIFReader()
{
	if (buffer) delete[] buffer;
	if (colortable) delete[] colortable;
	if (comment) delete[] comment;
	s.close();
}

GIFReader::GIFReader(std::string infilename):
infilename(infilename),buffer(NULL),
colortable(0),comment(0),
gifcols(0),gifrows(0),infofield(0),colortablesize(0)
{
	buffer = new char[3000];
}

//GIFReader::GIFReader(std::string infilename,std::stringstream &output):
//infilename(infilename),output(output),buffer(NULL),
//colortable(0),comment(0),
//gifcols(0),gifrows(0),infofield(0),colortablesize(0)
//{
//	buffer = new char[3000];
//}

void GIFReader::GetDimensions(int &width,int &height)
{
	s.open(infilename.c_str(),std::ios::binary);

	s.read(&(header[0]),6);
	s.read(&(logScreenDescr[0]),7);
	memcpy(&gifcols,&(logScreenDescr[0]),2);
	memcpy(&gifrows,&(logScreenDescr[2]),2);
	width=gifcols;
	height=gifrows;
	s.close();
}


void GIFReader::ReadGIF(double **outputgrid)
{
	grid=outputgrid;

	s.open(infilename.c_str(),std::ios::binary);

	s.read(&(header[0]),6);
	s.read(&(logScreenDescr[0]),7);
	memcpy(&gifcols,&(logScreenDescr[0]),2);
	memcpy(&gifrows,&(logScreenDescr[2]),2);
	memcpy(&infofield,&(logScreenDescr[4]),1);

	if (infofield&0x0080) // there is a global color table, but we do not need the table for reading GIF codes, just need overread it when we know how big it is
	{
		colortablesize= 0x1 << (((infofield & 0x70) >> 4) +1);
		colortable =new char[3*(colortablesize+1)];
		s.read(&(colortable[0]),3*colortablesize);
	}

	// no read the next byte in order to define if it is any extension or an image descriptor
	s.read(&(blockidentifyer),1);
	while (blockidentifyer!=(char)0x3B) // 0x3B would be the termination byte of a GIF file
	{
		if (blockidentifyer==(char)0x21)
		{
			s.read(&(extensionidentifyer),1);
			if (extensionidentifyer==(char)0x01) // Plain Textextension
			{
				char blocksize;
				char *sblock;
				// first block of probably 0x0C bytes will be overread
				s.read(&blocksize,1);
				sblock=new char[blocksize];
				s.read(&(sblock[0]),blocksize);
				delete[] sblock;
				// now subblocks are following
				s.read(&blocksize,1);
				while (blocksize!=(char)0x00)
				{	
					sblock=new char[(unsigned char)blocksize];
					s.read(&(sblock[0]),(unsigned char)blocksize);
					delete[] sblock;
					s.read(&blocksize,1);
				}
			} else if (extensionidentifyer==(char)0xf9) { // Graphics control extension
				char blocksize;
				char blockend;
				char *sblock;
				s.read(&blocksize,1);
				sblock=new char[(unsigned char)blocksize];
				s.read(&(sblock[0]),(unsigned char)blocksize);
				s.read(&blockend,1);
				if (blockend!=0x00)
				{
					std::stringstream errstream; 
					errstream << "Graphics COntrol extension found, but block end not found as 0x00 after "<<blocksize<<" bytes"<<std::endl;
					throw std::string(errstream.str());
				}
				delete[] sblock;
			} else if (extensionidentifyer==(char)0xfe) { // comment
				char blocksize;
				char *sblock;
				// now subblocks are following
				s.read(&blocksize,1);
				while (blocksize!=(char)0x00)
				{	
					sblock=new char[(unsigned char)blocksize];
					s.read(&(sblock[0]),(unsigned char)blocksize);
					delete[] sblock;
					s.read(&blocksize,1);
				}
			} else if (extensionidentifyer==(char)0xff) { // application extension
				char blocksize;
				char *sblock;
				// now subblocks are following. This is similar to comments: each subblock starts with its length indicator and then the indicated number of bytes must be read
				// until a length indicator of 0 is foun, which ends the block
				s.read(&blocksize,1);
				while (blocksize!=(char)0x00)
				{	
					sblock=new char[(unsigned char)blocksize];
					s.read(&(sblock[0]),(unsigned char)blocksize);
					delete[] sblock;
					s.read(&blocksize,1);
				}
			} else {
				// throw an error since we found an unknown etension ID
				unsigned short int code=0x00FF & extensionidentifyer;
				std::stringstream errstream; 
				errstream << "invalid extension type ID \"0x"<<std::hex<<std::uppercase<<code<<"\" found";
				throw std::string(errstream.str());
			}
		} else if (blockidentifyer==(char)0x2C) { // image with leading image descriptor
			// read parameters from image descriptor first - we ignore here left and top and only read for comparison width and height as well as the local color table flag
			unsigned int imgleft=0;
			unsigned int imgtop=0;
			unsigned int imgwidth=0;
			unsigned int imgheight=0;
			unsigned char packedField=0;
			rowcnt=1;
			colcnt=1;
			s.read(&(imgDescr[0]),9);
			
			memcpy(&imgleft,&(imgDescr[0]),2);
			memcpy(&imgtop,&(imgDescr[2]),2);
			memcpy(&imgwidth,&(imgDescr[4]),2);
			memcpy(&imgheight,&(imgDescr[6]),2);
			memcpy(&packedField,&(imgDescr[7]),1);

			if (packedField&0x80) // there is a local color table, but we do not need the table for reading GIF codes, just need overread it when we know how big it is
			{
				colortablesize= 0x1 << (((packedField & 0x70) >> 4) +1);
				if (colortable) delete[] colortable;
				colortable =new char[3*(colortablesize+1)];
				s.read(&(colortable[0]),3*colortablesize);
			}
			// now read the image stream and decompress it
			mincodesize=0;
			s.read(buffer,1);
			memcpy(&mincodesize,&(buffer[0]),1);
			mincodesize++; // one more than in file since clrear code and end of information code need to be coded as well
			char blocksize=0;
			char *sblock=NULL;

			decoderDWord=0;
			tempDWord=0;
			CodeMask=0;
			actCode=0;
			lastCode=0;
			clearCode=colortablesize;
			EOIcode=clearCode+1;

			// int startOffset=0;
			endOffset=0;
			actCodeLength=mincodesize;
			nextCodeExtensionIndex=0; // is the same as CodeMask after initializing the vector vCode later

			s.read(&blocksize,1);
			sblock=new char[(unsigned char)blocksize];
			int bufferPointer=0; // index of the actual valid input byte of sblock stream
			bool justinitialized=false;
			initializeCodeVector(colortablesize);
			justinitialized=true;
			while ((unsigned char)blocksize!=(unsigned char)0x00)
			{
				s.read(&(sblock[0]),(unsigned char)blocksize);
				bufferPointer=0; // reset index of actual valid input byte to 0

				// do something with the block --> LZW decryption must be done here and values written to the internal grid
				while ((unsigned char)bufferPointer<(unsigned char)blocksize)
				{
					while ((endOffset<actCodeLength)&&(unsigned char)bufferPointer<(unsigned char)blocksize)
					{
						tempDWord=0;
						memcpy(&tempDWord,&(sblock[bufferPointer]),1);
						tempDWord = tempDWord<<endOffset;
						decoderDWord = decoderDWord | tempDWord;
						endOffset+=8; // one byte added, so the offset was increased by 8 
						bufferPointer++;
					}
					while (endOffset>=actCodeLength)
					{
						lastCode = actCode; // save actCode in order to be able to create the code table 
						actCode = decoderDWord & CodeMask; // mask decoder word to get the actual code
						decoderDWord = decoderDWord >> actCodeLength; // move on to the next code
						endOffset -=actCodeLength;
						// call color index table or so and write value out to string stream
						if (actCode==clearCode)
						{
							initializeCodeVector(colortablesize);
							justinitialized=true;
						} else if (justinitialized && actCode<vCodes.size()-1) {
							writeOutputToStream(actCode);
							justinitialized=false; // important for first value after re-initialization of code table
						} else if (actCode==EOIcode)	{ // write out what's in the buffer ?
							// there should be nothing in the buffer at all since after every code the value is put to the output stream already.
							// So we need to make sure that both while loops finish here by setting their control variables to specific values that indicate the end of the loop
							endOffset=0;
							bufferPointer=blocksize; // next action will be to read the next blocsize. There shouldn't be any block after an EOI flag (at least not for the actual image)
						} else { // lookUp in Code-Table and write Values to stream and add new Codes to table
							if (actCode<vCodes.size())
							{ // code is present in table
								writeOutputToStream(actCode);
								char K = getFirstIndexOfCode(actCode);
								addNewCode(lastCode,K);
							}
							else // code is not present in table
							{
								char K = getFirstIndexOfCode(lastCode);
								addNewCode(lastCode,K);
								writeOutputToStream(vCodes.size()-1);
							}
						}
					}
				}
				s.read(&blocksize,1); // read the size of the next block (can be smaller, especially for the last block
			}
			delete[] sblock;
		} else	{
			unsigned short int code=0x00FF & blockidentifyer;
			std::stringstream errstream; 
			errstream << "invalid extension type ID \"0x"<<std::hex<<std::uppercase<<code<<"\" found";
			throw std::string(errstream.str());
		}
		s.read(&(blockidentifyer),1);
	}
}

// cutting out a model grid with different origins and different resolution
void GIFReader::CutModelGrid(double **GIFOrig,int mincolsgo,int maxcolsgo,int minrowsgo,int maxrowsgo,double xllgo,double yllgo,double csgo,
							double **ModelGrid,int mincolsmod,int maxcolsmod,int minrowsmod,int maxrowsmod,double xllmod,double yllmod,double csmod,double **maskgrid)
{
	// *go=GIF Original
	// *mod=model
	double ycoord=yllmod;
	double x[5],y[5],z[5],dist[5],u,v;
	int row[5],col[5];
//#pragma omp parallel for private(ycoord,x,y,z,u,v,row,col)
	for (int i=maxrowsmod;i>=minrowsmod;i--)
	{
		double xcoord=xllmod;
		for (int j=mincolsmod;j<=maxcolsmod;j++)
		{
			if (maskgrid && maskgrid[i][j]<=-9999)
			{
				ModelGrid[i][j]=maskgrid[i][j];
			}
			else
			{ // interpolating from GIFOrig
				// knwon values: ccordinates of the model grid cell: xcoord and ycoord
				// what we need: x- and y- coordinates of the GIFGrid of the 4 surrounding points and the values for that points
				bool nodata=false;
				row[1]=row[2]=maxrowsgo-floor((ycoord-(yllgo+csgo/2.0))/csgo);
				row[3]=row[4]=row[1]-1;
				col[1]=col[3]=mincolsgo+floor((xcoord-(xllgo+csgo/2.0))/csgo);
				col[2]=col[4]=col[1]+1;
				if (row[1]>maxrowsgo) { row[1]=row[2]=maxrowsgo; } else	if (row[1]<minrowsgo) {row[1]=row[2]=minrowsgo; }
				if (row[3]>maxrowsgo) { row[3]=row[4]=maxrowsgo; } else	if (row[3]<minrowsgo) {row[3]=row[4]=minrowsgo; }
				if (col[1]>maxcolsgo) { col[1]=col[3]=maxcolsgo; } else if (col[1]<mincolsgo) {col[1]=col[3]=mincolsgo; }
				if (col[2]>maxcolsgo) { col[2]=col[4]=maxcolsgo; } else if (col[2]<mincolsgo) {col[2]=col[4]=mincolsgo; }
				y[1]=y[2]=yllgo+(maxrowsgo-row[1]+0.5)*csgo;
				y[3]=y[4]=y[1]+csgo;
				x[1]=x[3]=xllgo+(col[1]-mincolsgo+0.5)*csgo;
				x[2]=x[4]=x[1]+csgo;
				for (int k=1;k<=4;k++)
				{
					z[k]=GIFOrig[row[k]][col[k]];
					if (z[k]<=-9999.0) nodata=true;
				}
				if (!nodata)
				{
					u=(xcoord-x[1])/csgo;
					v=(ycoord-y[1])/csgo;
					ModelGrid[i][j]=
						(1.0-u)*(1.0-v)*z[1]+
						(    u)*(1.0-v)*z[2]+
						(1.0-u)*(    v)*z[3]+
						(    u)*(    v)*z[4];
				}
				else
				{
					if ((z[1] != -9999) || (z[2] != -9999) || (z[3] != -9999) || (z[4] != -9999))
					{ // switch to IDW
						double dx, dy, sumdist=0;
						bool exactmatch = false;
						for (int k=1;k<=4;k++)
							//      for (pe=&(entf[nr][1]);pe<=&(entf[nr][(anzahl[nr])]);pe++)
						{
							if (z[k] >= -9999)
							{
								dx = (x[k] - xcoord);
								dy = (y[k] - ycoord);
								dist[k] = dx * dx + dy * dy;
								if (dist[k] <= 0.0)
								{
									ModelGrid[i][j] = z[k];
									exactmatch = true;
									break;
								}
								dist[k] = sqrt(dist[k]);
								sumdist += dist[k]; // weight is 1
								if (ModelGrid[i][j] == -9999)
								{
									ModelGrid[i][j] = z[k] * dist[k];
								}
								else
								{
									ModelGrid[i][j] += z[k] * dist[k];
								}
							}
						}
						if (!exactmatch)
						{
							ModelGrid[i][j] /= sumdist;
						}
					}
					else
					{
						ModelGrid[i][j] = -9999;
					}
				}
			}
			xcoord+=csmod;
		}
		ycoord+=csmod;
	}
}

void GIFReader::MapValues(double **ModelGrid, int mincolsmod, int maxcolsmod, int minrowsmod, int maxrowsmod, double f1, double x1, int cut, double f2, double x2)
{
//	double f1 = 8.916787e-5; // scaling parameter
//	double x1 = 0.0767528404; // exponent scaling
//	int cut = 127; // the new params f2 and x2 are valid for 135 and higher
//	double f2 = 6.00E-006; // scaling parameter for higher precipitations above ~2mm
//	double x2 = 0.098; // exponent scaling for higher precipitations above ~2mm

	for (int i=minrowsmod;i<=maxrowsmod;i++)
	{
		for (int j=mincolsmod;j<=maxcolsmod;j++)
		{
			if (ModelGrid[i][j]>-9999)
			{
				if (ModelGrid[i][j]==0) {
					ModelGrid[i][j]=0;
				}
				else if (ModelGrid[i][j] <= cut) {
					ModelGrid[i][j] = f1*exp(x1*ModelGrid[i][j]);
				}
				else if (ModelGrid[i][j] <= 250) {
					ModelGrid[i][j] = f2*exp(x2*ModelGrid[i][j]);
				}
				else {
					ModelGrid[i][j] = -9999;
				}
			}
		}
	}
}
