#include <cstdio>
#include <ctime>
#include <sstream>
#include "wsstream.h"

#ifdef _MPI
#include <mpi.h>
#endif

//#define _MPIDEBUG

using namespace wasim;

//mystream::mystream():std::ostream(NULL)
mystream::mystream():
newdiag(true),
level(0),
suppress(false),
iDefaultLevel(4),
mpi_rank(0),
mpi_size(1)
{
	ssc = &std::cout;
	this->sOutputfilename="WaSiM_Diag.xml";
#ifdef _MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	std::stringstream ss;
	ss << sOutputfilename.substr(0,sOutputfilename.rfind(".")) << "_" << mpi_rank << ".xml";
	sOutputfilename = ss.str();
#endif
#ifndef _MPIDEBUG
	if (mpi_rank>0) return;
#endif	
	os.open(sOutputfilename.c_str());
	if (os.good())
	{
		os << "<?xml version='1.0' encoding='UTF-8'?>\n"
		<< "<Diag version=\"1.2\" xmlns=\"http://www.wldelft.nl/fews/PI\" \n"
		<< "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \n"
		<< "xsi:schemaLocation=\"http://www.wldelft.nl/fews/PI http://fews.wldelft.nl/schemas/version1.0/pi-schemas/pi_diag.xsd\">\n";
		time_t ltime;
		time(&ltime);
		struct tm* loctime;
		loctime=localtime(&ltime);
		std::stringstream ss;
		std::string s;
		ss << std::setw(4) << std::setfill('0')
		<< loctime->tm_year+1900 << "-" << std::setw(2) 
		<< std::setw(2) << loctime->tm_mon+1  << "-" 
		<< std::setw(2) << loctime->tm_mday << "T" 
		<< std::setw(2) << loctime->tm_hour << ":" 
		<< std::setw(2) << loctime->tm_min  << ":" 
		<< std::setw(2) << loctime->tm_sec;
		ss >> s;
		os << "\t<line level=\"3\" description= \"WaSiM start at " << s << "\"/>\n";			
	}
	else
	{
		(*ssc) << "(mpi-rank "<<std::setw(2)<<mpi_rank<<") "<<"Error opening diagnostics output file " << sOutputfilename << std::endl; 
		throw 320;
	}
}

// for use within wasim itself, an outputfile is given as a parameter as -diagfile:<filenam> commandline parameter
mystream::mystream(std::string outputfilename,int defaultlevel):
newdiag(true),
level(0),
suppress(false),
iDefaultLevel(defaultlevel),
mpi_rank(0),
mpi_size(1)
{
	ssc = &std::cout;
	
	sOutputfilename=outputfilename;
	if (this->sOutputfilename.size()==0)
	{
		sOutputfilename="WaSiM_Diag.xml";
	}
#ifdef _MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	std::stringstream ss;
	ss << sOutputfilename.substr(0,sOutputfilename.rfind(".")) << "_" << mpi_rank << ".xml";
	sOutputfilename = ss.str();
#endif
#ifndef _MPIDEBUG
	if (mpi_rank>0) return;
#endif	
	os.open(sOutputfilename.c_str());
	if (!os.good())
	{
		os.clear();
		(*ssc) << "(mpi-rank " << std::setw(2) << mpi_rank << ") " << "\n ----- Warning -----\n error opening diagnostics output file " << sOutputfilename << std::endl;
		(*ssc) << "(mpi-rank " << std::setw(2) << mpi_rank << ") " << "\nusing default name WaSiM_Diag.xml in actual working directory\n -------------------\n\n" << std::endl;
		sOutputfilename="WaSiM_Diag.xml";
		os.open(sOutputfilename.c_str());
	}

	if (os.good())
	{
		os << "<?xml version='1.0' encoding='UTF-8'?>\n"
		<< "<Diag version=\"1.2\" xmlns=\"http://www.wldelft.nl/fews/PI\" \n"
		<< "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \n"
		<< "xsi:schemaLocation=\"http://www.wldelft.nl/fews/PI http://fews.wldelft.nl/schemas/version1.0/pi-schemas/pi_diag.xsd\">\n";
		time_t ltime;
		time(&ltime);
		struct tm* loctime;
		loctime=localtime(&ltime);
		std::stringstream ss;
		std::string s;
		ss << std::setw(4) << std::setfill('0')
		<< loctime->tm_year+1900 << "-" << std::setw(2) 
		<< std::setw(2) << loctime->tm_mon+1  << "-" 
		<< std::setw(2) << loctime->tm_mday << "T" 
		<< std::setw(2) << loctime->tm_hour << ":" 
		<< std::setw(2) << loctime->tm_min  << ":" 
		<< std::setw(2) << loctime->tm_sec;
		ss >> s;
		os << "\t<line level=\"3\" description= \"WaSiM start at " << s << "\"/>\n";			
	}
	else
	{
		(*ssc) << "(mpi-rank " << std::setw(2) << mpi_rank << ") " << "Error opening diagnostics output file " << sOutputfilename << std::endl;
		throw 320;
	}
}

// for use in the dll: the stream is some buffer for communication with the external caller
mystream::mystream(std::ostream *extstr,std::string outputfilename,int defaultlevel):
newdiag(true),
level(0),
suppress(false),
iDefaultLevel(defaultlevel)
{
	if (mpi_rank>0) return;
	ssc = extstr;
	sOutputfilename=outputfilename;
	if (this->sOutputfilename.size()==0)
	sOutputfilename="WaSiM_Diag.xml";
#ifdef _MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
//	std::stringstream ss;
//	ss << sOutputfilename.substr(0,sOutputfilename.rfind(".")) << "_" << mpi_rank << ".xml";
//	sOutputfilename = ss.str();
#endif
#ifndef _MPIDEBUG
	if (mpi_rank>0) return;
#endif	
	os.open(sOutputfilename.c_str());
	if (!os.good())
	{
		os.clear();
		(*ssc) << "(mpi-rank " << std::setw(2) << mpi_rank << ") " << "\n ----- Warning -----\n error opening diagnostics output file " << sOutputfilename << std::endl;
		(*ssc) << "(mpi-rank " << std::setw(2) << mpi_rank << ") " << "\nusing default name WaSiM_Diag.xml in actual working directory\n -------------------\n\n" << std::endl;
		sOutputfilename="WaSiM_Diag.xml";
		os.open(sOutputfilename.c_str());
	}
	if (os.good())
	{
		os << "<?xml version='1.0' encoding='UTF-8'?>\n"
		<< "<Diag version=\"1.2\" xmlns=\"http://www.wldelft.nl/fews/PI\" \n"
		<< "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \n"
		<< "xsi:schemaLocation=\"http://www.wldelft.nl/fews/PI http://fews.wldelft.nl/schemas/version1.0/pi-schemas/pi_diag.xsd\">\n";
		time_t ltime;
		time(&ltime);
		struct tm* loctime;
		loctime=localtime(&ltime);
		std::stringstream ss;
		std::string s;
		ss << std::setw(4) << std::setfill('0')
		<< loctime->tm_year+1900 << "-" << std::setw(2) 
		<< std::setw(2) << loctime->tm_mon+1  << "-" 
		<< std::setw(2) << loctime->tm_mday << "T" 
		<< std::setw(2) << loctime->tm_hour << ":" 
		<< std::setw(2) << loctime->tm_min  << ":" 
		<< std::setw(2) << loctime->tm_sec;
		ss >> s;
		os << "\t<line level=\"3\" description= \"WaSiM start at " << s << "\"/>\n";			
	}
	else
	{
		(*ssc) << "(mpi-rank " << std::setw(2) << mpi_rank << ") " << "Error opening diagnostics output file " << sOutputfilename << std::endl;
		throw 320;
	}
}

mystream::~mystream()
{
	if (os.good()) 
	{
		if (!newdiag)
		{
			// closing the open statement first
			os << "\"/>\n";
		}
		time_t ltime;
		time(&ltime);
		struct tm* loctime;
		loctime=localtime(&ltime);
		std::stringstream ss;
		std::string s;
		ss << std::setw(4) << std::setfill('0')
		<< loctime->tm_year+1900 << "-" << std::setw(2) 
		<< std::setw(2) << loctime->tm_mon+1  << "-" 
		<< std::setw(2) << loctime->tm_mday << "T" 
		<< std::setw(2) << loctime->tm_hour << ":" 
		<< std::setw(2) << loctime->tm_min  << ":" 
		<< std::setw(2) << loctime->tm_sec;
		ss >> s;
		os << "\t<line level=\"3\" description= \"WaSiM finished at " << s << "\"/>\n";			
		os << "</Diag>\n";	
		os.close();
	}
}

/*
mystream& mystream::endlist(mystream& _Ostr)
{
	(*ssc) << std::endl;
	// closing the line-tag
	if (os.good()&&(suppress==false)) 
	{
		if (newdiag==false)
		{
			os << "\"/>";
		}
		os << std::endl;
		newdiag=true;
	}
	suppress=false;
	return (*this);
}

mystream& mystream::endl(mystream& _Ostr)
{
	(*ssc) << std::endl;
	// closing the line-tag
	if (os.good()&&(suppress==false)) 
	{
		if (newdiag==false)
		{
			os << "\"/>";
		}
		os << std::endl;
		newdiag=true;
	}
	suppress=false;
	return (*this);
}

mystream& mystream::setw(std::streamsize width)
{
	(*ssc) << std::setw(width);
	if (os.good()&&(suppress==false)) os << std::setw(width);
	return (*this);
}

mystream& mystream::uppercase(mystream& _Ostr)
{
	(*ssc).uppercase;
	if (os.good()) os.uppercase;
	return (*this);
}
*/


mystream& mystream::operator<<(mystream& (*_Pfn)(mystream&))
{	// call basic_ostream manipulator
	return ((*_Pfn)(*this));
}

mystream& mystream::operator<<(std::basic_ios<char,std::char_traits<char> >& (*_Pfn)(std::basic_ios<char,std::char_traits<char> >&))
{	// call basic_ios manipulator
	(*_Pfn)(*(std::basic_ios<char,std::char_traits<char> >*)this);
	return (*this);
}

mystream& mystream::operator<<(std::ios_base& (*_Pfn)(std::ios_base&))
{	// call ios_base manipulator
	(*_Pfn)(*(std::ios_base *)this);
	return (*this);
}

mystream& mystream::operator<<(mystream& _str)
{	// call basic_ostream manipulator
	return (_str);
}

mystream& mystream::operator<<(const char _Val)
{
#ifndef _MPIDEBUG
	if (mpi_rank==0) (*ssc) << _Val;
#else	
	(*ssc) << _Val;
#endif	
	os_manipulations(NULL);
#ifndef _MPIDEBUG
	if (mpi_rank==0) if (os.good()&&(suppress==false)) os << _Val;
#else	
	if (os.good()&&(suppress==false)) os << _Val;
#endif	
	return *this;
}

mystream& mystream::operator<<(const char *_Val)
{
	std::string messg(_Val);
	bool isend=false;
	if (messg.find("[end-entry]")!=std::string::npos)
	{
		messg=messg.substr(0,messg.length()-11);
		isend=true;
	}
	if (messg.find("SKIP: ")==0) {
		messg=messg.substr(6);
	} else if (messg.find("DEBUG: ")==0) {
		messg=messg.substr(7);
	} else if (messg.find("INFO: ")==0) {
		messg=messg.substr(6);
	} else if (messg.find("WARNING: ")==0) {
		messg=messg.substr(9);
	} else if (messg.find("ERROR: ")==0) {
		messg=messg.substr(7);
	} else if (messg.find("FATAL: ")==0) {
		messg=messg.substr(7);
	}
	if (messg.size()>0)
	{
#ifndef _MPIDEBUG
		if (mpi_rank==0) (*ssc) << messg;
#else		
		(*ssc) << messg;
#endif
	}
	os_manipulations((char*)_Val);
#ifndef _MPIDEBUG
	if (mpi_rank==0) if (os.good()&&(suppress==false)) os << internalstr;
#else		
	if (os.good()&&(suppress==false)) os << internalstr;
#endif
	if (isend)
	{

#ifdef _MPI
#ifndef _MPIDEBUG
		if (mpi_rank==0) (*ssc) << std::endl << "(mpi-rank " << std::setw(2) << mpi_rank << ") ";
#else		
		(*ssc) << std::endl << "(mpi-rank " << std::setw(2) << mpi_rank << ") ";
#endif
#else
		(*ssc) << std::endl;
#endif
		if (!suppress) 
		{
#ifndef _MPIDEBUG
			if (mpi_rank==0) os << std::endl;
#else		
			os << std::endl;
#endif
		}
		this->newdiag=true;
		this->suppress=false;
	}
	return *this;
}

mystream& mystream::operator<<(const int _Val)
{
#ifndef _MPIDEBUG
	if (mpi_rank==0) (*ssc) << _Val;
#else
	(*ssc) << _Val;
#endif	
	os_manipulations(NULL);
#ifndef _MPIDEBUG
	if (mpi_rank==0) if (os.good()&&(suppress==false)) os << _Val;
#else
	if (os.good()&&(suppress==false)) os << _Val;
#endif	
	return *this;
}

mystream& mystream::operator<<(const unsigned int _Val)
{
#ifndef _MPIDEBUG
	if (mpi_rank==0) (*ssc) << _Val;
#else
	(*ssc) << _Val;
#endif	
	os_manipulations(NULL);
#ifndef _MPIDEBUG
	if (mpi_rank==0) if (os.good()&&(suppress==false)) os << _Val;
#else
	if (os.good()&&(suppress==false)) os << _Val;
#endif	
	return *this;
}

mystream& mystream::operator<<(	const long _Val)
{
#ifndef _MPIDEBUG
	if (mpi_rank==0) (*ssc) << _Val;
#else
	(*ssc) << _Val;
#endif	
	os_manipulations(NULL);
#ifndef _MPIDEBUG
	if (mpi_rank==0) if (os.good()&&(suppress==false)) os << _Val;
#else
	if (os.good()&&(suppress==false)) os << _Val;
#endif	
	return *this;
}

mystream& mystream::operator<<(	const unsigned long _Val)
{
#ifndef _MPIDEBUG
	if (mpi_rank==0) (*ssc) << _Val;
#else
	(*ssc) << _Val;
#endif	
	os_manipulations(NULL);
#ifndef _MPIDEBUG
	if (mpi_rank==0) if (os.good()&&(suppress==false)) os << _Val;
#else
	if (os.good()&&(suppress==false)) os << _Val;
#endif	
	return *this;
}

mystream& mystream::operator<<(	const float _Val)
{
#ifndef _MPIDEBUG
	if (mpi_rank==0) (*ssc) << _Val;
#else
	(*ssc) << _Val;
#endif	
	os_manipulations(NULL);
#ifndef _MPIDEBUG
	if (mpi_rank==0) if (os.good()&&(suppress==false)) os << _Val;
#else
	if (os.good()&&(suppress==false)) os << _Val;
#endif	
	return *this;
}

mystream& mystream::operator<<(	const double _Val)
{
#ifndef _MPIDEBUG
	if (mpi_rank==0) (*ssc) << _Val;
#else
	(*ssc) << _Val;
#endif	
	os_manipulations(NULL);
#ifndef _MPIDEBUG
	if (mpi_rank==0) if (os.good()&&(suppress==false)) os << _Val;
#else
	if (os.good()&&(suppress==false)) os << _Val;
#endif	
	return *this;
}

mystream& mystream::operator<<(const std::string _Val)
{
#ifndef _MPIDEBUG
	if (mpi_rank==0) (*ssc) << _Val;
#else
	(*ssc) << _Val;
#endif	
	os_manipulations((char*)_Val.c_str());
#ifndef _MPIDEBUG
	if (mpi_rank==0) if (os.good()&&(suppress==false)) os << _Val;
#else
	if (os.good()&&(suppress==false)) os << _Val;
#endif	
	return *this;
}

mystream& mystream::operator<<(const void *_Val)
{
#ifndef _MPIDEBUG
	if (mpi_rank==0) (*ssc) << _Val;
#else
	(*ssc) << _Val;
#endif	
	os_manipulations(NULL);
#ifndef _MPIDEBUG
	if (mpi_rank==0) if (os.good()&&(suppress==false)) os << _Val;
#else
	if (os.good()&&(suppress==false)) os << _Val;
#endif	
	return *this;
}

mystream& mystream::operator<<(const void **_Val)
{
#ifndef _MPIDEBUG
	if (mpi_rank==0) (*ssc) << _Val;
#else
	(*ssc) << _Val;
#endif	
	os_manipulations(NULL);
#ifndef _MPIDEBUG
	if (mpi_rank==0) if (os.good()&&(suppress==false)) os << _Val;
#else
	if (os.good()&&(suppress==false)) os << _Val;
#endif	
	return *this;
}

//mystream& mystream::operator<<(const std::_Smanip<std::ios_base::fmtflags>& _Manip)
mystream& mystream::operator<<(wasim::setw Manip)
{	// insert by calling function with output stream and argument
	//	(*ssc) << Manip;
#ifndef _MPIDEBUG
	if (mpi_rank == 0) Manip.Execute(*ssc);
#else
	Manip.Execute(*ssc);
#endif	
	os_manipulations(NULL);
	//	if (os.good()&&(suppress==false)) os << Manip;
#ifndef _MPIDEBUG
	if (mpi_rank == 0) if (os.good() && (suppress == false)) Manip.Execute(os);
#else
	if (os.good() && (suppress == false)) Manip.Execute(os);
#endif	
	return (*this);
}

mystream& mystream::operator<<(wasim::setfill Manip)
{	// insert by calling function with output stream and argument
	//	(*ssc) << Manip;
#ifndef _MPIDEBUG
	if (mpi_rank == 0) Manip.Execute(*ssc);
#else
	Manip.Execute(*ssc);
#endif	
	os_manipulations(NULL);
	//	if (os.good()&&(suppress==false)) os << Manip;
#ifndef _MPIDEBUG
	if (mpi_rank == 0) if (os.good() && (suppress == false)) Manip.Execute(os);
#else
	if (os.good() && (suppress == false)) Manip.Execute(os);
#endif	
	return (*this);
}

void mystream::os_manipulations(char *str)
{
	time_t ltime;
	time(&ltime);
	bool closefirst=false;
	std::string messg("");
	if (str)
	{
		messg = std::string(str);
		if (messg.find("FATAL: ")==0) {
			level=0;
			if (newdiag==false) 
			{
				closefirst=true;
				newdiag=true;
			}
			messg=messg.substr(7);
		} else if (messg.find("ERROR: ")==0) {
			level=1;
			if (newdiag==false)
			{
				closefirst=true;
				newdiag=true;
			}
			messg=messg.substr(7);
		} else if (messg.find("WARNING: ")==0) {
			level=2;
			if (newdiag==false) 
			{
				closefirst=true;
				newdiag=true;
			}
			messg=messg.substr(9);
		} else if (messg.find("INFO: ")==0) {
			level=3;
			if (newdiag==false) 
			{
				closefirst=true;
				newdiag=true;
			}
			messg=messg.substr(6);
		} else if (messg.find("DEBUG: ")==0) {
			level=4;
			if (newdiag==false) 
			{
				closefirst=true;
				newdiag=true;
			}
			messg=messg.substr(7);
		} else if (messg.find("SKIP: ")==0) {
			level=5;
			if (newdiag==false) 
			{
				closefirst=true;
				newdiag=true;
			}
			messg=messg.substr(6);
			suppress=true;
		} else {
			if (newdiag==true)
			{
				level=4; // 4 = DEBUG, set only if this is a new entry and no level was set because no level was given
			}
		}
	}
	else
	{
		if (newdiag==true)
		{
			level=4; // 4 = DEBUG, set only if this is a new entry and no level was set because no level was given
		}
	}
	if (closefirst && (!suppress)) // if the old tag may be open, close it, except when the old tag was suppressed
	{
		//		if (os.good()&&(suppress==false)) os << "\"/>" << std::endl;
#ifndef _MPIDEBUG
		if (mpi_rank==0) if (os.good()) os << "\"/>" << std::endl;
#else
		if (os.good()) os << "\"/>" << std::endl;
#endif	
		
		newdiag=true;
	}

	if (level>iDefaultLevel)
	suppress=true;
	else
	suppress=false;

	if (newdiag)
	{
		if (os.good()&&(suppress==false)) 
		{
			//			os << "<line level=\"" << level <<"\" description= \"["<< ctime(&ltime) << "]: ";			
#ifndef _MPIDEBUG
			if (mpi_rank==0) os << "\t<line level=\"" << level <<"\" description= \"";			
#else
			os << "\t<line level=\"" << level <<"\" description= \"";			
#endif		
		}
		newdiag=false;
	}

	// replacing non-XML conform substrings
	while (messg.find("&")!=std::string::npos)
	{
		messg.replace(messg.find("'"),1,"&amp;");
	}

	while (messg.find("\n")!=std::string::npos)
	{
		messg.replace(messg.find("\n"),1,";");
	}

	while (messg.find("<")!=std::string::npos)
	{
		messg.replace(messg.find("<"),1,"&lt;");
	}

	while (messg.find(">")!=std::string::npos)
	{
		messg.replace(messg.find(">"),1,"&gt;");
	}

	while (messg.find("'")!=std::string::npos)
	{
		messg.replace(messg.find("'"),1,"&apos;");
	}

	while (messg.find("\"")!=std::string::npos)
	{
		messg.replace(messg.find("\""),1,"&quot;");
	}

	while (messg.find("---")!=std::string::npos)
	{
		messg.replace(messg.find("---"),1,"");
	}
	// replacing end-tag
	if (messg.find("[end-entry]")!=std::string::npos)
	{
		messg=messg.substr(0,messg.length()-11)+std::string("\"/>");
		//		newdiag=true;
	}

	internalstr=messg;

}
