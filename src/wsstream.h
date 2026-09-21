#ifndef WSSTREAM
#define WSSTREAM

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include <string>

// to simulate setw for ostreams, we declare a class setw. Each time, setw is inserted into a stream, 
// an object of this class is created and passed to the wsstream-object. There, an inserter operator with 
// argument setw& executes the inserter command for setw on the underlying ostream (the argument given to the execute command in the wsstream class)

namespace wasim {

	class setw
	{
		std::streamsize _Param;
	public:

		setw(std::streamsize pParam) :_Param(pParam) {}
		~setw() {}

		template<typename T>
		void Execute(T& pTarget)
		{
			pTarget << std::setw(_Param);
		}
	};

	class setfill
	{
		char _Param;
	public:

		setfill(char pParam) :_Param(pParam) {}
		~setfill() {}

		template<typename T>
		void Execute(T& pTarget)
		{
			pTarget << std::setfill(_Param);
		}
	};

	class mystream
	{
	public:

		int mpi_rank;
		int mpi_size;

		std::ofstream os;
		std::ostream *ssc;

		bool newdiag;
		bool suppress;
		int level;
		std::string internalstr;
		std::string sOutputfilename;
		int iDefaultLevel;


		mystream(std::ostream *extstr,std::string outputfilename,int defaultlevel=4);
		mystream(std::string outputfilename,int defaultlevel);
		mystream();
		~mystream();

		// Manipulators
		/*	mystream& endlist(mystream& _Ostr);
	mystream& endl(mystream& _Ostr);
	mystream& setw(std::streamsize width);
	mystream& uppercase(mystream& _Ostr);
*/
		void os_manipulations(char *str);

		mystream& operator<<(mystream& (*_Pfn)(mystream&));
		mystream& operator<<(std::basic_ios<char,std::char_traits<char> >& (*_Pfn)(std::basic_ios<char,std::char_traits<char> >&));
		mystream& operator<<(std::ios_base& (*_Pfn)(std::ios_base&));
		mystream& operator<<(mystream& _str);
		mystream& operator<<(const char _Val);
		mystream& operator<<(const char *_Val);
		mystream& operator<<(const int _Val);
		mystream& operator<<(const unsigned int _Val);
		mystream& operator<<(const long _Val);
		mystream& operator<<(const unsigned long _Val);
		mystream& operator<<(const float _Val);
		mystream& operator<<(const double _Val);
		mystream& operator<<(const std::string _Val);
		mystream& operator<<(const void *_Val);
		mystream& operator<<(const void **_Val);
		//	mystream& operator<<(const std::_Smanip<std::ios_base::fmtflags>& _Manip);
		mystream& operator<<(wasim::setw Manip);
		mystream& operator<<(wasim::setfill Manip);
	};

}
#endif
