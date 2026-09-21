#ifndef CSPLINEH
#define CSPLINEH
#include "modef.h"

//   (c) 1999-2007 by Joerg Schulla
//   more version history: see source safe comments


class TSplines:
public memory_handling
{
public:

	double **yNRAbyNCE,**ycNRAbyNCE,**xcNRAbyNCE,**tempgrid; /* 2nd derivatives in y direction for all output row y values for each input column */
	double **CoordXe,**CoordYe,**CoordXa,**CoordYa,**zn;
	double *x,*y2,*xcol,*y2col,finefact,beta;
	double *y,*ycol;
	double xtempmin,xtempmax,ytempmin,ytempmax,xlltemp,ylltemp,xlltemprot,ylltemprot;
	long   nrtempmin,nrtempmax,nctempmin,nctempmax,nrtemp,nctemp;
	long   equalarea,nrowsein,ncolsein,nrowsaus,ncolsaus;
	double xllein,yllein,cseinx,cseiny,nodataein,xllaus,yllaus,csaus,nodataaus;
	long i,j,k,m,klo,khi,anf,end,r,r1,r2,c,rowindex,colindex;
	double rf,tx1,ty1,shy,shx;
	double CYa,CXe,CXa,cosbeta,sinbeta,yoffset,xlleinrot,ylleinrot,dx,dy,entf,alpha,gamma;
	std::FILE *tst;
	unsigned short **LookUpTable;

	/* Constructor */
	// int equalarea    :  if set to 1, then the last value in a row in repeated in ordwer to allow interpolation within single gridcells
	TSplines(double **CoordXe,double **CoordYe,double **CoordXa,double **CoordYa,
	long tnrowsein,long tncolsein,double txllein,double tyllein,double tcseinx,double tcseiny,double tnodataein,
	long tnrowsaus,long tncolsaus,double txllaus,double tyllaus,double tcsaus,double tnodataaus,
	double tbeta,double tfinefact,long tequalarea,double **tzn);
	/* Destructor */
	~TSplines();

	virtual void RotateCoordinates(double *x,double *y,double beta);

	// function ConstructCubicSpline
	// to construct the parameters for a 1-D-Cubic Spline (natural cubic spline if yp1 and ypn are >=1e30)
	// version 1, 27.7.1999, Joerg Schulla
	// parameters: 
	// double *x         :  input x values as a vector [0..n+1]
	// double *y         :  input y values as a vector [0..n+1]
	// int anf          :  first index of tabulated values
	// int end          :  last index of tabulated values
	// double yp1 and ypn:  initial boundary conditions for the second derivatives, if >=1e30
	//                     then the boundary conditions are set to 0 (natural spline)
	// double *y2        :  result of the constructed spline [0..n+1](second derivatives at the given coordinates in *x) 
	virtual void ConstructCubicSpline(double *x,double *y,long anf,long end,double yp1,double ypn,double *y2);

	// function EvalCubicSpline
	// to evaluate a 1-D-Spline (get the real interpolation results)
	// version 1, 28.7.1999, Joerg Schulla
	// parameters: 
	// double *x         :  input x values as a vector [0..n_ein+1], arranged in ascending order!!!
	// double *y         :  input y values as a vector [0..n_ein+1]  (must match the same order as the xa's!!!)
	// double *y2        :  input y2 (2nd derivatives) as a vector [0..n_ein+1] as result of ConstructCubicSpline
	// int n            :  number of tabulated input values
	// double xa         :  output-grid x value 
	// double ya         :  result of the spline evaluation (interpolation results for the coordinates given in xa parameter) 
	virtual void EvalCubicSpline(double *x,double *y,long *klo,long *khi,double *y2,long n,double xa,double *ya,double nodataein,double nodataaus);
	virtual void BicubicColumnSpline(double **ze,double **y2a);
	virtual void BicubicRowSpline(double **ze,double **y2a,double **za);
	virtual void GetValidSplineArea();
	virtual void ReRotateResults(double **za);
	virtual void LookUpTableInit();
};
#endif
