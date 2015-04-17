
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <algorithm>

#ifndef SPACE_H
#define SPACE_H
 

using namespace std;

class space{
public:
	int size[3];	
	//-----------------------------------------
	//model contains external layer of cells=>
	//(number of cells by axis) = (number of points in this axis) + 1
	//									
	//size[0]=nx+1, size[1]=ny+1; size[3]=nz+1;	
	//naxis = number of points in this axis
	//
	//cells numbered 
	//(0,0,0)(0,0,1)..(0,0,nz+1)(0,1,0)..(0,1,nz+1)(0,2,0)..(0,ny+1,nz+1)  (1,0,0)      ..  (nx+1,ny+1,nz+1)
	//	 0		1  ..	nz+1	 nz+2  .. 2*nz+1   2*nz+2 .. (ny+1)(nz+1)  (ny+1)(nz+1) ..  
	//------------------------------------------
	int* layers;	
	double *density, *xvelosity, *yvelosity, *zvelosity, *pressure, *source;
	double *xpoints, *ypoints, *zpoints;
	double deltat;
	space(const string& gridfname, const char* option);
	space(const space& s);
	~space();
	void showall(ofstream& logfile) const;
	void savedata(string& datafname);
	void savestate(string& statefname, double& t, double& tmax, int& nstep, int& nstepmax);
	void restore(string& statefname, double& t,  int& nstep);
	void fillspace(const string& ivfname);
	void getcellsize(double* cellsize, int ix, int iy, int iz);
	void getriemparam(double* param, const string& towards, int ix, int iy, int iz);
	void updatetimestep(const string& axis, int ix, double s1, double s3, double& deltat);
	void tanvel(double* Ures, double s2, int tan1, int tan2, const string&  towards, int ix, int iy, int iz);
	void stopheat();
};

double inenergy(double rho, double pressure);

#endif