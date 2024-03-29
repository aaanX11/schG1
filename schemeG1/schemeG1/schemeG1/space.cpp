#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <set>
#include "space.h"


space::space(const space& s){
	int nx, ny, nz;
	nx = size[0] = s.size[0];
	ny = size[1] = s.size[1];
	nz = size[2] = s.size[2];
	layers = new int[nx*ny*nz];
	for(int i = 0; i < nx*ny*nz; i++){
		layers[i] = s.layers[i];
		density[i] = s.density[i];
		xvelosity[i] = s.xvelosity[i];
		yvelosity[i] = s.yvelosity[i];
		zvelosity[i] = s.zvelosity[i];
		pressure[i] = s.pressure[i];
		source[i] = s.source[i];
	}
	for(int i = 0; i < nx+1; i++){		xpoints[i] = s.xpoints[i];	}
	for(int i = 0; i < ny+1; i++){		ypoints[i] = s.ypoints[i];	}
	for(int i = 0; i < nz+1; i++){		zpoints[i] = s.zpoints[i];	}
	deltat = s.deltat;
}

void readgrid(double* coord, const string& gridfname,const int& size,const string& axis){
	double tmp;
	string s;
	stringstream ss;
	ifstream grid;
	grid.open(gridfname.c_str(),std::ifstream::in);
	if(!grid){
		std::cerr<<"function 'readgrid' cannot open file "<<gridfname<<'\n';
		exit(1);
	}
	while(getline(grid, s) != NULL && s[0] != axis[0]){	}
	getline(grid, s);
	getline(grid, s);
	ss.str(s);
	for(int i = 0; i < size; i++){
		ss>>tmp;
		coord[i] = tmp;
	}	
	ss.clear();
	grid.close();
	return;
}

void readgridaxis(const string& gridfname, int* size){
	int n;
	string s;
	stringstream ss;
	ifstream grid(gridfname.c_str());
	if(!grid){
		std::cerr<<"function 'readgridaxis' cannot open file "<<gridfname<<'\n';
		exit(1);
	}
	while(getline(grid, s) != NULL && s[0] != 'X'){	}
	getline(grid, s);
	ss.str(s);
	ss>>n;
	size[0] = n+2;
	ss.clear();
	while(getline(grid, s) != NULL && s[0] != 'Y'){	}
	getline(grid, s);
	ss.str(s);
	ss>>n;
	size[1] = n+2;
	ss.clear();
	while(getline(grid, s) != NULL && s[0] != 'Z'){	}
	getline(grid, s);
	ss.str(s);
	ss>>n;
	size[2] = n+2;
	ss.clear();
	grid.close();
	return;
}

space::~space(){
	delete[]layers;
	delete[]density;
	delete[]xvelosity;
	delete[]yvelosity;
	delete[]zvelosity;
	delete[]pressure;
	delete[]source;
	delete[]xpoints;
	delete[]ypoints;
	delete[]zpoints;
}
space::space(const string& gridfname, const char* option){
	int tsize[3];
	std::cerr<<"creating space\n";
	readgridaxis(gridfname, tsize);
	if(tsize[0] > 1000 || tsize[0] < 1 || tsize[1] > 1000 || tsize[1] < 1 || tsize[2] > 1000 || tsize[2] < 1){exit(1);}
	std::cerr<<"space size: "<<tsize[0]<<" "<<tsize[1]<<" "<<tsize[2]<<'\n';
	size[0] = tsize[0]; size[1] = tsize[1]; size[2] = tsize[2];
	xpoints = new double[size[0]-1];
	ypoints = new double[size[1]-1];
	zpoints = new double[size[2]-1];
	deltat = 0.0;
	readgrid(xpoints, gridfname, size[0]-1, "X");
	readgrid(ypoints, gridfname, size[1]-1, "Y");
	readgrid(zpoints, gridfname, size[2]-1, "Z");
	layers = new int[size[0]*size[1]*size[2]];
	density = new double[size[0]*size[1]*size[2]];
	xvelosity = new double[size[0]*size[1]*size[2]];
	yvelosity = new double[size[0]*size[1]*size[2]];
	zvelosity = new double[size[0]*size[1]*size[2]];
	pressure = new double[size[0]*size[1]*size[2]];	
	source = new double[size[0]*size[1]*size[2]];
	for(int i = 0; i < size[0]*size[1]*size[2]; i++){
		layers[i] = -1;
		density[i] = 0.0;
		xvelosity[i] = 0.0;
		yvelosity[i] = 0.0;
		zvelosity[i] = 0.0;
		pressure[i] = 0.0;
		source[i] = 0;
	}
	return;
}
void space::fillspace(const string& ivfname){
	std::string s;
	stringstream ss;
	int nx, ny, nz, nny, nnx;
	nx = size[0];
	ny = size[1];
	nz = size[2];
	nny = nz;
	nnx = nz*ny;
	std::ifstream ivfile(ivfname.c_str());
	if(!ivfile){std::cerr<<"can't read "<<ivfname<<'\n'; exit(1);}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			getline(ivfile, s);
			ss.str(s);
			for(int iz = 0; iz < nz; iz++){
				ss>>layers[iz + iy*nny + ix*nnx];
				ss.clear();			}
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			getline(ivfile, s);
			ss.str(s);
			for(int iz = 0; iz < nz; iz++){
				ss>>density[iz + iy*nny + ix*nnx];
				ss.clear();			}
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			getline(ivfile, s);
			ss.str(s);
			for(int iz = 0; iz < nz; iz++){
				ss>>xvelosity[iz + iy*nny + ix*nnx];
				ss.clear();			}
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			getline(ivfile, s);
			ss.str(s);
			for(int iz = 0; iz < nz; iz++){
				ss>>yvelosity[iz + iy*nny + ix*nnx];
				ss.clear();			}
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			getline(ivfile, s);
			ss.str(s);
			for(int iz = 0; iz < nz; iz++){
				ss>>zvelosity[iz + iy*nny + ix*nnx];
				ss.clear();			}
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			getline(ivfile, s);
			ss.str(s);
			for(int iz = 0; iz < nz; iz++){
				ss>>pressure[iz + iy*nny + ix*nnx];
				ss.clear();			}
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			getline(ivfile, s);
			ss.str(s);
			for(int iz = 0; iz < nz; iz++){
				ss>>source[iz + iy*nny + ix*nnx];
				ss.clear();			}
		}
	}
	ivfile.close();
}
void space::getriemparam(double* param, const string&  towards, int ix, int iy, int iz){ 
	int nx, ny, nz, nny, nnx, index;
	nx = size[0];	ny = size[1];	nz = size[2];
	nny = nz; nnx = ny*nz;
	index = iz + iy*nny + ix*nnx;
	if(ix >= nx || iy >= ny || iz >= nz){		std::cerr<<"function 'getriemparam' index out of range\n";	}
	if(towards.length()>1){
		switch(towards[1]){
			case 'X':
				if(layers[index - nnx] != 0){
					//obstacle
					param[3] = density[index];
					param[4] = xvelosity[index];								
					param[5] = pressure[index];
					param[0] = density[index];
					param[1] = -xvelosity[index];
					param[2] = pressure[index];
				}
				else{
					//gas	
					param[3] = density[index];
					param[4] = xvelosity[index];
					param[5] = pressure[index];
					param[0] = density[index - nnx];
					param[1] = xvelosity[index - nnx];
					param[2] = pressure[index - nnx];
				}			
				break;
			case 'Y':
				if(layers[index - nny] != 0){
					//obstacle
					param[3] = density[index];
					param[4] = yvelosity[index];			
					param[5] = pressure[index];
					param[0] = density[index];
					param[1] = -yvelosity[index];
					param[2] = pressure[index];
				}
				else{
					//gas	
					param[3] = density[index];
					param[4] = yvelosity[index];
					param[5] = pressure[index];
					param[0] = density[index - nny];
					param[1] = yvelosity[index - nny];
					param[2] = pressure[index - nny];
				}	
				break;
			case 'Z':
				if(layers[index - 1] != 0){
					//obstacle
					param[3] = density[index];
					param[4] = zvelosity[index];
					param[5] = pressure[index];
					param[0] = density[index];
					param[1] = -zvelosity[index];
					param[2] = pressure[index];
				}
				else{
					//gas	
					param[3] = density[index];
					param[4] = zvelosity[index];
					param[5] = pressure[index];
					param[0] = density[index - 1];
					param[1] = zvelosity[index - 1];
					param[2] = pressure[index - 1];
				}		
				break;
			default:
				param[0] = density[index];
				param[1] = yvelosity[index];
				param[2] = pressure[index];
				param[3] = density[index];
				param[4] = yvelosity[index];
				param[5] = pressure[index];
		}
	}
	else{
		switch(towards[0]){
			case 'X':
				if(layers[index + nnx] != 0){
					//obstacle
					
					param[0] = density[index];
					param[1] = xvelosity[index];			
					param[2] = pressure[index];
					param[3] = density[index];
					param[4] = -xvelosity[index];
					param[5] = pressure[index];
				}
				else{
					//gas	
					param[0] = density[index];
					param[1] = xvelosity[index];
					param[2] = pressure[index];
					param[3] = density[index + nnx];
					param[4] = xvelosity[index + nnx];
					param[5] = pressure[index + nnx];
				}			
				break;
			case 'Y':
				if(layers[index + nny] != 0){
					//obstacle
					param[0] = density[index];
					param[1] = yvelosity[index];			
					param[2] = pressure[index];
					param[3] = density[index];
					param[4] = -yvelosity[index];
					param[5] = pressure[index];
				}
				else{
					//gas	
					param[0] = density[index];
					param[1] = yvelosity[index];
					param[2] = pressure[index];
					param[3] = density[index + nny];
					param[4] = yvelosity[index + nny];
					param[5] = pressure[index + nny];
				}	
				break;
			case 'Z':
				if(layers[index + 1] != 0){
					//obstacle
					param[0] = density[index];
					param[1] = zvelosity[index];
					param[2] = pressure[index];
					param[3] = density[index];
					param[4] = -zvelosity[index];
					param[5] = pressure[index];
				}
				else{
					//gas	
					param[0] = density[index];
					param[1] = zvelosity[index];
					param[2] = pressure[index];
					param[3] = density[index + 1];
					param[4] = zvelosity[index + 1];
					param[5] = pressure[index + 1];
				}		
				break;
		}		
	}
	return;
}
void space::getcellsize(double* lin, int ix, int iy, int iz){
	/*
	ix, iy, iz from 1 to nx-2, ny-2, nz-2	| n-2
	grid: from 0 to nx-2, ny-2, nz-2		| n-1
	i = between  i and i+1
	*/
	if(ix < 1 || ix > size[0]-2 || iy < 1 || iy > size[1]-2 || iz < 1 || iz > size[2]-2){		std::cerr<<"function 'getcellsize' index out of range\n";	}
 	lin[0] = xpoints[ix] - xpoints[ix-1];
	lin[1] = ypoints[iy] - ypoints[iy-1];
	lin[2] = zpoints[iz] - zpoints[iz-1];
	lin[3] = lin[0]*lin[1]*lin[2];
	return;
}
void space::tanvel(double* Ures, double s2, int tan1, int tan2, const string&  towards, int ix, int iy, int iz){
	int nx, ny, nz, nny, nnx, index;
	double *U = new double[6];
	nx = size[0];	ny = size[1];	nz = size[2];
	nny = nz; nnx = ny*nz;
	index = iz + iy*nny + ix*nnx;
	if(ix >= nx || iy >= ny || iz >= nz){		std::cerr<<"function 'tanvel' index out of range\n";	}
	U[0] = xvelosity[index];	U[1] = yvelosity[index];	U[2] = zvelosity[index];
	if(towards.length()>1){		//backward
		if(s2 <= 0){
				Ures[0] = U[tan1];
				Ures[1] = U[tan2];}
		else{
			switch(towards[1]){//minus direction			
				case 'X':
					U[3] = xvelosity[index-nnx];	U[4] = yvelosity[index -nnx];	U[5] = zvelosity[index -nnx];
					break;
				case 'Y':
					U[3] = xvelosity[index-nny];	U[4] = yvelosity[index -nny];	U[5] = zvelosity[index -nny];
					break;
				case 'Z':
					U[3] = xvelosity[index-1];	U[4] = yvelosity[index -1];	U[5] = zvelosity[index -1];
					break;
				default://case end of empty tube, towards = "Yend"
					U[3] = xvelosity[index];	U[4] = yvelosity[index];	U[5] = zvelosity[index];
			}
			Ures[0] = U[3+tan1];
			Ures[1] = U[3+tan2];
		}
	}
	else{		//forward
		if(s2 >= 0){		Ures[0] = U[tan1]; Ures[1] = U[tan2];}
		else{
			switch(towards[0]){
				case 'X':
					U[3] = xvelosity[index+nnx];	U[4] = yvelosity[index + nnx];	U[5] = zvelosity[index + nnx];
					break;
				case 'Y':
					U[3] = xvelosity[index + nny];	U[4] = yvelosity[index + nny];	U[5] = zvelosity[index + nny];
					break;
				case 'Z':
					U[3] = xvelosity[index + 1];	U[4] = yvelosity[index + 1];	U[5] = zvelosity[index + 1];
					break;
				
			}
			Ures[0] = U[3+tan1]; Ures[1] = U[3+tan2];
		}
	}
	delete[]U;
}
void space::updatetimestep(const string& axis, int i, double s1, double s3, double& ddeltat){
	/*
	ix, iy, iz from 1 to nx-2, ny-2, nz-2	| n-2
	grid: from 0 to nx-2, ny-2, nz-2		| n-1
	i = between  i and i+1
	*/
	int nx, ny, nz;
	nx = size[0];	ny = size[1];	nz = size[2];
	switch(axis[0]){
		case 'Z':
			if(i < 1 || i > nz-2){				std::cerr<<"function 'updatetimestep' index out of range\n";	}
			if(i > 1){				ddeltat = min(ddeltat, fabs((zpoints[i]-zpoints[i-1])/s1));			}
			if(i < nz-2){				ddeltat = min(ddeltat, fabs((zpoints[i+1] - zpoints[i])/s3));			}
			break;
		case 'Y':
			if(i < 1 || i > ny -2){				std::cerr<<"index out of range function getcellsize\n";			}			
			if(i > 1){				ddeltat = min(ddeltat, fabs((ypoints[i]-ypoints[i-1])/s1));			}
			if(i < ny-2){				ddeltat = min(ddeltat, fabs((ypoints[i+1] - ypoints[i])/s3));				}
			break;
		case 'X':
			if(i > 1){				ddeltat = min(ddeltat, fabs((xpoints[i]-xpoints[i-1])/s1));			}
			if(i < nx-2){				ddeltat = min(ddeltat, fabs((xpoints[i+1] - xpoints[i])/s3));				}			
			break;
	}
	if(ddeltat < 1e-2 || ddeltat > 2){		std::cerr<<axis<<"function 'updatetimestep' time step error deltat = "<<ddeltat<<'\n';}
	return;	
}
void space::stopheat(){
	int n = size[0]*size[1]*size[2];
	for(int i = 0; i < n; i++){ source[i] = 0;}
}
void space::showall(ofstream& logfile) const{
	int nx, ny, nz, nny, nnx;
	nx = size[0];	ny = size[1];	nz = size[2];
	nny = nz;	nnx = nz*ny;
	logfile<<"xpoints\n";
	for(int ix = 0; ix < nx-1; ix++){		logfile<<xpoints[ix]<<'\t';	}
	logfile<<"\nypoints\n";
	for(int iy = 0; iy < ny-1; iy++){		logfile<<ypoints[iy]<<'\t';	}
	logfile<<"\nzpoints\n";
	for(int iz = 0; iz < nz-1; iz++){		logfile<<zpoints[iz]<<'\t';	}
	logfile<<"density\n";
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				logfile<<density[iz+iy*nny+ix*nnx]<<'\t';			}
			logfile<<'\n';
		}
		logfile<<'\n';
	}
	
	logfile<<"x velosity\n";
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				logfile<<xvelosity[iz+iy*nny+ix*nnx]<<'\t';			}
			logfile<<'\n';
		}
		logfile<<'\n';
	}
	logfile<<"y velosity\n";
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				logfile<<yvelosity[iz+iy*nny+ix*nnx]<<'\t';			}
			logfile<<'\n';
		}
		logfile<<'\n';
	}
	logfile<<"z velosity\n";
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				logfile<<zvelosity[iz+iy*nny+ix*nnx]<<'\t';			}
			logfile<<'\n';
		}
		logfile<<'\n';
	}
	logfile<<"pressure\n";
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				logfile<<pressure[iz+iy*nny+ix*nnx]<<'\t';			}
			logfile<<'\n';
		}
		logfile<<'\n';
	}
}
void space::savedata(string& datafname){
	ofstream datafile(datafname.c_str());
	int nx, ny, nz, nny, nnx;
	nx = size[0];	ny = size[1];	nz = size[2];
	nny = nz;	nnx = nz*ny;
	datafile<<"VARIABLES = \"X\", \"Y\", \"Z\", \"R\", \"U\", \"V\",\"W\", \"P\", \"L\", \"Q\"\n";
	datafile<<"ZONE F=BLOCK I= "<<nx-1<<" J= "<<ny-1<<" K = "<<nz-1<<'\n';
	datafile<<"VARLOCATION=([4-10]=CELLCENTERED)\n";
	for(int iz = 1; iz < nz; iz++){	
		for(int iy = 1; iy < ny; iy++){
			for(int ix = 1; ix < nx; ix++){		
				datafile<<xpoints[ix-1]<<'\t';
				
			}
			datafile<<'\n';
		}
	}
	for(int iz = 1; iz < nz; iz++){
	
		for(int iy = 1; iy < ny; iy++){
				for(int ix = 1; ix < nx; ix++){		
				datafile<<ypoints[iy-1]<<'\t';
				
			}
			datafile<<'\n';
		}
	}
	for(int iz = 1; iz < nz; iz++){
	
		for(int iy = 1; iy < ny; iy++){
				for(int ix = 1; ix < nx; ix++){		
				datafile<<zpoints[iz-1]<<'\t';
				
			}
			datafile<<'\n';
		}
	}
	for(int iz = 1; iz < nz-1; iz++){
	
		for(int iy = 1; iy < ny-1; iy++){
				for(int ix = 1; ix < nx-1; ix++){
				datafile<<density[iz+iy*nny+ix*nnx]<<'\t';			}
			datafile<<'\n';
		}
	}
	for(int iz = 1; iz < nz-1; iz++){
	
		for(int iy = 1; iy < ny-1; iy++){
				for(int ix = 1; ix < nx-1; ix++){
				datafile<<xvelosity[iz+iy*nny+ix*nnx]<<'\t';			}
			datafile<<'\n';
		}
	}
	for(int iz = 1; iz < nz-1; iz++){
	
		for(int iy = 1; iy < ny-1; iy++){
				for(int ix = 1; ix < nx-1; ix++){
				datafile<<yvelosity[iz+iy*nny+ix*nnx]<<'\t';			}
			datafile<<'\n';
		}
	}
	for(int iz = 1; iz < nz-1; iz++){
	
		for(int iy = 1; iy < ny-1; iy++){
				for(int ix = 1; ix < nx-1; ix++){
				datafile<<zvelosity[iz+iy*nny+ix*nnx]<<'\t';			}
			datafile<<'\n';
		}
	}
	for(int iz = 1; iz < nz-1; iz++){
	
		for(int iy = 1; iy < ny-1; iy++){
				for(int ix = 1; ix < nx-1; ix++){
				datafile<<pressure[iz+iy*nny+ix*nnx]<<'\t';			}
			datafile<<'\n';
		}
	}
	for(int iz = 1; iz < nz-1; iz++){
	
		for(int iy = 1; iy < ny-1; iy++){
				for(int ix = 1; ix < nx-1; ix++){
				datafile<<layers[iz+iy*nny+ix*nnx]<<'\t';			}
			datafile<<'\n';
		}
	}
	for(int iz = 1; iz < nz-1; iz++){
	
		for(int iy = 1; iy < ny-1; iy++){
				for(int ix = 1; ix < nx-1; ix++){
				datafile<<source[iz+iy*nny+ix*nnx]<<'\t';			}
			datafile<<'\n';
		}
	}
	datafile.close();
}
void space::savestate(string& statefname, double& t, double& tmax, int& nstep, int& nstepmax){
	ofstream statef(statefname.c_str());
	statef<<t;			statef<<'\n';	
	statef<<tmax;		statef<<'\n';	
	statef<<nstep;		statef<<'\n';
	statef<<nstepmax;	statef<<'\n';
	statef<<t;			statef<<'\n';
	int nx, ny, nz, nny, nnx;
	nx = size[0];	ny = size[1];	nz = size[2];
	nny = nz;	nnx = nz*ny;
	//for(int ix = 0; ix < nx-1; ix++){	} grid from gridfile
	statef<<"layers\n";
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				statef<<layers[iz + iy*nny + ix*nnx];
				statef<<'\t';			}
			statef<<'\n';
		}
	}
	statef<<"density\n";
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				statef<<density[iz + iy*nny + ix*nnx];
				statef<<'\t';			}
			statef<<'\n';
		}
	}
	statef<<"x-vel\n";
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				statef<<xvelosity[iz + iy*nny + ix*nnx];
				statef<<'\t';			}
			statef<<'\n';
		}
	}
	statef<<"y-vel\n";
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				statef<<yvelosity[iz + iy*nny + ix*nnx];
				statef<<'\t';			}
			statef<<'\n';
		}
	}
	statef<<"z-vel\n";
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				statef<<zvelosity[iz + iy*nny + ix*nnx];
				statef<<'\t';			}
			statef<<'\n';
		}
	}
	statef<<"pressure\n";
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				statef<<pressure[iz + iy*nny + ix*nnx];
				statef<<'\t';			}
			statef<<'\n';
		}
	}
	statef<<"source\n";
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				statef<<source[iz + iy*nny + ix*nnx];
				statef<<'\t';			}
			statef<<'\n';
		}
	}
	statef.close();
}
void space::restore(string& statefname, double& t, int& nstep){
	double tmax, nstepmax;
	ifstream statef(statefname.c_str());
	string s;
	getline(statef, s);	stringstream ss(s);		ss>>t;	ss.clear();
	getline(statef, s);	ss.str(s);				ss>>tmax;	ss.clear();
	getline(statef, s);	ss.str(s);				ss>>nstep;	ss.clear();
	getline(statef, s);	ss.str(s);				ss>>nstepmax;	ss.clear();
	getline(statef, s);	ss.str(s);				ss>>deltat;	ss.clear();
	int nx, ny, nz, nny, nnx;
	nx = size[0];	ny = size[1];	nz = size[2];
	nny = nz;	nnx = nz*ny;
	getline(statef, s);
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			getline(statef, s);
			ss.str(s);
			for(int iz = 0; iz < nz; iz++){
				ss>>layers[iz + iy*nny + ix*nnx];
				ss.clear();	}
		}
	}
	getline(statef, s);
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			getline(statef, s);
			ss.str(s);
			for(int iz = 0; iz < nz; iz++){
				ss>>density[iz + iy*nny + ix*nnx];
				ss.clear();}
		}
	}
	getline(statef, s);
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			getline(statef, s);
			ss.str(s);
			for(int iz = 0; iz < nz; iz++){
				ss>>xvelosity[iz + iy*nny + ix*nnx];
				ss.clear();			}
		}
	}
	getline(statef, s);
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			getline(statef, s);
			ss.str(s);
			for(int iz = 0; iz < nz; iz++){
				ss>>yvelosity[iz + iy*nny + ix*nnx];
				ss.clear();			}
		}
	}
	getline(statef, s);
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			getline(statef, s);
			ss.str(s);
			for(int iz = 0; iz < nz; iz++){
				ss>>zvelosity[iz + iy*nny + ix*nnx];
				ss.clear();			}
		}
	}
	getline(statef, s);
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			getline(statef, s);
			ss.str(s);
			for(int iz = 0; iz < nz; iz++){
				ss>>pressure[iz + iy*nny + ix*nnx];
				ss.clear();			}
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			getline(statef, s);
			ss.str(s);
			for(int iz = 0; iz < nz; iz++){
				ss>>source[iz + iy*nny + ix*nnx];
				ss.clear();			}
		}
	}
	statef.close();
}