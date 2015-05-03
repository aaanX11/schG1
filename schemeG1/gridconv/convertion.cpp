#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>

void findfilenames(std::string &gridfname, std::string &coordfname, std::string &valfname){
	std::cout<<"Input grid file name:\n";
	std::cin>>gridfname;
	coordfname.assign("grid.lst");
	valfname.assign("grid.res");
	return;
}

void grid(const std::string& gridfile, int* size, double** x, double** y, double** z){
	int n;
	double *coordx, *coordy, *coordz;
	std::string s;
	std::ifstream gridf(gridfile.c_str());
	std::stringstream ss;
	if(!gridf){std::cout<<"Cannot open "<<gridfile<<".\n"; exit(1);}
	//X: size[0], x[]
	while(getline(gridf, s) != NULL && s[0] != 'X'){	}
	getline(gridf, s);	ss.str(s);	ss>>n; ss.clear();
	size[0] = n+1;		
	coordx = new double[n+1];
	getline(gridf, s);	ss.str(s);	for(int i = 0; i < size[0]; i++){ss>>coordx[i];	}	ss.clear();
	*x = coordx;
	//Y: size[1], y[]
	while(getline(gridf, s) != NULL && s[0] != 'Y'){	}
	getline(gridf, s);	ss.str(s);	ss>>n;	ss.clear();
	size[1] = n+1;
	coordy = new double[n+1];
	getline(gridf, s); 	ss.str(s);	for(int i = 0; i < size[1]; i++){ss>>coordy[i];	}	ss.clear();
	*y = coordy;
	//Z: size[2], z[]
	while(getline(gridf, s) != NULL && s[0] != 'Z'){	}
	getline(gridf, s);	ss.str(s);	ss>>n;	ss.clear();
	size[2] = n+1;
	coordz = new double[n+1];
	getline(gridf, s);	ss.str(s);	for(int i = 0; i < size[3]; i++){ss>>coordz[i];	}	ss.clear();
	*z = coordz;
	gridf.close();
}
/*
Binary search in the sorted array
Searches two given points:
x1 < x2
Answer: 
i1 <= i2 
ordinal numbers of grid points which are closest to the left of given point
*/
void binsearch(int size, double* x, double x1, double x2, int& i1, int& i2){
	int l;
	i1 = l = size/2;
	while(! (x[i1] <= x1 && x[i1+1] > x1)){
		l = l/2;
		if(x[i1] > x1){i1 -= l;}
		else{i1 += l;}
	}
	i2 = l = (size - i1)/2;
	while(! (x[i1 + i2] <= x2 && x[i1+i2+1] > x2)){
		l = l/2;
		if(x[i1 + i2] > x2){i2 -= l;}
		else{i2 += l;}
	}
	i2 = i1 + i2;
}
/*
searches ordinal numbers
6 given points: 
axis X:	factcell[0], factcell[3]
axis Y:	factcell[1], factcell[4]
axis Z:   factcell[2], factcell[5]

6 sought-for numbers: 
grid X:	coord[0], coord[3]
grid Y:	coord[1], coord[4]
grid Z:	coord[2], coord[5]

size[0] - size of grid x
size[1] - size of grid y
size[2] - size of grid z
*/
void findord(double* factcell, int* size, double* x, double* y, double* z, int* coord){
	int i1, i2;
	//X
	binsearch(size[0], x, factcell[0], factcell[3], i1, i2);
	coord[0] = i1; coord[3] = i2;
	//Y
	binsearch(size[1], y, factcell[1], factcell[4], i1, i2);
	coord[1] = i1; coord[4] = i2;
	//Z
	binsearch(size[2], z, factcell[2], factcell[5], i1, i2);
	coord[2] = i1; coord[5] = i2;
}

void fill(double* factcell, int* size, double* x, double* y, double* z, int* coord, double* values, double v){
	int dix, diy, diz;
	double *xv, *yv, *zv;
	dix = coord[3] - coord[0] + 1;	diy = coord[4] - coord[1] + 1;	diz = coord[5] - coord[2] + 1;
	xv = new double[dix+1];			yv = new double[diy+1];			zv = new double[diz+1];
	for(int i = 0; i < dix*diy*diz; i++){values[i] = v;}
	//X
	int ixv = 0;
	xv[ixv] = factcell[0];
	ixv++;
	for(int ix = coord[0] + 1; ix <= coord[3]; ix++){xv[ixv] = x[ix]; ixv++;}
	xv[ixv] = factcell[3];
	//Y
	int iyv = 0;
	yv[iyv] = factcell[1];
	iyv++;
	for(int iy = coord[1] + 1; iy <= coord[4]; iy++){yv[iyv] = y[iy]; iyv++;}
	yv[iyv] = factcell[4];
	//Z
	int izv = 0;
	zv[izv] = factcell[2];
	izv++;
	for(int iz = coord[2] + 1; iz <= coord[5]; iz++){zv[izv] = z[iz]; izv++;}
	zv[izv] = factcell[5];
	//fill
	for(izv = 0; izv < diz; izv++){
		for(iyv = 0; iyv < diy; iyv++){
			for(izv = 0; izv < diz; izv++){
				values[izv + iyv*diz + ixv*diz*diy] *= (xv[ixv+1]-xv[ixv])*(yv[iyv+1]-xv[iyv])*(zv[izv+1]-zv[izv]);
			}
		}
	}
}

void main(){
	std::string newgridfile, coordfile, valfile, s1, s2;
	int nx, ny, nz, nnx, nny, size[3], newsize[3], newcoord[6];
	double 		*newx, *newy, *newz,		//new grid for property distribution
		*edistr,
		factcell[6], v;
	newx = new double[];	newy = new double[];	newz = new double[];
	findfilenames(newgridfile, coordfile, valfile);
	grid(newgridfile, newsize, &newx, &newy, &newz);
	nny = newsize[2]; nnx = newsize[2]*newsize[1];
	edistr = new double[newsize[0]*newsize[1]*newsize[2]];
	for(int i = 0; i < newsize[0]*newsize[1]*newsize[2]; i++){edistr[i] = 0.0;}
	std::ifstream coord(coordfile.c_str());	
	std::ifstream val(valfile.c_str());
	if(!coord || !val){std::cout<<"Cannot open files 'grid.lst' and 'grid.tot'.\n"; return;}
	while(getline(val, s2) != NULL){
		std::stringstream ss(s2);
		while(ss && ss>>v){	
			if(getline(coord, s1) != NULL){
				std::stringstream ss1(s1);
				ss1>>factcell[0]>>factcell[1]>>factcell[2]>>factcell[3]>>factcell[4]>>factcell[5];	
				findord(factcell, newsize, newx, newy, newz, newcoord);
				int dix = (newcoord[3]-newcoord[0] + 1);
				int diy = (newcoord[4]-newcoord[1] + 1);
				int diz = (newcoord[5]-newcoord[2] + 1);
				double* values = new double[dix*diy*diz];
				fill(factcell, newsize, newx, newy, newz, newcoord, values, v);
				for(int ix = newcoord[0]; ix <= newcoord[3]; ix++){
					for(int iy = newcoord[1]; iy <= newcoord[4]; iy++){
						for(int iz = newcoord[2]; iz <= newcoord[5]; iz++){
							int ix_ = ix - newcoord[0];
							int iy_ = iy - newcoord[1];
							int iz_ = iz - newcoord[2];
							edistr[iz + iy*nny + ix*nnx] += values[iz_ + iy_*diz + ix_*diz*diy];
						}	
					}	
				}
				delete[] values;
			}
			else{std::cout<<"File 'grid.lst' doesn't match 'grid.tot': it is shorter.\n"; return;}
		}
	}
	coord.close();
	val.close();
}