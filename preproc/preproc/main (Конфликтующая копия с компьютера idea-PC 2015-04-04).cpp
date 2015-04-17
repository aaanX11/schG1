#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <crtdbg.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <algorithm>

void readgrid(double* coord, const std::string& gridfname,const int& size,const std::string& axis){
	double tmp;
	std::string s;
	std::stringstream ss;
	std::ifstream grid;
	grid.open(gridfname.c_str(),std::ifstream::in);
	if(!grid){
		std::cerr<<"function 'readgrid' cannot open file "<<gridfname<<'\n';
		exit(1);
	}
	while(getline(grid, s) != NULL && s[0] != axis[0]){	}
	getline(grid, s);
	getline(grid, s);
	ss.str(s);
	for(int i = 0; i < size + 1; i++){
		ss>>tmp;
		coord[i] = tmp;
	}	
	ss.clear();
	grid.close();
	return;
}

void readgridaxis(const std::string& gridfname, int* size){
	int n;
	std::string s;
	std::stringstream ss;
	std::ifstream grid(gridfname.c_str());
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

void findfilenames(std::string &gridfname, std::string &cellfname, std::string &initvalfname, std::string &schparamfname, std::string &logfname, std::string &datafname, std::string& ivfilename){
	std::ifstream perenos("perenos");
	if(!perenos){
		std::cerr<<"function 'findfilenames' cannot open file perenos\n";
		return;
	}
	int i = 0;
	std::string s;
	while(getline(perenos, s) != NULL){
		i++;
		if(i == 5){
			gridfname.assign(s);
		}
		if(i == 6){
			cellfname.assign(s);
		}
	}
	std::ifstream grid(gridfname.c_str());
	if(!grid){
		std::cout<<"function 'findfilenames' grid fails file "<<gridfname<<"\n";
		exit(1);
	}
	grid.close();
	perenos.close();
	initvalfname.assign("initval");
	schparamfname.assign("schparam");
	logfname.assign("log");
	datafname.assign("data");
	ivfilename.assign("state");
	return;
}

void f1(int* size, int* layers, double* density, double* xvelosity, double* yvelosity, double* zvelosity, double* pressure){
	//empty closed box, still air
	int nx, ny, nz, nnx, nny;
	nx = size[0];ny = size[1];nz = size[2]; nny = nz; nnx = nz*ny;
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				layers[iz + iy*nny + ix*nnx] = 0;
				density[iz + iy*nny + ix*nnx] = 1.0;
				xvelosity[iz + iy*nny + ix*nnx] = 0;
				yvelosity[iz + iy*nny + ix*nnx] = 0;
				zvelosity[iz + iy*nny + ix*nnx] = 0;
				pressure[iz + iy*nny + ix*nnx] = 1.0;
			}
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			{int iz = 0;
				//z1-wall 
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			{int iz = nz-1;
				//z2-wall 
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
	for(int ix = 0; ix < nx; ix++){
		{int iy = 0; 
			for(int iz = 0; iz < nz; iz++){
				//y1-wall 
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
	for(int ix = 0; ix < nx; ix++){
		{int iy = ny-1; 
			for(int iz = 0; iz < nz; iz++){
				//y2-wall 
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
	{int ix = 0;
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				//x1-wall 
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
	{int ix = nx-1;
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				//x2-wall 
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
}

void f2(int* size, int* layers, double* density, double* xvelosity, double* yvelosity, double* zvelosity, double* pressure){
	//empty box, front wall missing, Y-flux
	double vel;
	std::cout<<"velosity: v = ?\n";
	std::cin>>vel;
	int nx, ny, nz, nnx, nny;
	nx = size[0]; ny = size[1];nz = size[2]; nny = nz; nnx = nz*ny;
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				layers[iz + iy*nny + ix*nnx] = 0;
				density[iz + iy*nny + ix*nnx] = 1.0;
				xvelosity[iz + iy*nny + ix*nnx] = 0;
				yvelosity[iz + iy*nny + ix*nnx] = vel;
				zvelosity[iz + iy*nny + ix*nnx] = 0;
				pressure[iz + iy*nny + ix*nnx] = 1.0;
			}
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			{int iz = 0;
				//z1-wall 
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			{int iz = nz-1;
				//z2-wall 
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
	//no y1-wall
	for(int ix = 0; ix < nx; ix++){
		{int iy = ny-1; 
			for(int iz = 0; iz < nz; iz++){
				//y2-wall 
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
	{int ix = 0;
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				//x1-wall 
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
	{int ix = nx-1;
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				//x2-wall 
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
}

void f3(int* size, int* layers, double* density, double* xvelosity, double* yvelosity, double* zvelosity, double* pressure){
	//empty tube
	double vel;
	std::cout<<"velosity: v = ?\n";
	std::cin>>vel;
	int nx, ny, nz, nnx, nny;
	nx = size[0];ny = size[1]; nz= size[2]; nny = nz; nnx = nz*ny;
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				layers[iz + iy*nny + ix*nnx] = 0;
				density[iz + iy*nny + ix*nnx] = 1.0;
				xvelosity[iz + iy*nny + ix*nnx] = 0;
				yvelosity[iz + iy*nny + ix*nnx] = vel;
				zvelosity[iz + iy*nny + ix*nnx] = 0;
				pressure[iz + iy*nny + ix*nnx] = 1.0;
			}
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			{int iz = 0;
				//z1-wall 
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			{int iz = nz-1;
				//z2-wall 
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
	{int ix = 0;
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				//x1-wall 
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
	{int ix = nx-1;
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				//x2-wall 
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
}

void squareobstacle(int* size, int* layers){
	int nx, ny, nz, nnx, nny;
	nx = size[0]; ny = size[1]; nz = size[2]; nny = nz; nnx = nz*ny;
	int x1, x2, y1, y2, z1, z2;
	x1 = (nx-1)/2 - (nx-1)/20;
	x2 = (nx-1)/2 + (nx-1)/20;
	y1 = (ny-1)/2 - (ny-1)/20;
	y2 = (ny-1)/2 + (ny-1)/20;
	z1 = (nz-1)/2 - (nz-1)/20;
	z2 = (nz-1)/2 + (nz-1)/20;
	
	for(int ix = 1; ix < nx-1; ix++){		
		for(int iy = 1; iy < ny-1; iy++){
			for(int iz = 1; iz < nz-1; iz++){//19 20 21 22 | 37 38 39 40 41 42 43 44
				if(x1 < ix && ix <= x2 && y1 < iy && iy <= y2  && z1 < iz && iz <= z2 ){
					layers[iz+iy*nny+ix*nnx]=-1;}
			}
		}
	}
}

void squareobstacle(int* size,int* layers, int x1, int x2, int y1, int y2, int z1, int z2){
	int nx, ny, nz, nnx, nny;
	nx = size[0];ny = size[1];nz = size[2]; nny = nz; nnx = nz*ny;
	for(int ix = x1; ix < x2; ix++){
		for(int iy = y1; iy < y2; iy++){
			for(int iz = z1; iz < z2; iz++){
				layers[iz + iy*nny + ix*nnx] = -1;
			}
		}
	}
}


void f4(int* size, int* layers, double* density, double* xvelosity, double* yvelosity, double* zvelosity, double* pressure){
	//closed box, still air, obstacle
	int nx, ny, nz, nnx, nny;
	nx = size[0]; ny = size[1]; nz = size[2]; nny = nz; nnx = nz*ny;
	f1(size, layers, density, xvelosity, yvelosity, zvelosity, pressure);
	squareobstacle(size, layers);
}
void f5(int* size, int* layers, double* density, double* xvelosity, double* yvelosity, double* zvelosity, double* pressure){
	//obstacle in a tube
	f3(size, layers, density, xvelosity, yvelosity, zvelosity, pressure);
	squareobstacle(size, layers);
}

void f6(int* size, int* layers, double* density, double* xvelosity, double* yvelosity, double* zvelosity, double* pressure){
	//obstacle in a tube
	int x1, x2,y1, y2, z1, z2;
	f3(size, layers, density, xvelosity, yvelosity, zvelosity, pressure);
	std::cout<<"specify position: x1 < x2\ny1 < y2\nz1 < z2\n";
	std::cin>>x1>>x2;
	std::cin>>y1>>y2;
	std::cin>>z1>>z2;
	squareobstacle(size, layers, x1, x2,y1, y2, z1,z2);
}

void write(const std::string& initvalfname, int* size, int* layers, double* density, double* xvelosity, double* yvelosity, double* zvelosity, double* pressure){
	std::ofstream ivfile(initvalfname.c_str());
	if(!ivfile){
		std::cerr<<"Can't write to "<<initvalfname<<'\n';
		return;
	}
	std::string s;
	int nx, ny, nz, nny, nnx;
	nx = size[0];
	ny = size[1];
	nz = size[2];
	nny = nz;
	nnx = nz*ny;
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				ivfile<<layers[iz + iy*nny + ix*nnx];
				ivfile<<'\t';
			}
			ivfile<<'\n';
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				ivfile<<density[iz + iy*nny + ix*nnx];
				ivfile<<'\t';
			}
			ivfile<<'\n';
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				ivfile<<xvelosity[iz + iy*nny + ix*nnx];
				ivfile<<'\t';
			}
			ivfile<<'\n';
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				ivfile<<yvelosity[iz + iy*nny + ix*nnx];
				ivfile<<'\t';
			}
			ivfile<<'\n';
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				ivfile<<zvelosity[iz + iy*nny + ix*nnx];
				ivfile<<'\t';
			}
			ivfile<<'\n';
		}
	}
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				ivfile<<pressure[iz + iy*nny + ix*nnx];
				ivfile<<'\t';
			}
			ivfile<<'\n';
		}
	}
	ivfile.close();
}

typedef void (*pf)(int*, int*, double*, double*, double*, double*, double*);

void main(){
	int size[3], *layers, desc;
	double *density, *xvelosity, *yvelosity, *zvelosity, *pressure;
	std::string gridfname, datafname, cellfname, initvalfname, schparamfname, logfname, ivfilename;
	std::cout<<"Choose a model:\n";
	std::cout<<"1 - still air, closed box\n";
	std::cout<<"2 - shock wave in the box without front y-wall\n";
	std::cout<<"3 - y-flux through a tube\n";
	std::cout<<"4 - an obstacle in the middle of a box\n";
	std::cout<<"5 - an obstacle in the middle of a tube with y-flux\n";
	std::cout<<"6 - an obstacle with arbitrary position in a tube\n";
	std::cin>>desc;
	findfilenames(gridfname, cellfname, initvalfname, schparamfname, logfname, datafname, ivfilename);
	readgridaxis(gridfname, size);
	if(size[0] > 1000 || size[0] < 1 || size[1] > 1000 || size[1] < 1 || size[2] > 1000 || size[2] < 1){exit(1);}
	std::cout<<"space size: "<<size[0]<<" "<<size[1]<<" "<<size[2]<<'\n';
	pf f[]= {&f1, &f2, &f3, &f4, &f5, &f6};
	layers = new int[size[0]*size[1]*size[2]];
	density = new double[size[0]*size[1]*size[2]];
	xvelosity = new double[size[0]*size[1]*size[2]];
	yvelosity = new double[size[0]*size[1]*size[2]];
	zvelosity = new double[size[0]*size[1]*size[2]];
	pressure = new double[size[0]*size[1]*size[2]];
	(f[desc-1])(size, layers, density, xvelosity, yvelosity, zvelosity, pressure);
	write(initvalfname, size, layers, density, xvelosity, yvelosity, zvelosity, pressure);
}