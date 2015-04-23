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

void processgrid(const std::string& gridfname, double** x, double** y, double** z, int* size){
	int n;
	double tmp, *coord1, *coord2, *coord3;
	std::string s;
	std::stringstream ss;
	std::ifstream grid(gridfname.c_str());
	if(!grid){
		std::cerr<<"function 'readgridaxis' cannot open file "<<gridfname<<'\n';
		exit(1);
	}
	//----------------------X------------------------
	while(getline(grid, s) != NULL && s[0] != 'X'){	}
	getline(grid, s);	ss.str(s);	ss>>n;size[0] = n+2;	ss.clear();
	coord1 = new double[n+1];
	getline(grid, s);	ss.str(s);
	for(int i = 0; i < n + 1; i++){		ss>>tmp;		coord1[i] = tmp;	}	
	*x = coord1;	ss.clear();
	//----------------------Y------------------------
	while(getline(grid, s) != NULL && s[0] != 'Y'){	}
	getline(grid, s);	ss.str(s);	ss>>n;size[1] = n+2;	ss.clear();
	coord2 = new double[n+1];
	getline(grid, s);	ss.str(s);
	for(int i = 0; i < n + 1; i++){		ss>>tmp;		coord2[i] = tmp;	}	
	*y = coord2;	ss.clear();
	//----------------------Z------------------------
	while(getline(grid, s) != NULL && s[0] != 'Z'){	}
	getline(grid, s);	ss.str(s);	ss>>n;size[2] = n+2;	ss.clear();
	coord3 = new double[n+1];
	getline(grid, s);	ss.str(s);
	for(int i = 0; i < n + 1; i++){		ss>>tmp;		coord3[i] = tmp;	}	
	*z = coord3;	ss.clear();
	return;
}

void findfilenames(std::string &gridfname, std::string &cellfname, std::string &initvalfname, std::string &schparamfname, std::string &logfname, std::string &datafname, std::string& ivfilename, std::string &coord, std::string& val){
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
	coord.assign("grid.lst");
	val.assign("grid.res");
	return;
}

void f1(int* size, int* layers, double* density, double* xvelosity, double* yvelosity, double* zvelosity, double* pressure, double* hsource){
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
				hsource[iz + iy*nny + ix*nnx] = 0.0;
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

void f2(int* size, int* layers, double* density, double* xvelosity, double* yvelosity, double* zvelosity, double* pressure, double* hsource){
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
				hsource[iz + iy*nny + ix*nnx] = 0.0;
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

void f3(int* size, int* layers, double* density, double* xvelosity, double* yvelosity, double* zvelosity, double* pressure, double* hsource){
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
				hsource[iz + iy*nny + ix*nnx] = 0.0;
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


void f4(int* size, int* layers, double* density, double* xvelosity, double* yvelosity, double* zvelosity, double* pressure, double* hsource){
	//closed box, still air, obstacle
	int nx, ny, nz, nnx, nny;
	nx = size[0]; ny = size[1]; nz = size[2]; nny = nz; nnx = nz*ny;
	f1(size, layers, density, xvelosity, yvelosity, zvelosity, pressure, hsource);
	squareobstacle(size, layers);
}
void f5(int* size, int* layers, double* density, double* xvelosity, double* yvelosity, double* zvelosity, double* pressure, double* hsource){
	//obstacle in a tube
	f3(size, layers, density, xvelosity, yvelosity, zvelosity, pressure, hsource);
	squareobstacle(size, layers);
}

void f6(int* size, int* layers, double* density, double* xvelosity, double* yvelosity, double* zvelosity, double* pressure, double* hsource){
	//obstacle in a tube
	int x1, x2,y1, y2, z1, z2;
	f3(size, layers, density, xvelosity, yvelosity, zvelosity, pressure, hsource);
	std::cout<<"specify position: x1 < x2\ny1 < y2\nz1 < z2\n";
	std::cin>>x1>>x2;
	std::cin>>y1>>y2;
	std::cin>>z1>>z2;
	squareobstacle(size, layers, x1, x2,y1, y2, z1,z2);
}

void f7(int* size, int* layers, double* density, double* xvelosity, double* yvelosity, double* zvelosity, double* pressure, double* hsource){
	//processing.CEL file
	return;
}

void heat(int* size, double* hsource){
	int nx, ny, nz, nny, nnx;
	nx = size[0];
	ny = size[1];
	nz = size[2];
	nny = nz;
	nnx = nz*ny;
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < (int)(ny/4); iy++){
			for(int iz = 0; iz < nz; iz++){
				hsource[iz + iy*nny + ix*nnx] = 1.0;
			}
		}
	}
}

void writeinit(const std::string& initvalfname, int* size, int* layers, double* density, double* xvelosity, double* yvelosity, double* zvelosity, double* pressure, double* hsource){
	std::ofstream ivfile("initval");
	std::ofstream datafile("init.dat");
	std::string s;
	int nx, ny, nz, nny, nnx;
	nx = size[0];	ny = size[1];	nz = size[2];
	nny = nz;	nnx = nz*ny;
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
	for(int ix = 0; ix < nx; ix++){
		for(int iy = 0; iy < ny; iy++){
			for(int iz = 0; iz < nz; iz++){
				ivfile<<hsource[iz + iy*nny + ix*nnx];
				ivfile<<'\t';
			}
			ivfile<<'\n';
		}
	}
	ivfile.close();
}
void writedat(const std::string& datafname,const int* size, const int* layers,const double* density,const double* xvelosity,const double* yvelosity,const double* zvelosity,const double* pressure,const double* hsource,const double* xpoints,const double* ypoints,const double* zpoints){
	std::ofstream datafile(datafname.c_str());
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
				datafile<<hsource[iz+iy*nny+ix*nnx]<<'\t';			}
			datafile<<'\n';
		}
	}
	datafile.close();
}
typedef void (*pf)(int*, int*, double*, double*, double*, double*, double*, double*);

void main(){
	int size[3], *layers, desc;
	double *density, *xvelosity, *yvelosity, *zvelosity, *pressure, *hsource, *x, *y , *z;
	std::string gridfname, datafname, cellfname, initvalfname, schparamfname, logfname, ivfilename, scoordfname, svalfname;
	std::cout<<"Choose a model:\n";
	std::cout<<"1 - still air, closed box\n";
	std::cout<<"2 - shock wave in the box without front y-wall\n";
	std::cout<<"3 - y-flux through a tube\n";
	std::cout<<"4 - an obstacle in the middle of a box\n";
	std::cout<<"5 - an obstacle in the middle of a tube with y-flux\n";
	std::cout<<"6 - an obstacle with arbitrary position in a tube\n";
	//7 from cell file
	std::cin>>desc;
	findfilenames(gridfname, cellfname, initvalfname, schparamfname, logfname, datafname, ivfilename, scoordfname, svalfname);
	x = new double[]; y = new double[]; z = new double[];
	processgrid(gridfname, &x, &y, &z, size);
	if(size[0] > 1000 || size[0] < 1 || size[1] > 1000 || size[1] < 1 || size[2] > 1000 || size[2] < 1){exit(1);}
	std::cout<<"space size: "<<size[0]<<" "<<size[1]<<" "<<size[2]<<'\n';
	pf f[]= {&f1, &f2, &f3, &f4, &f5, &f6};
	layers = new int[size[0]*size[1]*size[2]];
	density = new double[size[0]*size[1]*size[2]];
	xvelosity = new double[size[0]*size[1]*size[2]];
	yvelosity = new double[size[0]*size[1]*size[2]];
	zvelosity = new double[size[0]*size[1]*size[2]];
	pressure = new double[size[0]*size[1]*size[2]];
	hsource = new double[size[0]*size[1]*size[2]];
	(f[desc-1])(size, layers, density, xvelosity, yvelosity, zvelosity, pressure, hsource);
	//std::cout<<"looking for heat source in files '"<<scoordfname<<", "<<svalfname<<"'\n";
	heat(size,hsource);
	writeinit("initval", size, layers, density, xvelosity, yvelosity, zvelosity, pressure, hsource);
	writedat("init.dat", size, layers, density, xvelosity, yvelosity, zvelosity, pressure, hsource, x, y, z);

}