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

#ifdef  _CRTDBG_MAP_ALLOC

#define malloc(s)      _malloc_dbg(s,_NORMAL_BLOCK,__FILE__,__LINE__)
#define calloc(c,s)    _calloc_dbg(c,s,_NORMAL_BLOCK,__FILE__,__LINE__)
#define realloc(p,s)   _realloc_dbg(p,s,_NORMAL_BLOCK,__FILE__,__LINE__)
#define _expand(p,s)   _expand_dbg(p,s,_NORMAL_BLOCK,__FILE__,__LINE__)
#define free(p)        _free_dbg(p,_NORMAL_BLOCK)
#define _msize(p)      _msize_dbg(p,_NORMAL_BLOCK)

#endif  /* _CRTDBG_MAP_ALLOC */

#ifdef _DEBUG
#ifdef _CRTDBG_MAP_ALLOC
#define new new(_NORMAL_BLOCK, __FILE__, __LINE__)
#endif /* _CRTDBG_MAP_ALLOC */
#endif /* _DEBUG */

#include <crtdbg.h>

using namespace std;

void showerror(string fname){
	cout<<"error reading file "<<fname<<"\n";
	exit(1);
	return;
}
void showerror(){
	printf("somewhere emerged an error\n");
	exit(1);
	return;
}
struct norma{
	int n1;
	int n2;
	int n3;
	int n;
	int normindex(){
		switch (n){
			case -1:
				return 0;
			case 1:
				return 1;
			case -2:
				return 2;
			case 2:
				return 3;
			case -3:
				return 4;
			case 3:
				return 5;
			default:
				showerror();
				return -1;
		}
	}
	bool coherence(){
		switch(n){
			case -1:
				if (!(n1 == -1 && n2 == 0 && n3 == 0)){
					return false;
				}
				break;
			case 1:
				if (!(n1 == 1 && n2 == 0 && n3 == 0)){
					return false;
				}
				break;
			case -2:
				if (!(n1 == 0 && n2 == -1 && n3 == 0)){
					return false;
				}
				break;
			case 2:
				if (!(n1 == 0 && n2 == 1 && n3 == 0)){
					return false;
				}
				break;
			case -3:
				if (!(n1 == 0 && n2 == 0 && n3 == -1)){
					return false;
				}

				break;
			case 3:
				if (!(n1 == 0 && n2 == 0 && n3 == 1)){
					return false;
				}
				break;
			default:
				return false;
		
		}
		return true;
	}
	norma(int x){
		n = x;
		switch (n){
			case -1:
				n1 = -1; n2 = 0; n3 = 0;
				break;
			case 1:
				n1 = 1; n2 = 0; n3 = 0;
				break;
			case -2:
				n1 = 0; n2 = -1; n3 = 0;
				break;
			case 2:
				n1 = 0; n2 = 1; n3 = 0;
				break;
			case -3:
				n1 = 0; n2 = 0; n3 = -1;
				break;
			case 3:
				n1 = 0; n2 = 0; n3 = 1;
				break;
			default:
				showerror();
		}
	}
	norma(int x1, int x2, int x3){
		n1 = x1;
		n2 = x2;
		n3 = x3;
		if(n1 == 0 && n2 == 0){
			if(n3 == 1){
				n = 3;
			}
			else{
				if(n3 == -1){
					n = -3;
				}
				else{
					showerror();
				}
			}
		}
		if(n1 == 0 && n3 == 0){
			if(n2 == 1){
				n = 2;
			}
			else{
				if(n2 == -1){
					n = -2;
				}
				else{
					showerror();
				}
			}
		}
		if(n2 == 0 && n3 == 0){
			if(n1 == 1){
				n = 1;
			}
			else{
				if(n1 == -1){
					n = -1;
				}
				else{
					showerror();
				}
			}
		}
		if((n1 != 0 && n2 != 0) || (n1 != 0 && n3 != 0) || (n2 != 0 && n3 != 0)){
			showerror();
		}
	}
	norma(string s){
		stringstream ss(s);

		ss>>n;
		ss>>n1>>n2>>n3;
	}
};


struct coord{
	int x; 
	int y;
	int z;
	void set(coord s){
		x = s.x;
		y = s.y;
		z = s.z;
		return;
	}
	void set(int q){
		x = q;
		y = q;
		z = q;
		return;
	}
	void set(int a, int b, int c){
		x = a;
		y = b;
		z = c;
		return;
	}
};


void toUpper(string& buf){
	//buf.assign(buf.ToUpper());
	return;
	/*
	int i, n;
	char buf1[32];
	int i = 0;
	while(buf[i] != '\0'){
		if(buf[i] <= 'z' && buf[i] >= 'a'){
			buf1[i] = buf[i] - 'a' + 'A';
		}
		else{
			buf1[i] = buf[i];
		}
		i++;
	}
	buf1[i] = '\0';
	return buf1;
	*/
}


bool isinvertion(int norm_prev, int norm){
	bool res;
	res = false;
	if(norm != norm_prev){
		if((norm == -1 && norm_prev == 1) || ((abs(norm) == 1)&& (1 < abs(norm_prev)))){
			res = true;
		}
		if((norm == -2 && norm_prev == 2) || (abs(norm) == 2 && abs(norm_prev) == 3)){
			res = true;
		}
		if(norm = -3 && norm_prev == 3){
			res = true;
		}
	}
	return res;
}

int createnorm(int n1, int n2, int n3){
	if(n1 == 0 && n2 == 0){
		if(n3 == 1){
			return 3;
		}
		else{
			if(n3 == - 1){
				return -3;
			}
		}
	}
	if(n1 == 0 && n3 == 0){
		if(n2 == 1){
			return 2;
		}
		else{
			if(n2 == - 1){
				return -2;
			}
		}
	}
	if(n2 == 0 && n3 == 0){
		if(n1 == 1){
			return 1;
		}
		else{
			if(n1 == - 1){
				return -1;
			}
		}
	}
	return 4;
}
bool ischanged(int norm_prev, int norm){
	if(norm == norm_prev){
		return false;
	}
	else{
		return true;
	}
}





void fillarray(int* a, int n){
	for(int i = 0; i < n; i++){
		a[i] = -1;
	}
}
void fillarray(double* a, int n){
	for(int i = 0; i < n; i++){
		a[i] = 0.0;
	}
}
coord getsizefromgrid(string gridfname){
	int x, y, z;
	ifstream grid(gridfname.c_str());
	if(grid.is_open()){
		string s;
		int i = 0;
		while(getline(grid, s) && i < 13){
			
			i++;
			if(i == 6){
				stringstream ss;
				ss.str(s);
				ss>>x;
			}
			if(i == 9){
				stringstream ss;
				ss.str(s);
				ss>>y;
			}
			if(i == 12){
				stringstream ss;
				ss.str(s);
				ss>>z;
			}
		}
		x = x + 2;
		y = y + 2;
		z = z + 2;
		grid.clear();
		coord q;
		q.set(x, y, z);
		return q;
	}
	else{
		showerror(gridfname);
		coord q;
		q.set(-1);
		return q;
	}
}


vector <int> readbadlayers(string optionfname){
	char buf[32];
	int layer;
	vector <int> badlayers;
	FILE* opt;
	opt = fopen(optionfname.c_str(), "r");
	if(opt){
		while(fgets(buf, 256, opt) != NULL){
			sscanf(buf, "%d", &layer);		
			if(find(badlayers.begin(), badlayers.end(), layer) == badlayers.end()){
				badlayers.push_back(layer);
			}
		}
		fclose(opt);
	}
	return badlayers;
}





class EC{
	int coaxis;
	coord size;
	double* valonplanes;
	double* axis1;
	double* axis2;
public:
	void setEC(int x, int y, int z, int n, double val);
	coord getsize() const{
		coord q;
		q.set(size.x, size.y, size.z);
		return q;
	}
	int getaxis() const{
		return coaxis;
	}
	void isnearborders(int x, int y, int z, int& side, int& corner) const;
	void nearside(int n, double* I) const;
	void nearcorner(int n, double* I) const;
	EC(int sizex, int sizey, int sizez, int coaxis);
	void badlayschangeI(int* laysdisposition, double* I) const;
	void getI(int x, int y, int z, double* I) const;
	void getdelta(int x, int y, int z, double* d) const;
	void show() const;
	void setgrid(const string& s1, const string& s2);
};

EC::EC(int sizex, int sizey, int sizez, int axis){
	size.set(sizex, sizey, sizez);
	valonplanes = new double[sizex*sizey*sizez];
	fillarray(valonplanes, sizex*sizey*sizez);
	coaxis = axis;
	switch(axis){
		case 'X':	
			axis1 = new double[sizey + 1];
			fillarray(axis1, sizey + 1);
			axis2 = new double[sizez + 1];
			fillarray(axis2, sizez + 1);
			break;
		case 'Y':
			axis1 = new double[sizex + 1];
			fillarray(axis1, sizex + 1);
			axis2 = new double[sizez + 1];
			fillarray(axis2, sizez + 1);
			break;
		case 'Z':
			axis1 = new double[sizex + 1];
			fillarray(axis1, sizex + 1);
			axis2 = new double[sizey + 1];
			fillarray(axis2, sizey + 1);
			break;
	}	
}

void EC::setEC(int x, int y, int z, int n, double val){
	int nny, nnx;
	nny = size.z;
	nnx = size.y*size.z;
	if(n < 0){val = -1.0*val;}
	switch (coaxis){
		case 'X':
			if(x < 1 || x > size.x || y < 1 || y >= size.y - 1 || z < 1 || z >= size.z - 1){
				char buf[32];
				sprintf(buf, "invalid coordinates of x-normal %d %d %d", x, y, z);
				throw invalid_argument(buf);
			}
			valonplanes[z + (y)*nny + (x - 1)*nnx] += val;
			break;
		case 'Y':
			if(y < 1 || y > size.y || x < 1 || x >= size.x - 1 || z < 1 || z >= size.z - 1){
				char buf[32];
				sprintf(buf, "invalid coordinates of y-normal %d %d %d", x, y, z);
				throw invalid_argument(buf);
			}
			valonplanes[z + (y - 1)*nny + (x)*nnx] += val;
			break;
		case 'Z':
			if(z < 1 || z > size.z || x < 1 || x >= size.x - 1 || y < 1 || y >= size.y - 1){
				char buf[32];
				sprintf(buf, "invalid coordinates of z-normal %d %d %d", x, y, z);
				throw invalid_argument(buf);
			}
			valonplanes[z - 1 + (y)*nny + (x)*nnx] += val;
			break;
	}
	return;
}
void EC::setgrid(const string& s1, const string& s2){
	int n1, n2;
	stringstream ss1, ss2;
	
	switch(coaxis){
		case 'X':
			n1 = size.y;
			n2 = size.z;
			break;
		case 'Y':
			n1 = size.x;
			n2 = size.z;
			break;
		case 'Z':
			n1 = size.x;
			n2 = size.y;
			break;
	}
	ss1.str(s1);
	for (int j = 1; j < n1; j++){
		ss1>>axis1[j];
	}	
	axis1[0] = axis1[1] - (axis1[2] - axis1[1]);
	axis1[n1] = axis1[n1 - 1] + (axis1[n1 - 1] - axis1[n1 - 2]);
	ss2.str(s2);
	for (int j = 1; j < n2; j++){
		ss2>>axis2[j];
	}
	axis2[0] = axis2[1] - (axis2[2] - axis2[1]);
	axis2[n2] = axis2[n2 - 1] + (axis2[n2 - 1] - axis2[n2 - 2]);
	return;	
}

void EC::badlayschangeI(int* laysdisposition, double* I) const{
	for(int i = 0; i < 4; i++){
		if(laysdisposition[i] == 1){
			I[i] = 0;
			I[i + 4] = 0;
		}
	}
	return;
}
void EC::getdelta(int x, int y, int z, double* d) const{
	double dax1_1, dax1_2, dax2_1, dax2_2;
	switch (coaxis){
		case 'X':
			dax1_1 = axis1[y] - axis1[y - 1];
			dax1_2 = axis1[y + 1] - axis1[y];
			dax2_1 = axis2[z] - axis2[z - 1];
			dax2_2 = axis2[z + 1] - axis2[z];
			break;
		case 'Y':
			dax1_1 = axis1[x] - axis1[x - 1];
			dax1_2 = axis1[x + 1] - axis1[x];
			dax2_1 = axis2[z] - axis2[z - 1];
			dax2_2 = axis2[z + 1] - axis2[z];
			break;
		case 'Z':
			dax1_1 = axis1[x] - axis1[x - 1];
			dax1_2 = axis1[x + 1] - axis1[x];
			dax2_1 = axis2[y] - axis2[y - 1];
			dax2_2 = axis2[y + 1] - axis2[y];
			break;
	}	
	d[0] = dax1_1;
	d[1] = dax1_2;
	d[2] = dax2_1;
	d[3] = dax2_2;
	return;

}
void EC::isnearborders(int x, int y, int z, int& side, int& corner) const{
	side = -1;
	corner = -1;
	switch(coaxis){
		case 'X':
			if(y == 1){
				if(z == 1){
					corner = 0;
				}
				else{
					if(z == size.z - 1){corner = 3;}
					else{side = 0;}
				}
			}
			else{
				if(y == size.y - 1){
					if(z == 1){
						corner = 1;
					}
					else{
						if(z == size.z - 1){corner = 2;}
						else{side = 2;}
					}
				}
				else{
					if(z == 1){side = 3;}
					else{
						if(z == size.z - 1){side = 1;}
					}
				}
			}
			break;
		case 'Y':
			if(x == 1){
				if(z == 1){
					corner = 0;
				}
				else{
					if(z == size.z - 1){corner = 3;}
					else{side = 0;}
				}
			}
			else{
				if(x == size.x - 1){
					if(z == 1){
						corner = 1;
					}
					else{
						if(z == size.z - 1){corner = 2;}
						else{side = 2;}
					}
				}
				else{
					if(z == 1){side = 3;}
					else{
						if(z == size.z - 1){side = 1;}
					}
				}
			}
			break;
		case 'Z':
			if(x == 1){
				if(y == 1){
					corner = 0;
				}
				else{
					if(y == size.y - 1){corner = 3;}
					else{side = 0;}
				}
			}
			else{
				if(x == size.x - 1){
					if(y == 1){
						corner = 1;
					}
					else{
						if(y == size.y - 1){corner = 2;}
						else{side = 2;}
					}
				}			
				else{
					if(y == 1){side = 3;}
					else{
						if(y == size.y - 1){side = 1;}
					}
				}
			}
			break;
	}
	return;
}
void EC::nearside(int n, double* I) const{
	switch (n){
		case 0:
			I[0] = I[1];I[2] = I[3];
			I[4] = I[5];I[6] = I[7];
			break;
		case 1:
			I[2] = I[0]; I[3] = I[1]; 
			I[6] = I[4]; I[7] = I[5];
			break;
		case 2:
			I[3] = I[2]; I[1] = I[0]; 
			I[7] = I[6]; I[5] = I[4]; 
			break;
		case 3:
			I[0] = I[2]; I[1] = I[3]; 
			I[4] = I[6]; I[5] = I[7];
	}
	return;
}
void EC::nearcorner(int n, double* I)const {
	switch (n){
		case 0:
			for(int i = 0; i < 4; i++){
				I[i] = I[3];
				I[i + 4] = I[7];
			}
			break;
		case 1:
			for(int i = 0; i < 4; i++){
				I[i] = I[2];
				I[i + 4] = I[6];
			}
			break;
		case 2:
			for(int i = 0; i < 4; i++){
				I[i] = I[0];
				I[i + 4] = I[4];
			}
			break;
		case 3:
			for(int i = 0; i < 4; i++){
				I[i] = I[1];
				I[i + 4] = I[5];
			}
			break;
	}
	return;
}
void EC::getI(int x, int y, int z, double* I) const{
	int side, corner;
	int nnx, nny;
	side = -1;
	corner = -1;
	nny = size.z;
	nnx = size.y*size.z;
	switch(coaxis){
		case 'X':
			I[0] = valonplanes[(z - 1) + (y - 1)*nny + (x - 1)*nnx];
			I[4] = valonplanes[(z - 1) + (y - 1)*nny + (x)*nnx];

			I[1] = valonplanes[(z - 1) + (y)*nny + (x - 1)*nnx];
			I[5] = valonplanes[(z - 1) + (y)*nny + (x)*nnx];

			I[2] = valonplanes[(z) + (y - 1)*nny + (x - 1)*nnx];
			I[6] = valonplanes[(z) + (y - 1)*nny + (x)*nnx];

			I[3] = valonplanes[(z) + (y)*nny + (x - 1)*nnx];
			I[7] = valonplanes[(z) + (y)*nny + (x)*nnx];
			break;
		case 'Y':
			I[0] = valonplanes[(z - 1) +(y - 1)*nny + (x - 1)*nnx];
			I[4] = valonplanes[(z - 1) + (y)*nny + (x - 1)*nnx];
			I[1] = valonplanes[(z - 1) +(y - 1)*nny + x*nnx];
			I[5] = valonplanes[(z - 1) +(y)*nny + x*nnx];
			I[2] = valonplanes[z +		(y - 1)*nny + (x - 1)*nnx];
			I[6] = valonplanes[z +		(y)*nny + (x - 1)*nnx];
			I[3] = valonplanes[z +		(y - 1)*nny + x*nnx];
			I[7] = valonplanes[z +		(y)*nny + x*nnx];
			break;
		case 'Z':
			I[0] = valonplanes[(z - 1) +	(y - 1)*nny +	(x - 1)*nnx];
			I[4] = valonplanes[(z) +		(y - 1)*nny +	(x - 1)*nnx];
			I[1] = valonplanes[(z - 1) +	(y - 1)*nny +	x*nnx];
			I[5] = valonplanes[(z) +		(y - 1)*nny +	x*nnx];
			I[2] = valonplanes[(z - 1) +	y*nny +		(x - 1)*nnx];
			I[6] = valonplanes[(z) +		y*nny +		(x - 1)*nnx];
			I[3] = valonplanes[(z - 1) +	y*nny +		x*nnx];
			I[7] = valonplanes[(z) +		y*nny +		x*nnx];
			break;
	}
	isnearborders(x, y, z, side, corner);
	nearside(side, I);
	nearcorner(corner, I);
	return;
}

void EC::show() const{
	int ix, iy, iz, nx, ny, nz, nnx, nny, n1, n2;
	
	nx = size.x;
	ny = size.y;
	nz = size.z;
	nnx = nz*ny;
	nny = nz;
	switch(coaxis){
		case 'X':
			n1 = size.y;
			n2 = size.z;
			break;
		case 'Y':
			n1 = size.x;
			n2 = size.z;
			break;
		case 'Z':
			n1 = size.x;
			n2 = size.y;
			break;
	}
	printf("faces %c\n", coaxis);
	for(ix = 0; ix < nx; ix++){
		for(iy = 0; iy < ny; iy++){
			for(iz = 0; iz < (nz ); iz++){
				printf("%.2lf\t", valonplanes[iz + iy*nny + ix*nnx]);
			}
			printf("\n");
		}
		printf("\n");
	}
	
	printf("grid: axis1\n");
	for(int i = 0; i <= n1; i++ ){
		printf("%lf\t", axis1[i]);
	}
	printf("\n");
	printf("grid: axis2\n");
	for(int i = 0; i <= n2; i++ ){
		printf("%lf\t", axis2[i]);
	}
	printf("\n");
	return;
}
class space{
	coord size;
	int* layers;
	bool containsbadlayers;
	int pointer;
	coord pointer3d;
public: 
	coord getsize(){
		coord q;
		q.set(size.x, size.y, size.z);
		return q;
	}
	space(const string& gridfname);
	space(const space& s);
	void findbadlaysnear(int x, int y, int z, int alongaxis, int* badlaysdisposition) const;
	bool filled(){
		if (pointer == (size.x - 2)*(size.y - 2)*(size.z - 2)){
			return true;
		}
		else{
			return false;
		}
	}
	int layerat(int x, int y, int z){
		return layers[pointer3d.z + pointer3d.y*size.z + pointer3d.x*size.z*size.y];
	}
	void fillsomecells(int cellnumber);
	void skipsomecells(int cellnumber);
	void show() const;
	void fillspace(const string& optlaysfname, const string& cellfname);
	void resetbadlays(int* badlaysdisposition) const;
};
space::space(const space& s){
	size.set(s.size);
	layers = new int[size.x*size.y*size.z];
	int nny, nnx;
	nny  = size.z;
	nnx = size.z*size.y;
	for(int i = 0; i < size.x*size.y*size.z; i++){
		layers[i] = s.layers[i];
	}
	pointer = 0;
	pointer3d.set(0);
	containsbadlayers = s.containsbadlayers;
}
space::space(const string& gridfname){	
	size.set(getsizefromgrid(gridfname));
	layers = new int[size.x*size.y*size.z];
	fillarray(layers, size.x*size.y*size.z);
	pointer = 0;
	pointer3d.set(0);		
	containsbadlayers = false;		
}
void space::resetbadlays(int* badlaysdisposition) const{
	fillarray(badlaysdisposition, 4);
	return;
}
void space::findbadlaysnear(int x, int y, int z, int alongaxis, int* badlaysdisposition) const{
	if(!containsbadlayers){return;}
	int nnx, nny;
	nny = size.z;
	nnx = size.z*size.y;
	switch(alongaxis){
		case 'X':
			if(layers[(z - 1)+	(y - 1)*nny + (x)*nnx] == -2){
				badlaysdisposition[3] = 1;}
			if(layers[(z - 1) + (y)*nny +		(x)*nnx] == -2){badlaysdisposition[2] = 1;}
			if(layers[z +		(y - 1)*nny + (x)*nnx] == -2){badlaysdisposition[1] = 1;}
			if(layers[z +		(y)*nny +		(x)*nnx] == -2){badlaysdisposition[0] = 1;}
			break;
		case 'Y':
			if(layers[z - 1 + (y)*nny + (x - 1)*nnx] == -2){badlaysdisposition[3] = 1;}
			if(layers[(z - 1) + (y)*nny + (x)*nnx] == -2){badlaysdisposition[2] = 1;}
			if(layers[(z) + (y)*nny + (x - 1)*nnx] == -2){badlaysdisposition[1] = 1;}
			if(layers[(z) + (y)*nny + (x)*nnx] == -2){badlaysdisposition[0] = 1;}
			
			break;
		case 'Z':
			if(layers[(z) +	(y - 1)*nny + (x - 1)*nnx] == -2){badlaysdisposition[3] = 1;}
			if(layers[(z) + (y - 1)*nny + (x)*nnx]		== -2){badlaysdisposition[2]= 1;}
			if(layers[(z) +	(y)*nny +		(x - 1)*nnx] == -2){badlaysdisposition[1] = 1;}
			if(layers[(z) + (y)*nny +		(x)*nnx]		== -2){badlaysdisposition[0] = 1;}
			
			break;
	}
	return;
}

void space::show() const{
		int ix, iy, iz, nnx, nny;
		nnx = size.z*size.y;
		nny = size.z;
		printf("layers in space\n");
		for(ix = 0; ix < size.x; ix++){
			for(iy = 0; iy < size.y; iy++){
				for(iz = 0; iz < size.z; iz++){
					printf("%d\t", layers[iz + iy*nny + ix*nnx]);
				}
				printf("\n");
			}
			printf("\n");
		}
	}
void space::skipsomecells(int cellnumber){
		int ix, iy, iz, nnx, nny;
		ix = pointer3d.x;
		iy = pointer3d.y;
		iz = pointer3d.z;
		nnx = size.z*size.y;
		nny = size.z;
		cellnumber = pointer + cellnumber;
		while(pointer < cellnumber){
			pointer++;			
			if(iz == size.z - 2){
				if(iy == size.y - 2){
					ix++;
					iy = 1;
					iz = 1;
				}
				else{
					iy++;
					iz = 1;
				}
			}
			else{
				iz++;
			}
		}
		pointer3d.x = ix;
		pointer3d.y = iy;
		pointer3d.z = iz;
	}
void space::fillsomecells(int cellnumber){
		int ix, iy, iz, nz, ny, nnx, nny;
		ix = pointer3d.x;
		iy = pointer3d.y;
		iz = pointer3d.z;
		nz = size.z;
		ny = size.y;
		nnx = size.z*size.y;
		nny = size.z;
		cellnumber = pointer + cellnumber;
		while(pointer < cellnumber){
			layers[iz + iy*nny + ix*nnx] = -2;
			pointer++;			
			if(iz == nz - 2){
				if(iy == ny - 2){
					ix++;
					iy = 1;
					iz = 1;
				}
				else{
					iy++;
					iz = 1;
				}
			}
			else{
				iz++;
			}
		}
		pointer3d.x = ix;
		pointer3d.y = iy;
		pointer3d.z = iz;
	}
void space::fillspace(const string& optlaysfname, const string& cellfname){
	int layer, lpoints;
	vector <int> badlayers;
	badlayers = readbadlayers(optlaysfname);
	ifstream cell(cellfname.c_str());
	if(!badlayers.empty()){
		if(cell.is_open()){
			pointer3d.set(1);
			string s;
			while(getline(cell, s)){
				stringstream ss1;
				ss1.str(s);
				ss1>>layer;
				getline(cell,s);
				stringstream ss2;
				ss2.str(s);
				ss2>>lpoints;
				if(find(badlayers.begin(), badlayers.end(), layer) == badlayers.end()){	
					skipsomecells(lpoints);
				}
				else{
					if(!containsbadlayers){containsbadlayers = true;}
					fillsomecells(lpoints);					
				}
				ss1.clear();
				ss2.clear();
			}
			cell.clear();
			if(!filled()){
				showerror(cellfname);
			}
		}
		else{
			showerror(cellfname);
		}
	}
	return;
}


void fillallECs(EC& ecx, EC& ecy, EC& ecz, const string& particlefname, const string& normalfname, int particle){
	if(particle == 3){return;}
	double val;
	
	string s1, s2, s3;
	coord sizex, sizey, sizez;

	double iplus;
	int n, ix,  iy, iz, xnny, xnnx, ynny, ynnx, znny, znnx;

	sizex.set(ecx.getsize());
	sizey.set(ecy.getsize());
	sizez.set(ecz.getsize());
	xnny = sizex.z;
	xnnx = sizex.z*sizex.y;
	ynny = sizey.z;
	ynnx = sizey.y*sizey.z;
	znny = sizez.z;
	znnx = sizez.y*sizez.z;

	ifstream par(particlefname.c_str());

	ifstream norm(normalfname.c_str());
	
	if(par && norm){
		while(getline(par, s1) && getline(norm, s3)) {
			stringstream ss1;
			
			stringstream ss3;
			ss1.str(s1);
			
			ss3.str(s3);
			ss1>>iplus;
			if(particle == 1){iplus = (-1.0)*iplus;}
			ss3>>n>>ix>>iy>>iz;
			switch(abs(n)){
				case 1:
					val = (iplus);
					try{
						ecx.setEC(ix, iy, iz, n, val);
					}
					catch(const invalid_argument& e){
						throw(e);
					}
					break;
				case 2:
					val = (iplus);
					try{
						ecy.setEC(ix, iy, iz, n, val);
					}
					catch(const invalid_argument& e){
						throw(e);
					}
					break;
				case 3:
					val = (iplus);
					try{
						ecz.setEC(ix, iy, iz, n, val);
					}
					catch(const invalid_argument& e){
						throw(e);
					}
					break;
			}
			ss1.clear();
			ss3.clear();
		}
	}
	else{
		if(!par){showerror(particlefname);}
		if(!norm){showerror(normalfname);}
	}
	return;
}

void fillECgrid(EC& ec, const string& gridfname){
	int axis;
	int i;
	ifstream grid(gridfname.c_str());
	if(grid){
		string s;
		string s1, s2;
		axis = ec.getaxis();
		switch(axis){
			case 'X':
				i = 0;
				while(getline(grid, s) && i < 13){
					i++;
					if(i == 10 ){
						s1.assign(s);
					}
					if(i == 13){
						s2.assign(s);
					}
				}
				break;
			case 'Y':
				i = 0;
				while(getline(grid, s) && i < 13){
					i++;
					if(i == 7){
						//s1.resi
						s1.assign(s);
					}
					if(i == 13){
						s2.assign(s);
					}
				}
				break;
			case 'Z':
				i = 0;
				while(getline(grid, s) && i < 10){
					i++;
					if(i == 7){
						s1.assign(s);					
					}
					if(i == 10 ){
						s2.assign(s);					
					}
				}
				break;

		}
		
		ec.setgrid(s1,s2);	
	}
	else{
		showerror(gridfname);
		return;
	}
	grid.clear();
	return;
}

void calccurrentonedges(const space& Sp1, const EC& ec){
	FILE* J;
	char fname[32];

	coord size;
	double tmp1, tmp2;
	int axis;
	int ix, iy, iz, nx, ny, nz, nny, nnx;
	double* d = new double[4];
	double* I = new double[8];
	int* disposition = new int[4];
	fillarray(d, 4);
	fillarray(I, 8);
	fillarray(disposition, 4);

	axis = ec.getaxis();
	sprintf(fname, "%s%c", "J", axis);
	J = fopen(fname, "w");
	fprintf(J, "ix,iy,iz,j%c[q/cm**2/c]\n", axis);

	size.set(ec.getsize());
	nny = size.z;
	nnx = size.z*size.y;
	nx = size.x;
	ny = size.y;
	nz = size.z;
	for (ix = 1; ix < nx; ix++){
		for(iy = 1; iy < ny; iy++){
			for(iz = 1; iz < nz; iz++){
				tmp1 = 0.0;
				tmp2 = 0.0;
				
				ec.getdelta(ix, iy, iz, d);
				ec.getI(ix, iy, iz, I);
				Sp1.resetbadlays(disposition);
				Sp1.findbadlaysnear(ix, iy, iz, axis, disposition);
				ec.badlayschangeI(disposition, I);
				
				tmp1 +=	(I[0] + I[4])*d[1]*d[3];
				tmp1 +=	(I[1] + I[5])*d[0]*d[3];
				tmp1 +=	(I[2] + I[6])*d[1]*d[2];
				tmp1 += (I[3] + I[7])*d[0]*d[2];
				tmp2 += 2.0*(d[0] + d[1])*(d[2] + d[3]);
				tmp1 = tmp1/tmp2;
				tmp1 = 4.80325021e-10*tmp1;
				switch (axis){
					case 'X':
						fprintf(J, "%d\t%d\t%d\t%e\n", ix, iy - 1, iz - 1, tmp1);
						break;
					case 'Y':
						fprintf(J, "%d\t%d\t%d\t%e\n", ix - 1, iy, iz - 1, tmp1);
						break;
					case 'Z':
						fprintf(J, "%d\t%d\t%d\t%e\n", ix - 1, iy - 1, iz, tmp1);
						break;
				}
			}
		}
	}
	
	fclose(J);
	delete[] d;
	delete[] I;
	delete[] disposition;
	return;
}


void createnormfname(const string& name, string& norm){
	norm.assign(name);
	norm.append(".lst");
	return;
}
void createresfname(const string& name, string& res){
	res.assign(name);
	res.append(".res");
	return;
}
void findname(string& s){
	int tmp1, tmp2;
	tmp1 = 65536;
	tmp2 = 65536;
	tmp1 = s.find_first_of(' ');
	tmp2 = s.find_first_of('\t');
	if(tmp1 > 0 && tmp2 > 0){
		if(tmp1 < tmp2){
			s.assign(s.substr(0, tmp1));
		}
		else{
			s.assign(s.substr(0, tmp2));
		}
	}
	else{
		if(tmp2 > 0){
			s.assign(s.substr(0, tmp2));
		}
		if(tmp1 > 0){
			s.assign(s.substr(0, tmp1));
		}
	}
	return;
}
void findfilenames(string &gridfname, string &cellfname, vector<string>& detectorsnames, vector<string>& normalsfnames, vector<string>& templatesfnames, vector<string>& resfnames, vector<int>& particlenames, vector<bool>& detecttypes){
	int i;
	bool old;
	string s, templatefname, normalfname;
	ifstream perenos("perenos");
	if(perenos){
		i = 0;
		while(getline(perenos, s) != NULL){
			i++;
			if(i == 5){
				gridfname.assign(s);
			}
			if(i == 6){
				cellfname.assign(s);
			}
		}
		/*
		normalfname.assign("faces.lst");
		
		particlenames.push_back(1);
		detectorsnames.push_back("EL");
		normalsfnames.push_back(normalfname);
		templatesfnames.push_back(templatefname);
		resfnames.push_back("faces.res");

		particlenames.push_back(2);
		detectorsnames.push_back("POS");
		normalsfnames.push_back(normalfname);
		templatesfnames.push_back(templatefname);
		resfnames.push_back("faces_pos.res");

		*/
	}

	ifstream detectors("detectors");
	if(detectors){
		bool old;
		i = 0;
		while(getline(detectors, s) != NULL){
			findname(s);
			string s1, s2;
			detectorsnames.push_back(s);
			createnormfname(s, s1);
			createresfname(s, s2);
			normalsfnames.push_back(s1);
			resfnames.push_back(s2);

			getline(detectors, s);
			{
				stringstream ss(s);
				int tmp;
				ss>>tmp;
				particlenames.push_back(tmp);
				ss.clear();
			}

			getline(detectors, s);
			{
				stringstream ss(s);				
				int tmp;
				ss>>tmp;
				if(tmp == 0){old = true;}
				if(tmp == 2){old = false;}
				detecttypes.push_back(old);
			}
			if(old){getline(detectors, s);}

			getline(detectors, s);
			getline(detectors, s);
			templatesfnames.push_back(s);
			getline(detectors, s);			
		}
		detectors.clear();
	}
	if(templatesfnames.empty()){throw invalid_argument("no detectors");}
	if(templatesfnames.size() != normalsfnames.size() || normalsfnames.size() != detectorsnames.size() || normalsfnames.size() != detectorsnames.size()) {
		throw invalid_argument("different numbers");
	}
	s.assign(string());
	templatefname.assign(string());
	normalfname.assign(string());
	return;
}

void findfilenames(vector<string>& detectorsnames, vector<string>& normalsfnames, vector<string>& templatesfnames, vector<string>& resfnames, vector<int>& particlenames){
	int i;
	
	string s, templatefname, normalfname;

	ifstream detectors("detectors");
	if(detectors){
		i = 0;
		while(getline(detectors, s) != NULL){
			i++;
			if(i%7 == 1){
				detectorsnames.push_back(s);
				s.append(".res");
				resfnames.push_back(s);
			}
			if(i%7 == 2){
				stringstream ss(s);
				int tmp;
				ss>>tmp;
				particlenames.push_back(tmp);
				ss.clear();
			}
			if(i%7 == 5){
				normalsfnames.push_back(s);
			}
			if(i%7 == 6){
				templatesfnames.push_back(s);
			}
		}
		detectors.clear();
	}
	if(templatesfnames.empty()){throw invalid_argument("no detectors");}
	if(templatesfnames.size() != normalsfnames.size() || normalsfnames.size() != detectorsnames.size() || normalsfnames.size() != detectorsnames.size()) {
		throw invalid_argument("different numbers");
	}
	s.assign(string());
	templatefname.assign(string());
	normalfname.assign(string());
	return;
}

void findfilenames(string &gridfname, string &cellfname){
	int i;
	string s;
	ifstream perenos("perenos");
	if(!perenos){showerror("perenos"); return;}

	i = 0;
	while(getline(perenos, s) != NULL){
		i++;
		if(i == 5){
			gridfname.assign(s);	
		}
		if(i == 6){
			cellfname.assign(s);
		}
	}
	s.assign(string());
	perenos.clear();
	return;
}

void findfilenames(string& optlaysfname){
	printf("input optional layers filename:\n");
	cin>>optlaysfname;
	return;
}

void readgridaxis(const string& gridfname, int& size,const string& axis){
	int n;
	string s;
	stringstream ss;
	ifstream grid(gridfname.c_str());
	while(getline(grid, s) != NULL && s[0] != axis[0]){	}
	getline(grid, s);
	ss.str(s);
	ss>>n;
	size = n;
	ss.clear();
	return;
}

void readgrid(double* coord, const string& gridfname,const int& size,const string& axis){
	double tmp;
	string s;
	stringstream ss;
	ifstream grid(gridfname.c_str());
	while(getline(grid, s) != NULL && s[0] != axis[0]){	}
	getline(grid, s);
	getline(grid, s);
	ss.str(s);
	for(int i = 0; i < size + 1; i++){
		ss>>tmp;
		coord[i] = tmp;
	}	
	ss.clear();
	return;
}

void convert(const string& oldnormfname, const string& gridfname, const string& normfname){
	string s;	
	double xx, yy, zz;
	int celx, cely, celz, nx, ny, nz, n;
	readgridaxis(gridfname, nx, "X");
	readgridaxis(gridfname, ny, "Y");
	readgridaxis(gridfname, nz, "Z");
	double* x = new double[nx + 1];
	double* y = new double[ny + 1];
	double* z = new double[nz + 1];
	readgrid(x, gridfname, nx, "X");
	readgrid(y, gridfname, ny, "Y");
	readgrid(z, gridfname, nz, "Z");
	ifstream oldnorm(oldnormfname.c_str());
	FILE* norm;
	norm = fopen(normfname.c_str(), "w");
	if(!oldnorm.is_open()){
		showerror();
	}
	if(!norm){
		showerror();
	}
	while(getline(oldnorm, s) != NULL){
		stringstream ss(s);
		ss>>n>>celx>>cely>>celz;
		norma nor(n);
		if(n == 1 || n == -1){
			xx = x[celx - 1];
		}
		else{
			xx = 0.5*(x[celx - 1] + x[celx]);
		}
		if(n == 2 || n == -2){
			yy = y[cely - 1];
		}
		else{
			yy = 0.5*(y[cely - 1] + y[cely]);
		}
		if(n == 3 || n == -3){
			zz = z[celz - 1];
		}
		else{
			zz = 0.5*(z[celz - 1] + z[celz]);
		}
		
		fprintf(norm, "%e\t", xx);
		fprintf(norm, "%e\t", yy);
		fprintf(norm, "%e\t", zz);
		fprintf(norm, "%d\t", nor.n1);
		fprintf(norm, "%d\t", nor.n2);
		fprintf(norm, "%d\t", nor.n3);
		fprintf(norm, "\n");
		ss.clear();
	}
	free(x);
	free(y);
	free(z);
	oldnorm.clear();
	fclose(norm);
	s.assign(string());
	return;
}


void readtempl(const string& templ, int& energyintervalsnum, bool& exists){
	int flag;
	string s;
	ifstream temp(templ.c_str());
	if(!temp){showerror(templ); return;}
	int i = 0;
	while(getline(temp, s) != NULL){
		i++;
		if(i == 2){
			stringstream ss(s);
			ss>>energyintervalsnum;
			ss.clear();
		}
		if(i ==  3){
			stringstream ss(s);
			ss>>flag;
			if(flag == 0){exists = true;}
			else{exists = false;}
			ss.clear();
		}
	}
	temp.clear();
	s.assign(string());
	return;
}



bool increase(int normprev, int normcur){
	if(normprev == normcur){return true;}
	if(abs(normprev) < abs(normcur)){return true;}
	if(normprev < normcur){return true;}
	return false;
}



int numberofnormal(int normal){
	switch(normal){
		case -1:
			return 0;
		case 1:
			return 1;
		case -2:
			return 2;
		case 2:
			return 3;
		case -3:
			return 4;
		case 3:
			return 5;
		default:
			return 0;
	}
}


int readdetect(const string& normfname, vector<int>& normchanges, bool old, int* fsizes){
	int k;
	ifstream norm(normfname.c_str());
	if(!norm){
		throw invalid_argument("no norm file");
		return 0;
	}
	string s;
	int i, iprev;
	int n1, n2, n3, ncur, nprev,n;
	double x, y, z;
	i = 0;
	k = 0;
	nprev = 0;
	while(getline(norm, s) != NULL){
		stringstream ss1(s);
		ss1>>x>>y>>z;
		ss1>>n1>>n2>>n3;
		norma nor(n1, n2, n3);
		ncur = nor.n;		
		if(nprev != ncur){
			
			if(i > 0){
				fsizes[numberofnormal(nprev)] += i - normchanges[k - 1];
			}	
			normchanges.push_back(i);
			k++;
		}		
		nprev = ncur;
		i++;
		ss1.clear();
	}
	normchanges.push_back(i);
	fsizes[numberofnormal(nprev)] += i - normchanges[k - 1];
	return i;
}

void askuser(const string& name, const int& detectnumber, bool& split){
	char answer1;

	cout<<"detector "<<name<<" has "<<detectnumber<<" detectors";
	cout<<'\n';
	cout<<"split output file?(y/n)\n";
	cin>>answer1;
	
	if(answer1 == 'y'){
		split = true;
	}
	else{		
		if(answer1 != 'n'){
			cout<<"no split";
		}
		split = false;
	}
	
}
void askuser(bool& calcEC){
	char answer1;

	cout<<"Calculate currents? (y/n)\n";
	cin>>answer1;
	
	if(answer1 == 'y'){
		calcEC = true;
	}
	else{		
		if(answer1 != 'n'){
			cout<<"not calculating";
		}
		calcEC = false;
	}
	
}

class detector{
	bool old;
	int familynumber;
	int detectnumber;
	int particle;

	int fsizes[6];

	vector <int> normchanges;

	string name;
	string resname;
	string norms;

	int energyintervalsnum;
	
public:
	bool exists;
	
	detector(int i, string detect, string norm, string templ, string res, int particle);
	detector(int i, string detect, string norm, string templ, string res, int particle, bool type);
	int getnumnormals();
	void calculateandwrite(ofstream& common);

};

detector::detector(int i, string detect, string norm, string templ, string res, int par){
	readtempl(templ, energyintervalsnum, exists);
	if(exists){
		old = true;
		familynumber = i;
		particle = par;
		name.assign(detect);
		resname.assign(res);
		norms.assign(norm);
		for(int q = 0; q < 6; q++){
			fsizes[q] = 0;
		}
		try{detectnumber = readdetect(norm, normchanges, old, fsizes);}
		catch(const invalid_argument& e){
			throw e;
		}
		//askuser(name, detectnumber, split);
	}
}
detector::detector(int i, string detect, string norm, string templ, string res, int par, bool type){
	readtempl(templ, energyintervalsnum, exists);
	if(exists){
		old = type;
		familynumber = i;
		particle = par;
		name.assign(detect);
		resname.assign(res);
		norms.assign(norm);
		for(int q = 0; q < 6; q++){
			fsizes[q] = 0;
		}
		try{detectnumber = readdetect(norm, normchanges, old, fsizes);}
		catch(const invalid_argument& e){
			throw e;
		}
	}
}



void makecap(ofstream& dest,int detectsnumber, int particle, int fnumber, int normalnumber,  int energyintervalsnum){
	string part_name;
	int i;
	if(particle == 1){part_name.assign("Электроны"); i = 1;}
	else{
		if(particle == 2){part_name.assign("Позитроны"); i = 0;}
		else{
			if(particle == 3){part_name.assign("Кванты"); i = 3;}
			else{part_name.assign("Strange particle");}
		}
	}
	dest<<part_name;
	dest<<'\n';
	dest<<"Номер спектра";
	dest<<'\n';
	dest<<100*fnumber+10*normalnumber + particle;
	dest<<'\n';
	dest<<"Мощность спектра (шт/см**2/с)- для поверхностных, (шт/см**3/с)- для объемных";
	dest<<'\n';
	dest<<1;
	dest<<'\n';
	dest<<"Тип спектра (0-фиксированный, 1-разыгрывание, 2-список, 3-детекторы)";
	dest<<'\n';
	dest<<3;
	dest<<'\n';
	dest<<"Число частиц";
	dest<<'\n';
	dest<<energyintervalsnum;
	dest<<'\n';
	dest<<"Количество детекторов";
	dest<<'\n';
	dest<<detectsnumber;
	dest<<'\n';
	dest<<"Энергия+нормаль";
	dest<<'\n';
	part_name.assign(string());
	return;
}

void makecap(ofstream& dest, int detectsnumber, int particle, int  j, int energyintervalsnum){
	string part_name;
	int i;
	if(particle == 1){part_name.assign("Электроны"); i = 1;}
	else{
		if(particle == 2){part_name.assign("Позитроны"); i = 0;}
		else{
			if(particle == 3){part_name.assign("Кванты"); i = 3;}
			else{part_name.assign("Strange particle");}
		}
	}
	dest<<part_name;
	dest<<'\n';
	dest<<"Номер спектра";
	dest<<'\n';
	dest<<60 + particle;
	dest<<'\n';
	dest<<"Мощность спектра (шт/см**2/с)- для поверхностных, (шт/см**3/с)- для объемных";
	dest<<'\n';
	dest<<1;
	dest<<'\n';
	dest<<"Тип спектра (0-фиксированный, 1-разыгрывание, 2-список, 3-детекторы)";
	dest<<'\n';
	dest<<3;
	dest<<'\n';
	dest<<"Число частиц";
	dest<<'\n';
	dest<<energyintervalsnum;
	dest<<'\n';
	dest<<"Количество детекторов";
	dest<<'\n';
	dest<<detectsnumber;
	dest<<'\n';
	dest<<"Энергия+нормаль";
	dest<<'\n';
	part_name.assign(string());
	return;
}
bool digital(const string& s){
	if (s.compare("0") == 0 ||
        s.compare("1") == 0 ||
        s.compare("2") == 0 ||
        s.compare("3") == 0 ||
        s.compare("4") == 0 ||
        s.compare("5") == 0 ||
        s.compare("6") == 0 ||
        s.compare("7") == 0 ||
        s.compare("8") == 0 ||
        s.compare("9") == 0)
    {
        return true;
    }

    return false;
}
bool hastwoindex(const string& name){
	int size = name.size();
	string s1, s2;
	s1.assign(name.substr(size - 1, 1));
	s2.assign(name.substr(size - 3, 1));
	if(digital(s1)){
		if(digital(s2)){
			return true;
		}
		else{
			return false;
		}
	}
	else{
		return false; //and throw exeption
	}
}
void createname(string& outputfname1, const string& name){
	outputfname1.assign(name);
	outputfname1.append(".spc");
	transform(outputfname1.begin(), outputfname1.end(), outputfname1.begin(), toupper);
	return;
}
void createname(string& outputfname, const string& name, int numb){
	outputfname.assign(name);
	stringstream tmpss2;
	tmpss2<<"_" << numb;
	string str2 = tmpss2.str();
	outputfname.append(str2);
	outputfname.append(".spc");
	transform(outputfname.begin(), outputfname.end(), outputfname.begin(), toupper);				
	return;
}


void writelog(ofstream& common,const string& fname, int findex, int nindex, int particle, int damount){
	common<<fname<<"\t";
	common<<findex<<nindex<<particle<<"\t";
	common<<damount<<" detectors\t";
	common<<"\n";
}
void detector::calculateandwrite(ofstream& log){
	double x, y, z, d1, d2, d3, d4, d5;
	int nnor, j, i, k, n1, n2, n3;
	string s, s1, s2;
	string outputfname1, outputfname;

	bool firsttimemeetnorm[6];
	for(int q = 0; q < 6; q++){
		firsttimemeetnorm[q] = true;
	}
	
	ifstream res(resname.c_str());
	ifstream norm(norms.c_str());
	ofstream dest;
	dest.open("nevermind");
	
	j = 0;
	i = 0;
	k = 0;

	for(int a = 0; a < 6; a++){
		createname(outputfname, name, a);
		remove(outputfname.c_str());
		outputfname.assign("");
	}

	log<<name<<'\n';
	int normchangesrange = normchanges.size();
	while(getline(norm, s1) != NULL){		
		stringstream ss1(s1);
		ss1>>x>>y>>z;
		ss1>>n1>>n2>>n3;
		norma nor(n1, n2, n3);
		nnor = nor.normindex();
		if(k < normchangesrange && j == normchanges[k]){			
			dest.clear();
			dest.close();
			
			createname(outputfname, name, nnor);
			dest.open(outputfname.c_str(), fstream::app|fstream::out);
			
			if(!dest.is_open()){cout<<"ERROR!!\n"; }
			dest.precision(10);
			k++;
			int detectsinfile;
			if(!old){
				detectsinfile = normchanges[k] - j;
				makecap(dest, detectsinfile, particle, familynumber, nnor, energyintervalsnum);
				writelog(log, outputfname, familynumber, nnor, particle, detectsinfile);
			}
			else{
				if(firsttimemeetnorm[nnor]){
					firsttimemeetnorm[nnor] = false;
					detectsinfile = fsizes[nnor]; 
					makecap(dest, detectsinfile, particle, familynumber, nnor, energyintervalsnum);
					writelog(log, outputfname, familynumber, nnor, particle, detectsinfile);
				}
			}	
			
		}		
		
		dest<<x<<'\t'<<y<<'\t'<<z<<'\t'<<nor.n<<'\n';
		for(int i = 0; i < energyintervalsnum; i++){
			getline(res, s2);
			stringstream ss2(s2);
			ss2>>d1>>d2>>d3>>d4>>d5;
			d1 = 0.001*d1;
			dest<<d1<<'\t'<<d2<<'\t'<<d3<<'\t'<<d4<<'\t'<<d5<<'\n';
		}
		getline(res, s2);
		j++;
	}
	dest.clear();
	dest.close();
	s.assign(string()); 
	remove(outputfname1.c_str());
	remove("nevermind");
	s1.assign(string());
	s2.assign(string());
	outputfname.assign(string());	
	return;
}

void findschemetype(bool& old){
	int i;
	string s;
	ifstream perenos("perenos");
	if(perenos){
		i = 0;
		while(getline(perenos, s) != NULL && i < 3){
			i++;
		}
		
		stringstream ss(s);
		int tmp;
		ss>>tmp;
		if(tmp == 0){
			old = true;
		}
		else{
			old = false;
		}
	}
	else{
		showerror("perenos");
	}
	return;
}
void newscheme(){
	string gridfname, cellfname, normalfname, templatefname, optlaysfname;	
	vector <string> detectorsfnames;
	vector <string> templatesfnames;
	vector <string> normalsfnames;
	vector <string> resfnames;
	vector <int> particlenames;
	vector <bool> detecttypes;

	try{
		findfilenames(gridfname, cellfname, detectorsfnames, normalsfnames, templatesfnames, resfnames, particlenames, detecttypes);
	}
	catch(const invalid_argument& e){
		cerr<<e.what()<<'\n';
		return;
	}
	
	int i, detectorsnumber;
	i = 0;
	detectorsnumber = detectorsfnames.size();
	vector <string>::iterator pdetect = detectorsfnames.begin();
	vector <string>::iterator pnorm = normalsfnames.begin();
	vector <string>::iterator ptempl = templatesfnames.begin();
	vector <string>::iterator pres = resfnames.begin();
	vector <int>::iterator pparticle = particlenames.begin();
	vector <bool>::iterator pdtype = detecttypes.begin();
	
	space Sp1(gridfname);	
	bool calcEC;
	askuser(calcEC);
	if(calcEC){
		findfilenames(optlaysfname);
	}
	Sp1.fillspace(optlaysfname, cellfname);	
	
	coord size;
	size.set(Sp1.getsize());
	EC ecx(size.x - 1,	size.y,		size.z,		'X');
	EC ecy(size.x,		size.y - 1, size.z,		'Y');
	EC ecz(size.x,		size.y,		size.z - 1, 'Z');	

	fillECgrid(ecx, gridfname);
	fillECgrid(ecy, gridfname);
	fillECgrid(ecz, gridfname);	
	
	ofstream common("commonfile");

	
	while(i < detectorsnumber){
		
		i++;
		cout<<i<<" of "<< detectorsnumber<<'\n';
		string totfname;
		totfname.assign(*pdetect);
		totfname.append(".tot");
		if(*pdtype){
			normalfname.assign( *pnorm);			
		}
		else{
			normalfname.assign(*pdetect);
			normalfname.append("normals");
			convert(*pnorm, gridfname, normalfname);
			if(calcEC){	
				fillallECs(ecx, ecy, ecz, totfname, *pnorm, *pparticle);
			}
		}					
		try{
			detector d(i, *pdetect, normalfname, *ptempl, *pres, *pparticle, *pdtype);
			if(d.exists){
				d.calculateandwrite(common);
			}
		}
		catch(const invalid_argument& e){
			cerr<<e.what()<<'\n';
			return;
		}
		if(!*pdtype){
			remove(normalfname.c_str());
		}
		++pdetect;
		++pnorm;
		++ptempl;
		++pres;
		++pparticle;
		++pdtype;
		
	}
	
	/*
	Sp1.show();
	ecx.show();
	ecy.show();
	ecz.show();
	*/
	if(calcEC){
		calccurrentonedges(Sp1, ecx);	
		calccurrentonedges(Sp1, ecy);
		calccurrentonedges(Sp1, ecz);
		
		common<<"JX"<<" "<<"JY"<<" "<<"JZ";
	}
	common.clear();
	common.close();
	
	return;
}


void oldscheme(){
	string gridfname, cellfname, templatefname, optlaysfname;	
	vector <string> detectorsfnames;
	vector <string> templatesfnames;
	vector <string> normalsfnames;
	vector <string> resfnames;
	vector <int> particlenames;

	try{
		findfilenames(detectorsfnames, normalsfnames, templatesfnames, resfnames, particlenames);
	}
	catch(const invalid_argument& e){
		cerr<<e.what()<<'\n';
		return;
	}
	
	ofstream common("commonfile");

	int i, detectorsnumber;
	i = 0;
	detectorsnumber = detectorsfnames.size();
	vector <string>::iterator pdetect = detectorsfnames.begin();
	vector <string>::iterator pnorm = normalsfnames.begin();
	vector <string>::iterator ptempl = templatesfnames.begin();
	vector <string>::iterator pres = resfnames.begin();
	vector <int>::iterator pparticle = particlenames.begin();	
	
	while(i < detectorsnumber){
		
		i++;
		cout<<i<<" of "<< detectorsnumber<<'\n';		
		try{
			detector d(i, *pdetect, *pnorm, *ptempl, *pres, *pparticle);
			if(d.exists){
				d.calculateandwrite(common);
			}
		}
		catch(const invalid_argument& e){
			cerr<<e.what()<<'\n';
			return;
		}
		++pdetect;
		++pnorm;
		++ptempl;
		++pres;
		++pparticle;
	}
	common.clear();
	common.close();
	return;
}


void scheme(){
	bool old;
	findschemetype(old);
	if(old){
		oldscheme();
	}
	else{
		newscheme();
	}

}


int main(){
	_CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF);
	_CrtSetDbgFlag(_CRTDBG_LEAK_CHECK_DF);
	scheme();
	_CrtDumpMemoryLeaks();
	return 0;
}