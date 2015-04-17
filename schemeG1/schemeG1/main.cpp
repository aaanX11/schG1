
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
#include "space.h"
#include "riemi.h"

using namespace std;

void showerror(string& fname){
	cout<<"error reading file "<<fname<<"\n";
	exit(1);
	return;
}

void findfilenames(string &gridfname, string &cellfname, string &initvalfname, string &schparamfname, string &logfname, string &datafname, string& statefname){
	ifstream perenos("perenos");
	if(perenos){
		int i = 0;
		string s;
		while(getline(perenos, s) != NULL){
			i++;
			if(i == 5){
				gridfname.assign(s);
			}
			if(i == 6){
				cellfname.assign(s);
			}
		}
	}
	perenos.close();
	initvalfname.assign("initval");
	schparamfname.assign("schparam");
	logfname.assign("log");
	datafname.assign("data.dat");
	statefname.assign("state");
	return;
}

double inenergy(double rho, double pressure){
	double kappa = 1.4;
	return pressure/((kappa - 1)*rho);
}
double energyconverter(double rho, double u, double v, double w, double e){
	double kappa = 1.4;	
	return (kappa - 1)*(e - rho*0.5*(u*u + v*v +w*w));
}
void fluxes(double* fl, char* direction, double rho1, double u1, double p1, double rho2, double u2, double p2, double tan1, double tan2){
	fl[0] = rho1*u1;
	fl[5+0] = u2*rho2;
	//cout<<"in  fluxes\n";
	switch(direction[0]){
		case 'X':			
			fl[1] = rho1*u1*u1 + p1;
			fl[2] = rho1*u1*tan1;
			fl[3] = rho1*u1*tan2;

			fl[5+1] = rho2*u2*u2 + p2;
			fl[5+2] = rho2*u2*tan1;
			fl[5+3] = rho2*u2*tan2;
			break;
		case 'Y':
			fl[1] = rho1*u1*tan1;
			fl[2] = rho1*u1*u1 + p1;
			fl[3] = rho1*u1*tan2;

			fl[5+1] = rho2*u2*tan1;
			fl[5+2] = rho2*u2*u2 + p2;
			fl[5+3] = rho2*u2*tan2;
			break;
		case 'Z':
			fl[1] = rho1*u1*tan1;
			fl[2] = rho1*u1*tan2;
			fl[3] = rho1*u1*u1 + p1;

			fl[5+1] = rho2*u2*tan1;
			fl[5+2] = rho2*u2*tan2;
			fl[5+3] = rho2*u2*u2 + p2;
			break;
	}
	fl[4] = u1*(rho1*(inenergy(rho1,p1) + 0.5*(u1*u1+tan1*tan1+tan2*tan2)) + p1);
	fl[5+4] = u2*(rho2*(inenergy(rho2,p2) + 0.5*(u2*u2+tan1*tan1+tan2*tan2)) + p2);	
	return;
}

void step(space& sp, string& logfname, string& datafname){
	ofstream logfile(logfname.c_str(), std::ios_base::app);
	int nnx, nny, nx, ny, nz;
	nx = sp.size[0];
	ny = sp.size[1];
	nz = sp.size[2];
	nny = nz; nnx = ny*nz;
	//-----auxiliary big value storage------
	double *prevZ, *prevY, *prevX;
	//----------------fluxes----------------
	double *Qz, *Qy, *Qx, *sigma;
	double *par, r_o, u_o, p_o, s1, s2, s3, *lin;
	int alarm;
	double deltatx, deltaty, deltatz;
	
	lin = new double[6];
	par = new double[12];
	Qz = new double[10];
	Qy = new double[10];
	Qx = new double[10];
	sigma = new double[5];
	prevZ = new double[3];
	prevY = new double[3*(nz-2)];
	prevX = new double[3*(nz-2)*(ny-2)];
	//-------------towards -X----------------
	for(int iy = 1; iy < ny-1; iy++){
		for(int iz = 1; iz < nz-1; iz++){
			sp.getriemparam(par, "-X", 1, iy, iz);
			riemi(par, &alarm);
			prevX[3*(iz-1+(iy-1)*(nz-2))]= r_o = par[6];
			prevX[3*(iz-1+(iy-1)*(nz-2))+1] = u_o = par[7];
			prevX[3*(iz-1+(iy-1)*(nz-2))+2] = p_o = par[8] ;
			s1 =par[9];  s2 =par[10];   s3 =par[11];
		}
	}
	//----------------main cycle-------------
	deltatx = deltaty = deltatz = 10.;
	sp.savedata(datafname);
	for(int ix = 1; ix < nx-1; ix++){
		//-----------towards -Y--------------
		for(int iz = 1; iz < nz-1; iz++){
			sp.getriemparam(par, "-Y", ix, 1, iz);
			riemi(par, &alarm);
			s1 =par[9];  s2 =par[10];   s3 =par[11];
			prevY[3*(iz-1)] = r_o = par[6];
			prevY[3*(iz-1)+1] = u_o = par[7];
			prevY[3*(iz-1)+2] = p_o = par[8] ;
		}
		for(int iy = 1; iy < ny-1; iy++){
			//---------towards -Z------------
			sp.getriemparam(par, "-Z", ix, iy, 1);
			riemi(par, &alarm);
			prevZ[0] = r_o = par[6];
			prevZ[1] = u_o = par[7];
			prevZ[2] = p_o = par[8] ;
			for(int iz = 1; iz < nz-1; iz++){
				int index = iz + iy*nny + ix*nnx;
				if(sp.layers[index] == 0/*non obstacle cell (ix, iy, iz) proccessing*/){
					if(sp.layers[index - 1] == -1){	//cell (ix,iy, iz-1) belongs to an obstacle area 
													//so flux through -Z-face needs to be calculated
						sp.getriemparam(par, "-Z", ix, iy, iz);
						riemi(par, &alarm);
						prevZ[0] = r_o = par[6];
						prevZ[1] = u_o = par[7];
						prevZ[2] = p_o = par[8] ;}
					double RHO, U, V, W, P;
					RHO = sp.density[index];
					U = sp.xvelosity[index];
					V = sp.yvelosity[index];
					W = sp.zvelosity[index];
					P = sp.pressure[index];
					sp.getcellsize(lin, ix, iy, iz);
					sp.getriemparam(par, "Z", ix, iy, iz);
					riemi(par, &alarm);
					s1 =par[9];  s2 =par[10];   s3 =par[11];
					r_o = par[6];u_o = par[7]; p_o = par[8];
					sp.updatetimestep("Z", iz, s1, s3, deltatz);
					//fluxes
					fluxes(Qz, "Z", 
						prevZ[0], prevZ[1], prevZ[2], 
						r_o, u_o, p_o, 
						U, V);
					//auxiliary array update
					prevZ[0] = r_o;
					prevZ[1] = u_o;
					prevZ[2] = p_o;
					//-------towards Y-------
					if(sp.layers[index-nz] == -1){	//cell (ix,iy-1, iz) belongs to an obstacle area 
													//flux through -Y-face
						sp.getriemparam(par, "-Y", ix, iy, iz);
						riemi(par, &alarm);
						prevY[3*(iz-1)] = r_o = par[6];
						prevY[3*(iz-1)+1] = u_o = par[7];
						prevY[3*(iz-1)+2] = p_o = par[8] ;}
						sp.getriemparam(par, "Y", ix, iy, iz);
					riemi(par, &alarm);
					r_o=par[6];  u_o=par[7] ;   p_o=par[8] ;
					s1 =par[9];  s2 =par[10];   s3 =par[11];	
					sp.updatetimestep("Y", iy, s1, s3, deltaty);
					fluxes(Qy, "Y", 
						prevY[3*(iz-1)], prevY[3*(iz-1)+1], prevY[3*(iz-1)+2], 
						r_o, u_o, p_o, 
						U, W);
					prevY[3*(iz-1)] = r_o;
					prevY[3*(iz-1)+1] = u_o;
					prevY[3*(iz-1)+2] = p_o;
					//-------towards X-------
					if(sp.layers[index-nnx] == -1){//cell (ix-1,iy, iz) belongs to an obstacle area
						sp.getriemparam(par, "-X", ix, iy, iz);
						riemi(par, &alarm);
						prevX[3*(iz-1+(iy-1)*(nz-2))]= r_o = par[6];
						prevX[3*(iz-1+(iy-1)*(nz-2))+1] = u_o = par[7];
						prevX[3*(iz-1+(iy-1)*(nz-2))+2] = p_o = par[8];}
					sp.getriemparam(par, "X", ix, iy, iz);
					riemi(par, &alarm);
					r_o=par[6];  u_o=par[7] ;   p_o=par[8] ;
					s1 =par[9];  s2 =par[10];   s3 =par[11];
					sp.updatetimestep("X", ix, s1, s3, deltatx);
					fluxes(Qx, "X", 
						prevX[3*(iz-1+(iy-1)*(nz-2))], prevX[3*(iz-1+(iy-1)*(nz-2))+1], prevX[3*(iz-1+(iy-1)*(nz-2))+2], 
						r_o, u_o, p_o, 
						V, W);	
					prevX[3*(iz-1+(iy-1)*(nz-2))] = r_o;
					prevX[3*(iz-1+(iy-1)*(nz-2))+1] = u_o;
					prevX[3*(iz-1+(iy-1)*(nz-2))+2] = p_o;		
					//--------scheme---------
					sigma[0] = RHO;
					sigma[1] = RHO*U;
					sigma[2] = RHO*V;
					sigma[3] = RHO*W;
					sigma[4] = RHO*(inenergy(RHO,P) + 0.5*(U*U + V*V + W*W));
					for(int isig = 0; isig < 5; isig++){					
						sigma[isig] = 
							sigma[isig] + sp.deltat*(	
								(Qx[isig]-Qx[5+isig])/lin[0] 
								+ (Qy[isig]-Qy[5+isig])/lin[1] 
								+ (Qz[isig]-Qz[5+isig])/lin[2]		);}
					//-gas properties update-
					sp.density[index] = sigma[0];
					sp.xvelosity[index] = sigma[1]/sigma[0];
					sp.yvelosity[index] = sigma[2]/sigma[0];
					sp.zvelosity[index] = sigma[3]/sigma[0];
					sp.pressure[index] = energyconverter(sigma[0], sigma[1]/sigma[0], sigma[2]/sigma[0], sigma[3]/sigma[0], sigma[4]);}
			}
		}
	}
	sp.deltat = 1./(1./deltatx + 1./deltaty +1./deltatz);
	//sp.deltat = 0.09;
	//sp.showlog(logfile);
	//---------------------------------------
	delete[](lin);
	delete[](par);
	delete[](Qz);
	delete[](Qy);
	delete[](Qx);
	delete[](sigma);
	delete[](prevZ);
	delete[](prevY);
	delete[](prevX);
	logfile.close();
	return;

}

void getschparam(string schparamfname, int& nstepmax, double& tmax){
	ifstream schparam(schparamfname.c_str());
	if(schparam.is_open()){
		string s;
		getline(schparam, s);
		stringstream ss2(s);
		ss2>>tmax;
		getline(schparam, s);
		stringstream ss1(s);
		ss1>>nstepmax;
		schparam.close();
	}
	return;
}

int main(int argc, char* argv[]){
	double t, tmax;
	int nstep, nstepmax;
	string gridfname, cellfname, schparamfname, initvalfname, logfname, datafname, statefname;
	
	findfilenames(gridfname, cellfname, initvalfname, schparamfname, logfname, datafname, statefname);
	space sp(gridfname, argv[1]);
	if(argv[1][0] == 'c'){
		sp.restore(statefname, t, nstep);
	}
	else{
		sp.fillspace(cellfname, initvalfname, argv[1]);
		nstep = 0;
		t = 0.;
		
	}
	getschparam(schparamfname, nstepmax, tmax);	
	ofstream logfile(logfname.c_str());
	sp.showall(logfile);
	logfile.close();
	sp.deltat = 0.1;
	while(t < tmax && nstep < nstepmax){
		cout<<"step "<<nstep<<'\n';
		t = t + sp.deltat;
		step(sp, logfname, datafname);
		cout<<"time step "<<sp.deltat<<'\n';
		
		if(t < 1e-5){
			cout<<"t damaged\n";
		}
		if(sp.deltat < 1e-5){
			cout<<"deltast damaged\n";
		}
		nstep = nstep + 1;
		sp.saveslice(sp.size[0]/2, sp.size[2]/2, nstep);
		sp.savestate(statefname, t, tmax, nstep, nstepmax);
	}
	
	system("pause");
	return 0;
}