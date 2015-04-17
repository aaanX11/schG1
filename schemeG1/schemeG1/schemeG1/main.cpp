
#include <mpi.h>
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
	std::cerr<<"error reading file "<<fname<<"\n";
	exit(1);
	return;
}

void findfilenames(string &gridfname, string &cellfname, string &initvalfname, string &schparamfname, string &logfname, string &datafname, string& statefname){
	ifstream perenos("perenos");
	if(!perenos){
		std::cerr<<"function 'findfilenames' cannot open file perenos\n";
		return;
	}
	int i = 0;
	string s;
	while(getline(perenos, s) != NULL){
		i++;
		if(i == 5){
			s.erase( remove(s.begin(), s.end(), '\r'), s.end() );
			gridfname.assign(s);
		}
		if(i == 6){
			s.erase( remove(s.begin(), s.end(), '\r'), s.end() );
			cellfname.assign(s);
		}
	}
	ifstream grid(gridfname.c_str());
	if(!grid){
		std::cerr<<"function 'findfilenames' grid fails file "<<gridfname<<"\n";
		exit(1);
	}
	grid.close();
	perenos.close();
	initvalfname.assign("initval");
	schparamfname.assign("schparam");
	logfname.assign("log");
	datafname.assign("data");
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
void axisattr(int iaxis, int ix, int iy, int iz, int nx, int ny, int nz, int& stepbehind, int& stepforward, string& minusdirect, int& facebehind, string& direct, int& iter, int& tan1, int& tan2){
	int nny, nnx;
	nny = nz;
	nnx = nz*ny;
	switch(iaxis + 'X'){
		case 'X':
			stepbehind = -nnx;
			stepforward = nnx;
			direct.assign("X");
			minusdirect.assign("-X");
			facebehind = 5 + 5*(nz-2) + 5*(iz-1+(iy-1)*(nz-2));
			iter = ix;
			tan1 = 1;
			tan2 = 2;
			break;
		case 'Y':
			stepbehind = -nny;
			stepforward = nny;
			direct.assign("Y");
			minusdirect.assign("-Y");
			facebehind = 5 + 5*(iz-1);
			iter = iy;
			tan1 = 0;
			tan2 = 2;
			break;
		case 'Z':
			stepbehind = -1;
			stepforward =  1;
			direct.assign("Z");
			minusdirect.assign("-Z");
			facebehind = 0;
			iter = iz;
			tan1 = 0;
			tan2 = 1;
			break;
	}
	return;
}
void fluxes(double* fl, string& direction, double rho1, double u1, double p1, double ltan1, double ltan2, double rho2, double u2, double p2, double rtan1, double rtan2){
	fl[0] = rho1*u1;
	fl[5+0] = u2*rho2;
	//std::cerr<<"in  fluxes\n";
	switch(direction[0]){
		case 'X':			
			fl[1] = rho1*u1*u1 + p1;
			fl[2] = rho1*u1*ltan1;
			fl[3] = rho1*u1*ltan2;

			fl[5+1] = rho2*u2*u2 + p2;
			fl[5+2] = rho2*u2*rtan1;
			fl[5+3] = rho2*u2*rtan2;
			break;
		case 'Y':
			fl[1] = rho1*u1*ltan1;
			fl[2] = rho1*u1*u1 + p1;
			fl[3] = rho1*u1*ltan2;

			fl[5+1] = rho2*u2*rtan1;
			fl[5+2] = rho2*u2*u2 + p2;
			fl[5+3] = rho2*u2*rtan2;
			break;
		case 'Z':
			fl[1] = rho1*u1*ltan1;
			fl[2] = rho1*u1*ltan2;
			fl[3] = rho1*u1*u1 + p1;

			fl[5+1] = rho2*u2*rtan1;
			fl[5+2] = rho2*u2*rtan2;
			fl[5+3] = rho2*u2*u2 + p2;
			break;
	}
	fl[4] = u1*(rho1*(inenergy(rho1,p1) + 0.5*(u1*u1+ltan1*ltan1+ltan2*ltan2)) + p1);
	fl[5+4] = u2*(rho2*(inenergy(rho2,p2) + 0.5*(u2*u2+rtan1*rtan1+rtan2*rtan2)) + p2);	
	return;
}

void step(space& sp, string& logfname, string& datafname){
	ofstream logfile(logfname.c_str(), std::ios_base::app);
	string minusdirect, direct; 
	int nnx, nny, nx, ny, nz, 
		index, stepbehind, stepforward, facebehind, iter, tan1, tan2,
		alarm;
	double *bigval,
		*sigma, *Q,
		*par, r_o, u_o, p_o, s1, s2, s3,
		*lin, *deltat,
		RHO, *U, P, *U2;
	nx = sp.size[0]; ny = sp.size[1]; nz = sp.size[2];
	lin = new double[6];
	deltat = new double[3];
	par = new double[12];
	Q = new double[30];
	sigma = new double[5];
	bigval = new double[5 + 5*(nz-2) + 5*(nz-2)*(ny-2)];
	U = new double[6];
	U2 = new double[2];
	//----------------main cycle-------------
	deltat[0] = deltat[1] = deltat[2] = 10.;
	nny = nz; nnx = ny*nz;
	for(int ix = 1; ix < nx-1; ix++){
		for(int iy = 1; iy < ny-1; iy++){
			for(int iz = 1; iz < nz-1; iz++){
				index = iz + iy*nny + ix*nnx;
				if(sp.layers[index] == 0/*non obstacle cell (ix, iy, iz) proccessing*/){					
					RHO = sp.density[index];
					U[0] = sp.xvelosity[index];	U[1] = sp.yvelosity[index];	U[2] = sp.zvelosity[index];
					P = sp.pressure[index];
					sp.getcellsize(lin, ix, iy, iz);
					for(int iaxis = 0; iaxis < 3; iaxis++){
						axisattr(iaxis, ix, iy, iz, nx, ny, nz, stepbehind, stepforward, minusdirect,facebehind, direct, iter, tan1, tan2);
						if(sp.layers[index + stepbehind] == -1 || (iaxis == 1 && iy == 1)){	
							//cell behind belongs to an obstacle area 
							sp.getriemparam(par, minusdirect, ix, iy, iz);
							riemi(par, &alarm);
							if(alarm){
								std::cerr<<"function 'riemi' sends an alarm at ("<<ix<<" "<<iy<<" "<<iz<<'\n';}
							s1 =par[9];  s2 =par[10];   s3 =par[11];
							bigval[facebehind] = r_o = par[6];
							bigval[facebehind+1] = u_o = par[7];
							bigval[facebehind+2] = p_o = par[8] ;
							//sp.tanvel(&bigval[facebehind+3], s2, tan1, tan2, minusdirect, ix, iy, iz);
							bigval[facebehind+3] = U[tan1];
							bigval[facebehind+4] = U[tan2];							
						}
						if(iaxis == 1 && iy == ny -2 && sp.layers[index + stepforward] != -1){sp.getriemparam(par, "Yend", ix, iy, iz);}
						else{
						sp.getriemparam(par, direct, ix, iy, iz);}
						riemi(par, &alarm);
						if(alarm){
							std::cerr<<"function 'riemi' sends an alarm at ("<<ix<<" "<<iy<<" "<<iz<<'\n';}
						s1 =par[9];  s2 =par[10];   s3 =par[11];
						r_o = par[6];u_o = par[7]; p_o = par[8];
						if(iaxis == 1 && iy == ny -2){sp.tanvel(U2, s2, tan1, tan2, "Yend", ix, iy, iz);}
						else{
						sp.tanvel(U2, s2, tan1, tan2, direct, ix, iy, iz);}
						sp.updatetimestep(direct, iter, s1, s3, deltat[iaxis]);
						//fluxes
						fluxes(&Q[10*iaxis], direct, 
							bigval[facebehind], bigval[facebehind+1], bigval[facebehind+2], bigval[facebehind+3], bigval[facebehind+4],
							r_o, u_o, p_o, U2[0], U2[1]);
						//auxiliary array update
						bigval[facebehind] = r_o;
						bigval[facebehind+1] = u_o;
						bigval[facebehind+2] = p_o;
						bigval[facebehind+3] = U2[0];
						bigval[facebehind+4] = U2[1];
					}
					//--------scheme---------
					sigma[0] = RHO;
					sigma[1] = RHO*U[0];
					sigma[2] = RHO*U[1];
					sigma[3] = RHO*U[2];
					sigma[4] = RHO*(inenergy(RHO,P) + 0.5*(U[0]*U[0] + U[1]*U[1] + U[2]*U[2]));
					for(int isig = 0; isig < 5; isig++){					
						sigma[isig] = 
							sigma[isig] + sp.deltat*(	
								(Q[isig]-Q[5+isig])/lin[0] 
								+ (Q[10+isig]-Q[10+5+isig])/lin[1] 
								+ (Q[20+isig]-Q[20+5+isig])/lin[2]	)	+ sp.source[index];}
					//-gas properties update-
					sp.density[index] = sigma[0];
					sp.xvelosity[index] = sigma[1]/sigma[0];
					sp.yvelosity[index] = sigma[2]/sigma[0];
					sp.zvelosity[index] = sigma[3]/sigma[0];
					sp.pressure[index] = energyconverter(sigma[0], sigma[1]/sigma[0], sigma[2]/sigma[0], sigma[3]/sigma[0], sigma[4]);}
			}
		}
	}
	sp.deltat = 1./(1./deltat[0] + 1./deltat[1] +1./deltat[2]);
	logfile<<"deltat = "<<sp.deltat<<'\n';
	//---------------------------------------
	delete[](lin);
	delete[](U);
	delete[](U2);
	delete[](deltat);
	delete[](par);
	delete[](Q);
	delete[](sigma);
	delete[](bigval);
	logfile.close();
	return;

}

void getschparam(string schparamfname, int& nstepmax, double& tmax, int& frequency){
	ifstream schparam(schparamfname.c_str());
	if(!schparam){		std::cerr<<"function 'getchparam' cannot open file "<<schparamfname<<'\n';		return;	}
	string s;	
	getline(schparam, s);	stringstream ss2(s);	ss2>>tmax;		ss2.clear();	
	getline(schparam, s);	ss2.str(s);				ss2>>nstepmax;	ss2.clear();	
	getline(schparam, s);	ss2.str(s);				ss2>>frequency;
	schparam.close();
	std::cerr<<tmax<<'\t'<<nstepmax<<'\t'<<frequency<<'\n';
	return;
}

int main(int argc, char* argv[]){
	std::cerr<<"started\n";
	double t, tmax;
	int nstep, nstepmax, frequency, npro, myid;
	std::string gridfname, cellfname, schparamfname, initvalfname, logfname, datafname, statefname;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &npro);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	
	std::cerr<<"reading files\n";
	findfilenames(gridfname, cellfname, initvalfname, schparamfname, logfname, datafname, statefname);
	std::cerr<<"reading grid\n";
	space sp(gridfname, "ob");
	std::cerr<<"creating space\n";
	sp.fillspace(initvalfname);
	nstep = 0;
	t = 0.;
	getschparam(schparamfname, nstepmax, tmax,frequency);	
	sp.deltat = 0.1;
	while(t < tmax && nstep < nstepmax){
		std::cerr<<"step "<<nstep<<'\n';
		t = t + sp.deltat;
		if(nstep%frequency == 0){
			stringstream ss;
			string s1, s2;
			ss<<nstep;
			s1.append(datafname);
			s1.append(ss.str()); 
			s1.append(".dat");
			sp.savedata(s1);
			s2.assign(statefname);
			s2.append(ss.str());
			sp.savestate(s2, t, tmax, nstep, nstepmax);
		}
		step(sp, logfname, datafname);
		if(nstep == 0){sp.stopheat();}
		std::cerr<<"time step "<<sp.deltat<<'\n';		
		if(t < 1e-5){			std::cerr<<"t damaged\n";		}
		if(sp.deltat < 1e-5){			std::cerr<<"deltast damaged\n";		}
		nstep = nstep + 1;		
	}
	MPI_Finalize();
	return 0;
}