/*
CONTENTS
void inp_mesh(void);
void inp_gdf(void);
void deriv_lr(void);
void deriv_du(void);
void step(void);
void flow_lr(void);
void flow_du(void);
void new_gdf(void);
*/


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "time.h"

#include "func.h"

#define IN(i,j)   ( (i)*(JT+1)+(j) )
#define IC(i,j)   ( (i)*JT+(j) )

#define IC_3(i,j,k)   ( (k)*(IT+1)*(JT+1)+(i)*JT+(j) )
#define Ilr_3(i,j,k)  ( (k)*(IT+1)*(JT+1)+(i)*JT+(j) )
#define Idu_3(i,j,k)  ( (k)*(IT+1)*(JT+1)+(i)*(JT+1)+(j) )



#define  GAM  1.4


double *gdf,*x,*y,*fllr,*fldu,*nlr,*tlr,*ndu,*tdu,*vol,*dxflr,*dyflr,*dxfdu,*dyfdu;
double t,dt,t_max,coord,turb,Mach,Re_ref,u_ref,mu_ref;
double lxmin,lymin,mu_max,H_Z;
int n_stp,max_stp,inp_par,IT,JT;



//*****************************************************

int main(void)
{
 FILE *fpr;
 clock_t start_time, end_time;
 double t_run;

 fpr=fopen("par_i.dat","r");
 fscanf(fpr,"%d %d",&IT,&JT);
 fscanf(fpr,"%lf",&t_max);
 fscanf(fpr,"%d %d",&max_stp,&inp_par);
 fscanf(fpr,"%lf %lf",&Mach,&Re_ref);
 fscanf(fpr,"%lf %lf",&coord,&turb);
 fscanf(fpr,"%lf",&H_Z);
 fclose(fpr);

 t=0.;
 dt=1.e-6;

 inp_mesh();
 inp_gdf();
 n_stp=1;

 start_time = clock();

 while( t < t_max )
   {
   double lcell,ki,bi;

   printf(" N_step=%d t=%e dt=%e \n",n_stp,t,dt);
   step();

   lcell=min(lxmin,lymin);
   ki=lcell/(u_ref+u_ref/Mach);
   bi=0.25*lcell*lcell/mu_max;
   dt=min(bi,ki);
   dt=H_Z*dt;
   t+=dt;
   n_stp++;
   if(n_stp > max_stp) break;
   }

   end_time=clock();
   t_run=((double)end_time-(double)start_time)/CLOCKS_PER_SEC;
   printf("Time of run %e \n",t_run);

 out();
 out_tec();
}

//*****************************************************
void step(void)
{
 int i,j,k;

// BOUNDARY CONDITION:  ghost cells
// down
   for(i=2; i < IT-2; i++)
    {
     gdf[IC_3(i,0,0)]=gdf[IC_3(i,1,0)]= gdf[IC_3(i,2,0)];
     gdf[IC_3(i,0,1)]=gdf[IC_3(i,1,1)]=-gdf[IC_3(i,2,1)];
     gdf[IC_3(i,0,2)]=gdf[IC_3(i,1,2)]=-gdf[IC_3(i,2,2)];
     gdf[IC_3(i,0,3)]=gdf[IC_3(i,1,3)]= gdf[IC_3(i,2,3)];
     gdf[IC_3(i,0,4)]=gdf[IC_3(i,1,4)]= gdf[IC_3(i,2,4)];
    }
// right
   for(j=0; j < JT; j++)
    for(k=0; k < 5; k++)
     gdf[IC_3(IT-1,j,k)]=gdf[IC_3(IT-2,j,k)]=gdf[IC_3(IT-3,j,k)];
// left
   for(j=0; j < JT; j++)
    {
    for(k=0; k < 5; k++)
     gdf[IC_3(0,j,k)]=gdf[IC_3(1,j,k)]=gdf[IC_3(2,j,k)];
    gdf[IC_3(0,j,2)]=gdf[IC_3(1,j,2)]=-gdf[IC_3(2,j,2)];
    }


 flow_lr();
 flow_du();
 new_gdf();

}

//*****************************************************

void inp_mesh(void)
{
 int i,j;
 double dx,dy;
 FILE *fpr;

 fpr=fopen("mesh.dat","r");
// inp coord 
 fscanf(fpr,"%d %d",&IT,&JT);
 memo_on();
 for(i=2; i <= IT-2; i++)
  for(j=2; j <= JT-2; j++)
   fscanf(fpr,"%lf %lf",&x[IN(i,j)],&y[IN(i,j)]);

 lxmin=9.9e+9; lymin=9.9e+9;

// norms left_right
 for(i=2; i <= IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    dx=x[IN(i,j+1)]-x[IN(i,j)];
    dy=y[IN(i,j+1)]-y[IN(i,j)];
    nlr[Ilr_3(i,j,0)]=dy;
    nlr[Ilr_3(i,j,1)]=-dx;
    nlr[Ilr_3(i,j,2)]=sqrt(dx*dx+dy*dy);
    lxmin=min(lxmin,nlr[Ilr_3(i,j,2)]);
    tlr[Ilr_3(i,j,0)]=dx;
    tlr[Ilr_3(i,j,1)]=dy;
   }
// norms down_up
 for(i=2; i < IT-2; i++)
  for(j=2; j <= JT-2; j++)
   {
    dx=x[IN(i+1,j)]-x[IN(i,j)];
    dy=y[IN(i+1,j)]-y[IN(i,j)];
    ndu[Idu_3(i,j,0)]=-dy;
    ndu[Idu_3(i,j,1)]=dx;
    ndu[Idu_3(i,j,2)]=sqrt(dx*dx+dy*dy);
    lymin=min(lymin,ndu[Idu_3(i,j,2)]);
    tdu[Idu_3(i,j,0)]=dx;
    tdu[Idu_3(i,j,1)]=dy;
   }
// volumes
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    dx=x[IN(i+1,j+1)]-x[IN(i,j)];
    dy=y[IN(i,j+1)]-y[IN(i+1,j)];
    vol[IC(i,j)]=0.5*dx*dy;
    dx=x[IN(i,j+1)]-x[IN(i+1,j)];
    dy=y[IN(i+1,j+1)]-y[IN(i,j)];
    vol[IC(i,j)]=vol[IC(i,j)]-0.5*dx*dy;
   }
 fclose(fpr);
}

//*****************************************************
void inp_gdf(void)
{
 int i,j,k;
 double d,u,v,p,e;
 double L_ref;
 FILE *fpr;
 d=1.; p=1.;  v=0.;
 u_ref=u=Mach*sqrt(GAM*p/d);
// !!!!!!!!!!!!!!!!
 L_ref=1.;
 mu_ref=u*L_ref*d/Re_ref;
 for(i=0; i < IT; i++)
  for(j=0; j < JT; j++)
   {
    gdf[IC_3(i,j,0)]=d;
    gdf[IC_3(i,j,1)]=u;
    gdf[IC_3(i,j,2)]=v;
    gdf[IC_3(i,j,3)]=p;
    gdf[IC_3(i,j,4)]=p/d;
    gdf[IC_3(i,j,5)]=0.;
   }
 if(inp_par == 1)
  {
   fpr=fopen("reg_i.dat","rb");
   fread(&IT,sizeof(int),1,fpr);
   fread(&JT,sizeof(int),1,fpr);
   fread(&t,sizeof(double),1,fpr);
   fread(&dt,sizeof(double),1,fpr);
   for(i=2; i < IT-2; i++)
    for(j=2; j < JT-2; j++)
     for(k=0; k < 5; k++)
      fread(&gdf[IC_3(i,j,k)],sizeof(double),1,fpr);
   fclose(fpr);  
  }
// BOUNDARY CONDITION:  ghost cells
// down
   for(i=2; i < IT-2; i++)
    {
     gdf[IC_3(i,0,0)]=gdf[IC_3(i,1,0)]= gdf[IC_3(i,2,0)];
     gdf[IC_3(i,0,1)]=gdf[IC_3(i,1,1)]=-gdf[IC_3(i,2,1)];
     gdf[IC_3(i,0,2)]=gdf[IC_3(i,1,2)]=-gdf[IC_3(i,2,2)];
     gdf[IC_3(i,0,3)]=gdf[IC_3(i,1,3)]= gdf[IC_3(i,2,3)];
     gdf[IC_3(i,0,4)]=gdf[IC_3(i,1,4)]= gdf[IC_3(i,2,4)];
    }
// right
   for(j=0; j < JT; j++)
    for(k=0; k < 5; k++)
     gdf[IC_3(IT-1,j,k)]=gdf[IC_3(IT-2,j,k)]=gdf[IC_3(IT-3,j,k)];
// left
   for(j=0; j < JT; j++)
    {
    for(k=0; k < 5; k++)
     gdf[IC_3(0,j,k)]=gdf[IC_3(1,j,k)]=gdf[IC_3(2,j,k)];
    gdf[IC_3(0,j,2)]=gdf[IC_3(1,j,2)]=-gdf[IC_3(2,j,2)];
    }
}


//*****************************************************

void flow_lr(void)
{
 int i,j,k,l,lom,alarm;
 double dx,dy,db,ub,vb,pb,e,NV1,TV1,NV2,TV2;
 double t11,t12,t21,t22,q1,q2;
 double mum,mut,muef,kapm,kapt,kapef;
 double dfx[5],dfy[5],ff[5];
 double xc1,yc1,xc2,yc2;
 double par[32],rqp[3][2],r_o,u_o,v_o,p_o,om,s1,s2,s3;
 
 lom=0;  om=0.; alarm=0;
 mu_max=0.;
               
 for(i=2; i <= IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
// left  values for riem
    for(k=0; k < 5; k++)
     {
       dfx[k]=0.;
       dfy[k]=0.;
     ff[k]=gdf[IC_3(i-1,j,k)]+0.5*dx*dfx[k]+0.5*dy*dfy[k];
     }
     NV1=(nlr[Ilr_3(i,j,0)]*ff[1]+nlr[Ilr_3(i,j,1)]*ff[2])/nlr[Ilr_3(i,j,2)];
     TV1=(tlr[Ilr_3(i,j,0)]*ff[1]+tlr[Ilr_3(i,j,1)]*ff[2])/nlr[Ilr_3(i,j,2)];
     rqp[0][0]=ff[0];
     rqp[1][0]=NV1;
     rqp[2][0]=ff[3];
// right values for riem
    for(k=0; k < 5; k++)
     {
       dfx[k]=0.;
       dfy[k]=0.;
     ff[k]=gdf[IC_3(i,j,k)]-0.5*dx*dfx[k]-0.5*dy*dfy[k];
     }
     NV2=(nlr[Ilr_3(i,j,0)]*ff[1]+nlr[Ilr_3(i,j,1)]*ff[2])/nlr[Ilr_3(i,j,2)];
     TV2=(tlr[Ilr_3(i,j,0)]*ff[1]+tlr[Ilr_3(i,j,1)]*ff[2])/nlr[Ilr_3(i,j,2)];
     rqp[0][1]=ff[0];
     rqp[1][1]=NV2;
     rqp[2][1]=ff[3];
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(i==2)
      {
       rqp[0][0]=rqp[0][1];
       rqp[1][0]=-rqp[1][1];
       rqp[2][0]=rqp[2][1];
      }



       for(l=0; l < 2; l++)
        for(k=0; k < 3; k++)
         par[k+3*l]=rqp[k][l]; 

       riemi(par,&alarm);
        if(alarm )
          {
          printf("ALARM_RIEM FROM l_r : i=%d j=%d  \n",i,j);
          printf("L: %e %e %e \n",rqp[0][0],rqp[1][0],rqp[2][0]);
          printf("R: %e %e %e \n",rqp[0][1],rqp[1][1],rqp[2][1]);
          exit(-1);
          }

     r_o=par[6];  u_o=par[7] ;   p_o=par[8] ;
     s1 =par[9];  s2 =par[10];   s3 =par[11]; 

     v_o=TV1;
     if(s2 < om) v_o=TV2; 

     db=r_o; pb=p_o;
     ub=(u_o*nlr[Ilr_3(i,j,0)]+v_o*tlr[Ilr_3(i,j,0)])/nlr[Ilr_3(i,j,2)];
     vb=(u_o*nlr[Ilr_3(i,j,1)]+v_o*tlr[Ilr_3(i,j,1)])/nlr[Ilr_3(i,j,2)];

     e=pb/(GAM-1.)+0.5*db*(ub*ub+vb*vb);

// Eu_flux
     fllr[Ilr_3(i,j,0)]=db*ub*nlr[Ilr_3(i,j,0)]+db*vb*nlr[Ilr_3(i,j,1)];
     fllr[Ilr_3(i,j,1)]=(db*ub*ub+pb)*nlr[Ilr_3(i,j,0)]+db*ub*vb*nlr[Ilr_3(i,j,1)];
     fllr[Ilr_3(i,j,2)]=db*ub*vb*nlr[Ilr_3(i,j,0)]+(db*vb*vb+pb)*nlr[Ilr_3(i,j,1)];
     fllr[Ilr_3(i,j,3)]=ub*(e+pb)*nlr[Ilr_3(i,j,0)]+vb*(e+pb)*nlr[Ilr_3(i,j,1)];

   }
}


//*****************************************************

void flow_du(void)
{
 int i,j,k,l,lom,alarm;
 double dx,dy,db,ub,vb,pb,e,NV1,TV1,NV2,TV2;
 double t11,t12,t21,t22,q1,q2;
 double mum,mut,muef,kapm,kapt,kapef;
 static double dfx[5],dfy[5],ff[5];
 double xc3,yc3,xc4,yc4;
 static double par[32],rqp[3][2],r_o,u_o,v_o,p_o,om,s1,s2,s3;

 lom=0;  om=0.;
 mu_max=0.;

 for(i=2; i < IT-2; i++)
  for(j=2; j <= JT-2; j++)
   {
// left(down)  values for riem
    for(k=0; k < 5; k++)
     {
       dfx[k]=0.;
       dfy[k]=0.;
     ff[k]=gdf[IC_3(i,j-1,k)]+0.5*dx*dfx[k]+0.5*dy*dfy[k];
     }
     NV1=(ndu[Idu_3(i,j,0)]*ff[1]+ndu[Idu_3(i,j,1)]*ff[2])/ndu[Idu_3(i,j,2)];
     TV1=(tdu[Idu_3(i,j,0)]*ff[1]+tdu[Idu_3(i,j,1)]*ff[2])/ndu[Idu_3(i,j,2)];
     rqp[0][0]=ff[0];
     rqp[1][0]=NV1;
     rqp[2][0]=ff[3];
// right (up) values for riem
    for(k=0; k < 5; k++)
     {
       dfx[k]=0.;
       dfy[k]=0.;
     ff[k]=gdf[IC_3(i,j,k)]-0.5*dx*dfx[k]-0.5*dy*dfy[k];
     }

     NV2=(ndu[Idu_3(i,j,0)]*ff[1]+ndu[Idu_3(i,j,1)]*ff[2])/ndu[Idu_3(i,j,2)];
     TV2=(tdu[Idu_3(i,j,0)]*ff[1]+tdu[Idu_3(i,j,1)]*ff[2])/ndu[Idu_3(i,j,2)];
     rqp[0][1]=ff[0];
     rqp[1][1]=NV2;
     rqp[2][1]=ff[3];
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(j==2)
      {
      NV1=-NV2;
      TV1=-TV2;
      rqp[0][0]=rqp[0][1];
      rqp[1][0]=-rqp[1][1];
      rqp[2][0]=rqp[2][1];
      }

       for(l=0; l < 2; l++)
        for(k=0; k < 3; k++)
         par[k+3*l]=rqp[k][l]; 

       riemi(par,&alarm);
        if(alarm )
          {
          printf("ALARM_RIEM FROM l_r : i=%d j=%d  \n",i,j);
          printf("L: %e %e %e \n",rqp[0][0],rqp[1][0],rqp[2][0]);
          printf("R: %e %e %e \n",rqp[0][1],rqp[1][1],rqp[2][1]);
          exit(-1);
          }

     r_o=par[6];  u_o=par[7] ;   p_o=par[8] ;
     s1 =par[9];  s2 =par[10];   s3 =par[11]; 

     v_o=TV1;
     if(s2  < om) v_o=TV2; 
     db=r_o; pb=p_o;
     ub=(u_o*ndu[Idu_3(i,j,0)]+v_o*tdu[Idu_3(i,j,0)])/ndu[Idu_3(i,j,2)];
     vb=(u_o*ndu[Idu_3(i,j,1)]+v_o*tdu[Idu_3(i,j,1)])/ndu[Idu_3(i,j,2)];
     e=pb/(GAM-1.)+0.5*db*(ub*ub+vb*vb);
 

// Eu_flux
     fldu[Idu_3(i,j,0)]=db*ub*ndu[Idu_3(i,j,0)]+db*vb*ndu[Idu_3(i,j,1)];
     fldu[Idu_3(i,j,1)]=(db*ub*ub+pb)*ndu[Idu_3(i,j,0)]+db*ub*vb*ndu[Idu_3(i,j,1)];
     fldu[Idu_3(i,j,2)]=db*ub*vb*ndu[Idu_3(i,j,0)]+(db*vb*vb+pb)*ndu[Idu_3(i,j,1)];
     fldu[Idu_3(i,j,3)]=ub*(e+pb)*ndu[Idu_3(i,j,0)]+vb*(e+pb)*ndu[Idu_3(i,j,1)];

   }
}


//*****************************************************

void new_gdf(void)
{
 int i,j,n,k;
 double p,e,yc;
 static double rhs[5],a0[5],a1[5];


 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    a0[0]=gdf[IC_3(i,j,0)];
    a0[1]=gdf[IC_3(i,j,0)]*gdf[IC_3(i,j,1)];
    a0[2]=gdf[IC_3(i,j,0)]*gdf[IC_3(i,j,2)];
    e=gdf[IC_3(i,j,3)]/(GAM-1.)+0.5*gdf[IC_3(i,j,0)]*
      (gdf[IC_3(i,j,1)]*gdf[IC_3(i,j,1)]+gdf[IC_3(i,j,2)]*gdf[IC_3(i,j,2)]);
    a0[3]=e;

//!!!!!!!!!!!!!!!!!!! for axesymm
    
    yc=0.25*(y[IN(i,j)]+y[IN(i+1,j)]+y[IN(i+1,j+1)]+y[IN(i,j+1)]);     

    rhs[0]=-a0[2];
    rhs[1]=-a0[2]*gdf[IC_3(i,j,1)];
    rhs[2]=-a0[2]*gdf[IC_3(i,j,2)];
    rhs[3]=-(e+gdf[IC_3(i,j,3)])*gdf[IC_3(i,j,2)];

    for(k=0; k < 4; k++)
     a1[k]=a0[k]-dt*(fllr[Ilr_3(i+1,j,k)]-fllr[Ilr_3(i,j,k)]+
                     fldu[Idu_3(i,j+1,k)]-fldu[Idu_3(i,j,k)])/vol[IC(i,j)]
           +dt*coord*rhs[k]/yc;

    gdf[IC_3(i,j,0)]=a1[0];
    gdf[IC_3(i,j,1)]=a1[1]/a1[0];
    gdf[IC_3(i,j,2)]=a1[2]/a1[0];
    p=(GAM-1.)*(a1[3]-0.5*(a1[1]*gdf[IC_3(i,j,1)]+a1[2]*gdf[IC_3(i,j,2)]));
    gdf[IC_3(i,j,3)]=p;
    gdf[IC_3(i,j,4)]=p/a1[0];
   }
}

