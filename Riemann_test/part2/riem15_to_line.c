#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define maxit 20
#define ee    1.e-7
#define estr  (1.+1.e-7)

#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))


FILE *fpr,*fpw;

       double rqp[3][2],r_o,u_o,p_o,om,s1,s2,s3;
       int iwpar;
       double g[2],ga[2],gb[2],gc[2],gg[2],c[2],rc[2],du,ff,df;

       void ff_df ( double p );
       void riemi ( int lom );
/*----------------------------main-------------------------------*/
    void main(void)
    {
    int l,nvar,n,lom,np;
    double x,dx,x0,xn,t;
    fpr=fopen("inp.dat","r");
    fpw=fopen("out.dat","w");
     fscanf(fpr,"%d ",&(nvar));
//     iwpar=0;
     for( n=1; n <= nvar; n++ )
       {
       fprintf(fpw,"Initial data \n");
       for ( l=0; l <=1; l++ )
         {
         fscanf(fpr,"%lf  %lf  %lf",&(rqp[0][l]),&(rqp[1][l]),&(rqp[2][l]) );
         fprintf(fpw,"%15.7e  %15.7e  %15.7e \n",rqp[0][l],rqp[1][l],rqp[2][l] );
         }
fscanf(fpr," %lf %lf %lf ",&x0,&xn,&t);
fscanf(fpr,"%d",&np);
         lom=0;
dx=(xn-x0)/(double)np;
for(l=0; l <= np; l++)
  {         
   x=x0+l*dx;
   om=x/t;
   riemi(lom);
   fprintf(fpw," %15.7e  %15.7e  %15.7e  %15.7e  \n",x,r_o,u_o,p_o);
  }
//       fprintf(fpw," \n  \n");
       }
     fclose(fpr);
     fclose(fpw);
     return;
     }
/*----------------------------riemi------------------------------*/
       void riemi (int lom )
       {
       static double ss[2][2],rr[2],cc[2];
       int i,it;
       double e_pmin,u_vc,p0,p1,pp,uu,ai,tmp1,tmp2;
       double del_x;

          for ( i=0; i <= 1; i++ )
          {
          g[i]=1.4;
          ga[i]=0.5*(g[i]+1.)/g[i];
          gb[i]=0.5*(g[i]-1.)/g[i];
          gc[i]=2./(g[i]-1.);
          gg[i]=1./g[i];
          }

       for ( i=0; i <= 1; i++ )
       {
       c[i]=sqrt(g[i]*rqp[2][i]/rqp[0][i]);
       rc[i]=rqp[0][i]*c[i];
       }
       e_pmin=min(rqp[2][0],rqp[2][1])*ee;

/*       check existance    */

       u_vc=-c[0]*gc[0]-c[1]*gc[1];
       du=rqp[1][0]-rqp[1][1];
       if(du <= u_vc)
        {
        printf(" vacuem - no solution  \n");
        exit (-1);
        }

/*      inicial guess  */

       p0=(rqp[2][0]*rc[1]+rqp[2][1]*rc[0]+du*rc[0]*rc[1])/(rc[0]+rc[1]);
       p0=max(p0, e_pmin);
       it=1;
start:
       if( it > maxit ) goto finish;
       ff_df(p0);
       p1=p0-ff/df;
//       fprintf(fpw,"p1=%lf  p0=%lf  ff=%lf  df=%lf  \n",p1,p0,ff,df);   
       if( fabs(p1-p0) > e_pmin )
                                 { p0=p1;  goto start;  }
finish:
       pp=p0;  uu=s2;
/*     wave coeff. and dens   */
       for ( i=0; i <= 1; i++ )
         {
         ai=2.*( (double) i -0.5);
         if(p0 < rqp[2][i]*estr)
           {
           ss[0][i]=rqp[1][i]+ai*c[i];
           cc[i]=c[i]-0.5*ai*(g[i]-1.)*(rqp[1][i]-uu);
           ss[1][i]=uu+ai*cc[i];
           rr[i]=g[i]*pp/(cc[i]*cc[i]);
           }
         else
           {
           tmp1=(g[i]+1.)*pp+(g[i]-1.)*rqp[2][i];
           tmp2=(g[i]-1.)*pp+(g[i]+1.)*rqp[2][i];
           rr[i]=rqp[0][i]*tmp1/tmp2;
           tmp1=rqp[0][i]-rr[i];
           tmp2=rqp[0][i]*rqp[1][i]-rr[i]*uu;
           ss[0][i]=tmp2/tmp1;
           ss[1][i]=ss[0][i];
           }
         }
  
/*     data on line  */

       if(lom == 1) om=ss[0][0];
       if(lom == 2) om=uu;
       if(lom == 3) om=ss[0][1];
       s1=ss[0][0];
       s3=ss[0][1];
/*           1          */
       if( om <= ss[0][0] )
         {   r_o=rqp[0][0];  u_o=rqp[1][0]; p_o=rqp[2][0];  goto rtn_o;  }
/*           2          */
       if( ( om > ss[0][0] ) && ( om <= ss[1][0]) )
         {
         tmp1=1./(g[0]+1.);
         cc[0]=(2.*c[0]+(g[0]-1.)*(rqp[1][0]-om))*tmp1;
         u_o=om+cc[0];
         r_o=rqp[0][0]*pow( (cc[0]/c[0]),gc[0]);
         p_o=cc[0]*cc[0]*r_o/g[0];
         goto rtn_o;
         }
/*          3           */
       if( (om > ss[1][0]) && (om <= s2) )
         {  r_o=rr[0]; u_o=uu; p_o=pp; goto rtn_o;  }
/*          4           */
       if( (om > s2) && (om <= ss[1][1]) )
         {  r_o=rr[1]; u_o=uu; p_o=pp;  goto rtn_o; }
/*          5           */
       if( (om > ss[1][1]) && (om <= ss[0][1]) )
         {
          tmp2=1./(g[1]+1.);
          cc[1]=(2.*c[1]+(g[1]-1.)*(om-rqp[1][1]))*tmp2;
          u_o=om-cc[1];
          r_o=rqp[0][1]*pow( (cc[1]/c[1]),gc[1] );
          p_o=cc[1]*cc[1]*r_o/g[1];
          goto rtn_o;
         }
/*        6             */
       if( om > ss[0][1] )
         { r_o=rqp[0][1]; u_o=rqp[1][1]; p_o=rqp[2][1];  goto rtn_o; }
rtn_o:
       return; 
       }
/*----------------------------ff_df------------------------------*/
       void  ff_df( double p )
       {
       int i;
       double pp,rc_1,tmp1,tmp2;
       double f[2];
       df=0.;
       for ( i=0; i <= 1; i++ )
         {
          pp=p/rqp[2][i];
          if(pp > 1.)
            {
             rc_1=1./rc[i];
             tmp1=1./sqrt(ga[i]*pp+gb[i]);
             f[i]=(p-rqp[2][i])*tmp1*rc_1;
             tmp2=(g[i]+1.)*pp+3.*g[i]-1.;
             df=df+0.25*tmp2*(tmp1*tmp1*tmp1)*rc_1*gg[i];
             }
          else
             {
             tmp1=pow(pp,gb[i]);
             f[i]=gc[i]*c[i]*(tmp1-1.);
             df=df+c[i]*tmp1*gg[i]/p;
             }
         }
       ff=f[0]+f[1]-du;
       s2=rqp[1][0]-f[0];
       return;
       }


