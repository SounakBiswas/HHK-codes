#include "global.h"
void make_lattice();
void over_relaxation();
void heatbath();
void initialize();
void bin();
void autocorr(double, double, double);
void fft3d(double *, double*,int);
void measure_sq();
void iterative_gs();
void copyspins(){
  int i;
  for(i=0;i<NSITES;i++){
    sxqr[i]=sx[i];
    syqr[i]=sy[i];
    szqr[i]=sz[i];
    sxqi[i]=0;
    syqi[i]=0;
    szqi[i]=0;
  
  }

}


void main(){
  initialize();
  make_lattice();
  int i;
  int j,k,flag;
  int idx;
  int x,y,z,cell,mcell,mx,my,mz;
  int strides[12]={0,12,23,33,42,50,57,63,68,72,75,77};
  FILE *fp;

  for(i=0;i<MCSTEPS;i++){
    for(j=0;j<RELAX_STEPS;j++)
      over_relaxation();
    iterative_gs();
    if(i%200==0)
      printf("%d %.16f\n",i,energy);
  }
    copyspins();
    fft3d(sxqr,sxqi,1);
    fft3d(syqr,syqi,1);
    fft3d(szqr,szqi,1);
    measure_sq();
    for(i=0;i<12;i++){
      for(j=i;j<12;j++){
	idx=strides[i]+j-i;
	sprintf(sfacname,"%si%dj%d.dat",sfacnamepref,i,j);
	fp=fopen(sfacname,"w");
	for(x=0;x<=lx;x++){
	  for(y=0;y<=ly;y++){
	    for(z=0;z<=lz;z++) { 
	      cell=x+y*lx+z*lx*ly;
	      mx=(mx==lx)?lx/2:(x+lx/2)%lx;
	      my=(my==ly)?ly/2:(y+ly/2)%ly;
	      mz=(mz==lz)?lz/2:(z+lz/2)%lz;
	      mcell=mx+lx*my+mz*lx*ly;

	      fprintf(fp,"%f %f\n",sqasqbre[idx*NCELLS+mcell],sqasqbim[idx*NCELLS+mcell]);
	    }
	  }
	}
      fclose(fp);

      }
    }

}
