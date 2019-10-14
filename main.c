#include "global.h"
void make_lattice();
void over_relaxation();
void heatbath();
void initialize();
void bin();
void autocorr(double, double, double);
void fft3d(double *, double*,int);
void measure_sq();
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
  int j;
  for(i=0;i<WARMUP;i++){
    for(j=0;j<RELAX_STEPS;j++)
      over_relaxation();
    heatbath();
  }
  printf("warmed up \n");
  getchar();

  for(i=0;i<MCSTEPS;i++){
    for(j=0;j<RELAX_STEPS;j++)
      over_relaxation();
    heatbath();
    bin();
    autocorr(energy,energy2,0);
    if(i%SFACM==0){
            copyspins();
	    fft3d(sxqr,sxqi,1);
	    fft3d(syqr,syqi,1);
	    fft3d(szqr,szqi,1);
	    measure_sq();

    }
    //printf("%d\n",i);
  }
}
