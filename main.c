#include "global.h"
void make_lattice();
void over_relaxation();
void heatbath();
void initialize();
void bin();
void autocorr(double, double, double);

void main(){
  initialize();
  make_lattice();
  int i;
  int j;
  for(i=0;i<WARMUP;i++){
    for(j=0;j<RELAX_STEPS;j++)
      over_relaxation();
    heatbath();
    printf("%d\n",i);
  }
  printf("warmed up \n");
  getchar();

  for(i=0;i<MCSTEPS;i++){
    for(j=0;j<RELAX_STEPS;j++)
      over_relaxation();
    heatbath();
    bin();
    autocorr(energy,energy2,0);
    printf("%d\n",i);
  }



}
