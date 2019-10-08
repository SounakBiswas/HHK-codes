#include "global.h"
void make_lattice();
void over_relaxation();
void heatbath();
void initialize();

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

  for(i=0;i<MCSTEPS;i++){
    for(j=0;j<RELAX_STEPS;j++)
      over_relaxation();
    heatbath();
  }



}
