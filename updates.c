#include "global.h" 
#include "assert.h" 
void over_relaxation(){
  int i;
  double hx,hy,hz;
  double sdoth,hdoth;
  int site;
  for(site=0;site<nsites;site++){
    i=(int)((1.0*rand())/((double)(RAND_MAX+1.0))*nsites);
    hx=hy=hz=0;
    for(j=0; j<4; j++){
      hx+= sx[neigh[6*i+j]];
      hy+= sy[neigh[6*i+j]];
      hz+= sz[neigh[6*i+j]];
    }

    for(j=5; j<6; j++){
      hx+= j2*sx[neigh[6*i+j]];
      hy+= j2*sy[neigh[6*i+j]];
      hz+= j2*sz[neigh[6*i+j]];
    }
    hdoth=hx*hx+hy*hy+hz*hz;
    sdoth=sx[i]*hx+sy[i]*hy+sz[i]*hz;
    assert(hdoth>10e-8);
    sx[i]=2*sdoth*hx/(hdoth)-sx[i];
    sy[i]=2*sdoth*hy/(hdoth)-sy[i];
    sz[i]=2*sdoth*hz/(hdoth)-sz[i];
  }

}

void heatbath(){
  int i;
  double hx,hy,hz;
  double sdoth,hdoth;
  int site;
  for(site=0;site<nsites;site++){
    i=(int)((1.0*rand())/((double)(RAND_MAX+1.0))*nsites);
    hx=hy=hz=0;
    for(j=0; j<4; j++){
      hx+= sx[neigh[6*i+j]];
      hy+= sy[neigh[6*i+j]];
      hz+= sz[neigh[6*i+j]];
    }

    for(j=5; j<6; j++){
      hx+= j2*sx[neigh[6*i+j]];
      hy+= j2*sy[neigh[6*i+j]];
      hz+= j2*sz[neigh[6*i+j]];
    }

  }

}
