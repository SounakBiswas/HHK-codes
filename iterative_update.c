#include "global.h" 
double getenergy(){
  int i,j;
  double e=0;
  double hx,hy,hz;
  for(i=0;i<nsites;i++){ 
    hx=hy=hz=0;
    for(j=0; j<4; j++){
      hx+= -sx[neigh[6*i+j]];
      hy+= -sy[neigh[6*i+j]];
      hz+= -sz[neigh[6*i+j]];
    }

    for(j=4; j<6; j++){
      hx+= -j2*sx[neigh[6*i+j]];
      hy+= -j2*sy[neigh[6*i+j]];
      hz+= -j2*sz[neigh[6*i+j]];
    }
    e+= (sx[i]*hx+sy[i]*hy+sz[i]*hz)/2.0;
  }
  return e;


}
void over_relaxation(){
  int i;
  int j;
  double hx,hy,hz;
  double sdoth,hdoth;
  int site;
  FILE *fp;
  for(site=0;site<nsites;site++){
    i=(int)((1.0*rand())/((double)(RAND_MAX+1.0))*nsites);
    i=site;
    hx=hy=hz=0;
    for(j=0; j<4; j++){
      hx+= -sx[neigh[6*i+j]];
      hy+= -sy[neigh[6*i+j]];
      hz+= -sz[neigh[6*i+j]];
    }

    for(j=4; j<6; j++){
      hx+= -j2*sx[neigh[6*i+j]];
      hy+= -j2*sy[neigh[6*i+j]];
      hz+= -j2*sz[neigh[6*i+j]];
    }
    //printf("%d %f %f %f\n",i,hx,hy,hz);
    //getchar();
    hdoth=hx*hx+hy*hy+hz*hz;
    sdoth=sx[i]*hx+sy[i]*hy+sz[i]*hz;
    //BE CAREFUL
    //if(hdoth>10e-10){
    sx[i]=2*sdoth*hx/hdoth-sx[i];
    sy[i]=2*sdoth*hy/hdoth-sy[i];
    sz[i]=2*sdoth*hz/hdoth-sz[i];


   // }
  }

}
void iterative_gs(){
  int i,j;
  double hx,hy,hz;
  double sdoth,hdoth;
  int site;
  double phiH,thetaH;
  double costheta,phi,modH;
  double costhetaH,sinthetaH, cosphiH, sinphiH;
  double sintheta,sinphi,cosphi;
  double tempsx,tempsy;
  double h_xy;
  for(site=0;site<nsites;site++){
    i=(int)((1.0*rand())/((double)(RAND_MAX+1.0))*nsites);
    //printf("%d \n",i);
    //assert(i<nsites);
    //i=site;
    hx=hy=hz=0;
    for(j=0; j<4; j++){
      hx+= -sx[neigh[6*i+j]];
      hy+= -sy[neigh[6*i+j]];
      hz+= -sz[neigh[6*i+j]];
    }

    for(j=4; j<6; j++){
      hx+= -j2*sx[neigh[6*i+j]];
      hy+= -j2*sy[neigh[6*i+j]];
      hz+= -j2*sz[neigh[6*i+j]];
    }
    h_xy=hx*hx+hy*hy;
    hdoth=h_xy+hz*hz;
    modH=sqrt(hdoth);
    sx[i]=hx/modH;
    sy[i]=hy/modH;
    sz[i]=hz/modH;
    //h_xy=sqrt(h_xy);
    ////Be Careful
    //
    //if(hdoth>10e-10 && h_xy>10e-10){

    //thetaH=acos(hz/sqrt(hdoth));
    //sinthetaH=h_xy/modH;
    //costhetaH=hz/modH;
    //cosphiH=hx/h_xy;
    //sinphiH=hy/h_xy;


    //phi=(double)(rand()/(1.0*RAND_MAX))*2*M_PI;
    //double r=(double)(rand()/(1.0*RAND_MAX))*1;
    ////costheta=(1/(beta*modH))*log(1+r*(exp(2*beta*modH)-1.0)) -1.0;
    //costheta=1;
    //double sintheta=0;
    //sinthetaH=sin(thetaH);
    //sinphi=sin(phi);
    //cosphi=cos(phi);

    //tempsx=costheta*sinthetaH+sintheta*costhetaH+cosphi;
    //tempsy=sintheta*sinphi;

    //sx[i]=tempsx*cosphiH-tempsy*sinphiH;
    //sy[i]=tempsx*sinphiH+tempsy*cosphiH;
    //sz[i]=costheta*costhetaH-sintheta*sinthetaH*cosphi;
    //}


  }
  energy=0;
  for(i=0;i<nsites;i++){ 
    hx=hy=hz=0;
    for(j=0; j<4; j++){
      hx+= -sx[neigh[6*i+j]];
      hy+= -sy[neigh[6*i+j]];
      hz+= -sz[neigh[6*i+j]];
    }

    for(j=4; j<6; j++){
      hx+= -j2*sx[neigh[6*i+j]];
      hy+= -j2*sy[neigh[6*i+j]];
      hz+= -j2*sz[neigh[6*i+j]];
    }
    energy+= (sx[i]*hx+sy[i]*hy+sz[i]*hz)/2.0;
  }
  energy2=energy*energy;

}
