#include "global.h" 
#include "assert.h" 
double genrand64_real2(void); 
void over_relaxation(){
  int i;
  int j;
  double hx,hy,hz;
  double sdoth,hdoth;
  int site;
  for(site=0;site<nsites;site++){
    i=(int)(genrand64_real2()*nsites);
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
    if(hdoth>10e-10){
    sx[i]=2*sdoth*hx/(hdoth)-sx[i];
    sy[i]=2*sdoth*hy/(hdoth)-sy[i];
    sz[i]=2*sdoth*hz/(hdoth)-sz[i];
    }
  }

}

void heatbath(){
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
    i=(int)(genrand64_real2()*nsites);
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
    h_xy=sqrt(h_xy);
    //Be Careful
    
    if(hdoth>10e-10 ){

    thetaH=acos(hz/sqrt(hdoth));
    sinthetaH=h_xy/modH;
    costhetaH=hz/modH;
    if(fabs(h_xy)>10e-10){
      cosphiH=hx/h_xy;
      sinphiH=hy/h_xy;
    }
    else{
      cosphiH=1;//does not really matter
      sinphiH=0;
    
    }
	


    phi=(double)(genrand64_real2()*2*M_PI);
    double r=(genrand64_real2());
    costheta=(1/(beta*modH))*logl(1+r*(expl(2*beta*modH)-1.0)) -1.0;
    if(isinf(costheta))
      costheta=(1/(beta*modH))*( log(r)+2*beta*modH)-1.0;
    
    double sintheta=sqrt(1-costheta*costheta);
    if(isnan(sintheta))
      sintheta=0;

    sinthetaH=sin(thetaH);
    sinphi=sin(phi);
    cosphi=cos(phi);

    tempsx=costheta*sinthetaH+sintheta*costhetaH+cosphi;
    tempsy=sintheta*sinphi;

    sx[i]=tempsx*cosphiH-tempsy*sinphiH;
    sy[i]=tempsx*sinphiH+tempsy*cosphiH;
    sz[i]=costheta*costhetaH-sintheta*sinthetaH*cosphi;
    }


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

