#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>
#include "global.h"
#define LEN 8
void fft3d(double *ar,double *ai,int sign){
  //3d fft mod variables
  int stride[3];
  int lengths[3]={LX,LY,LZ};
  int sublat;
  int dir,dim1,dim2,dsize0,dsize1,dsize2;

  //1d fft variables definition
  int m,k,i,j;
  double wr,wi,tempr,tempi,termr,termi;
  int dummy;
  int sitei,sitej;
  for(dir=0; dir <3; dir++){
    switch(dir){
      case 0: 
	dsize0=LX;
	dsize1=LY;
	dsize2=LZ;
	stride[0]=1;
	stride[1]=LX;
	stride[2]=LX*LY;
	break;
      case 1: 
	dsize0=LY;
	dsize1=LZ;
	dsize2=LX;
	stride[0]=LX;
	stride[1]=LX*LY;
	stride[2]=1;
	break;
      case 2: 
	dsize0=LZ;
	dsize1=LX;
	dsize2=LY;
	stride[0]=LX*LY;
	stride[1]=1;
	stride[2]=LX;
    }

    for(sublat=0; sublat<12; sublat++){
      //start transform
      for(dim1=0;dim1<dsize1;dim1++){
	for(dim2=0;dim2<dsize2;dim2++){

	  //bit reverse
	  j=0;
	  for(i=0;i<dsize0;i++){
	    sitei=12*(i*stride[0]+dim1*stride[1] + dim2*stride[2]) + sublat;
	    if(j>i){
	      sitej=12*(j*stride[0]+dim1*stride[1] + dim2*stride[2]) + sublat;
	      tempr=ar[sitei];
	      ar[sitei]=ar[sitej];
	      ar[sitej]=tempr;

	      tempr=ai[sitei];
	      ai[sitei]=ai[sitej];
	      ai[sitej]=tempr;
	    }
	    m=dsize0/2;
	    while(m>=2 && j>=m){
	      j-=m;
	      m/=2;
	    }
	    j+=m;
	  }

	  k=1;
	  //FFT recursions

	  while(k<dsize0){
	    dummy=0;
	    wr=1.0; wi=0.0;
	    for(i=0;i<dsize0;i++){
	      j=i+k;
	      sitei=12*(i*stride[0]+dim1*stride[1] + dim2*stride[2]) + sublat;
	      sitej=12*(j*stride[0]+dim1*stride[1] + dim2*stride[2]) + sublat;
	      tempr=ar[sitei];tempi=ai[sitei];
	      termr=wr*ar[sitej]-wi*ai[sitej];
	      termi=wi*ar[sitej]+wr*ai[sitej];
	      ar[sitei]=ar[sitei]+termr;
	      ai[sitei]=ai[sitei]+termi;
	      ar[sitej]=tempr-termr;
	      ai[sitej]=tempi-termi;
	      tempr=wr;
	      wr=   wr*cos(M_PI/(1.0*k)) -wi*sin(M_PI/(1.0*k));
	      wi=tempr*sin(M_PI/(1.0*k)) +wi*cos(M_PI/(1.0*k));
	      dummy++;
	      if(dummy==k){
		i+=k;
		dummy=0;
		wr=1;
		wi=0;
		//	printf("reset\n");
	      }
	    }
	    k=k*2; 
	  } //recursion loop


	} //dim2 loop

      } //dim1 loop



    }//sublat loop
  }//direction loop(x,y,z)

}//main

void main(){
  ;

}
