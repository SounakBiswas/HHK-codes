#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>
#define LEN 64
void fft(double *ar,double *ai,int n,int sign){
  int i,j,k;
  double tempr,tempi;
  int iprime=0;
  int m;
  FILE *fp;
  for(i=0;i<n;i++){
    if(iprime>i){
      tempr=ar[i];
      ar[i]=ar[iprime];
      ar[iprime]=tempr;

      tempr=ai[i];
      ai[i]=ai[iprime];
      ai[iprime]=tempr;
    }
    m=n/2;
    while(m>=2 && iprime>=m){
      iprime-=m;
      m/=2;
    }
    iprime+=m;
  }
  k=1; //k is 2**level
  double wr;
  double wi;
  int dummy;
  double termr,termi;

  while(k<n){
    dummy=0;
    wr=1.0; wi=0.0;
    for(i=0;i<n;i++){
      j=i+k;
      tempr=ar[i];tempi=ai[i];
      termr=wr*ar[j]-wi*ai[j];
      termi=wi*ar[j]+wr*ai[j];
      ar[i]=ar[i]+termr;
      ai[i]=ai[i]+termi;
      //a[j]=a[i]-w*a[j]
      ar[j]=tempr-termr;
      ai[j]=tempi-termi;
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
  }

	
}
void main(){
  double array[LEN],arrayi[LEN];
  int i;
  FILE *fp;
  for(i=0;i<LEN;i++){
    array[i]=cos(2*M_PI*i*10.0/(1.0*LEN));
   arrayi[i]=-sin(2*M_PI*i*10.0/(1.0*LEN));
  }
  fft(array,arrayi,LEN,1);
  fp=fopen("testdata.dat","w");
  for(i=0;i<LEN;i++)
    fprintf(fp,"%f %f %f\n",2*M_PI*(1.0*i)/((double)LEN),array[i],arrayi[i]);
  fclose(fp);
}
