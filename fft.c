#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>
#define LEN 512
void fft(double *ar,double *ai,int n,int sign){
  int i,j,k;
  double tempr,tempi;
  int iprime=0;
  int m;
  FILE *fp;
  for(i=0;i<n;i++){
    printf("%d %d\n",i,iprime);
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
  unsigned int index=(unsigned int)n;
  int log2n=0;
  k=1; //k is 2**level
  while(index>>=1) log2n++;
  double wr;
  double wi;
  int dummy;

  while(k<n){
    dummy=0;
    wr=1.0; wi=0.0;
    for(i=0;i<n;i++){
      j=i+k;
      printf("lv %d i %d j %d\n",k,i,j);
      tempr=ar[i];tempi=ai[i];
      //a[i]=a[i]+w*a[j]
      ar[i]=ar[i]+wr*ar[j]-wi*ai[j];
      ai[i]=ai[i]+wi*ar[j]+wr*ai[j];
      //a[j]=a[i]-w*a[j]
      ar[j]=tempr-wr*ar[j]+wi*ai[j];
      ai[j]=tempi-wi*ar[j]-wr*ai[j];
      tempr=wr;
      wr=   wr*cos(M_PI/(1.0*k)) -wi*sin(M_PI/(1.0*k));
      wi=tempr*sin(M_PI/(1.0*k)) +wi*cos(M_PI/(1.0*k));

      dummy++;
      if(dummy==k){
	i+=k;
	dummy=0;
	wr=1;
	wi=0;
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
    array[i]=cos(2*M_PI*i*100/(1.0*LEN));
   arrayi[i]=sin(2*M_PI*i*100/(1.0*LEN));
  }
  fft(array,arrayi,LEN,1);
  fp=fopen("testdata.dat","w");
  for(i=0;i<LEN;i++)
    fprintf(fp,"%f %f %f\n",2*M_PI*(1.0*i)/((double)LEN),array[i],arrayi[i]);
  fclose(fp);
}
