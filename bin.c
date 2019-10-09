#include "global.h"
double min(double a,double b){
	return a*(a<b)+b*(b<=a);
}
void bin(){
  if(nmeasure%binsize==0){
    energy_bin=0;
    energy2_bin=0;

  }
  energy_bin+=energy;
  energy2_bin+=energy2;

  if(nmeasure%binsize==(binsize-1)){
    FILE *binfp;
    binfp=fopen(binfname,"a");
    fprintf(binfp, "%d %lld %.16f %.16f\n",binno,nmeasure,energy_bin/((double)binsize), energy2_bin/((double)binsize));
    binno++;
    fclose(binfp);

  }
  nmeasure++;

}

void autocorr(double obs1, double obs2,double obs3){
  static double o1,o2,o3,o1err,o2err,o3err;
  static double o1list[TAUMAX],o2list[TAUMAX],o3list[TAUMAX];
  static double ao1[TAUMAX],ao2[TAUMAX],ao3[TAUMAX];
  static double aeo1[TAUMAX],aeo2[TAUMAX],aeo3[TAUMAX];
  static double eo1,eo2,eo3;
  static double avo1,avo2,avo3;
  static double avo1sq,avo2sq,avo3sq;
  double autocorr1[TAUMAX],autocorr2[TAUMAX],autocorr3[TAUMAX],autoerr1[TAUMAX],autoerr2[TAUMAX],autoerr3[TAUMAX];
  static int nsamples[TAUMAX];
  int i;
  int posi;
  int pos=nautocorr%(taumax);
  o1list[pos]=obs1;
  o2list[pos]=obs2;
  o3list[pos]=obs3;
  avo1+=obs1;
  avo2+=obs2;
  avo3+=obs3;

  avo1sq+=obs1*obs1;
  avo2sq+=obs2*obs2;
  avo3sq+=obs3*obs3;


  nautocorr++;

  for(i=0;i<min(nautocorr,taumax);i++){

    ao1[i]+=o1list[pos]*o1list[(taumax+pos-i)%taumax];
    ao2[i]+=o2list[pos]*o2list[(taumax+pos-i)%taumax];
    ao3[i]+=o3list[pos]*o3list[(taumax+pos-i)%taumax];
    aeo1[i]+=o1list[pos]*o1list[(taumax+pos-i)%taumax]*o1list[pos]*o1list[(taumax+pos-i)%taumax];
    aeo2[i]+=o2list[pos]*o2list[(taumax+pos-i)%taumax]*o2list[pos]*o2list[(taumax+pos-i)%taumax];
    aeo3[i]+=o3list[pos]*o3list[(taumax+pos-i)%taumax]*o3list[pos]*o3list[(taumax+pos-i)%taumax];
    nsamples[i]++;
  }



  if(nautocorr%(MCSTEPS/5)==0){
    FILE *fp;
    fp=fopen(autofname,"w");
    eo1=avo1sq/(1.0*nautocorr)-avo1*avo1/(1.0*nautocorr*nautocorr);
    eo2=avo2sq/(1.0*nautocorr)-avo2*avo2/(1.0*nautocorr*nautocorr);
    eo3=avo3sq/(1.0*nautocorr)-avo3*avo3/(1.0*nautocorr*nautocorr);
    for(i=0;i<taumax;i++){ 
      autocorr1[i]=(ao1[i]/(1.0*nsamples[i]))-avo1*avo1/(1.0*nautocorr*nautocorr);
      autocorr1[i]=autocorr1[i]/(eo1);

      autocorr2[i]=(ao2[i]/(1.0*nsamples[i]))-avo2*avo2/(1.0*nautocorr*nautocorr);
      autocorr2[i]=autocorr2[i]/eo2;

      autocorr3[i]=(ao3[i]/(1.0*nsamples[i]))-avo3*avo3/(1.0*nautocorr*nautocorr);
      autocorr3[i]=autocorr3[i]/eo3;

      autoerr1[i]=(aeo1[i]/(1.0*nsamples[i])-ao1[i]*ao1[i]/(1.0*nsamples[i]*nsamples[i]))/(1.0*eo1);
      autoerr2[i]=(aeo2[i]/(1.0*nsamples[i])-ao2[i]*ao2[i]/(1.0*nsamples[i]*nsamples[i]))/(1.0*eo2);
      autoerr3[i]=(aeo3[i]/(1.0*nsamples[i])-ao3[i]*ao3[i]/(1.0*nsamples[i]*nsamples[i]))/(1.0*eo3);

      fprintf(fp,"%d %.16f %.16f %.16f %.16f %.16f %.16f\n",i,autocorr1[i], sqrt(autoerr1[i]/(1.0*nsamples[i]))+sqrt(eo1/nautocorr),autocorr2[i],sqrt(autoerr2[i]/(1.0*nsamples[i]))+sqrt(eo2/nautocorr), autocorr3[i], sqrt(autoerr3[i]/(1.0*nsamples[i]))+sqrt(eo3/nautocorr));
    }
    fclose(fp);

  }
}
