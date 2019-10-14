#include "global.h"
double min(double a,double b){
	return a*(a<b)+b*(b<=a);
}
void measure_sq(){
  sq_measure++;
  int i,j;
  int cell;
  int idx;// coordinates of lower triangular matrix
  int x,y,z;
  int mx,my,mz;
  int idx1,idx2;// coordinates of q and minus q
  int mcell;// unit cell at -q
  int strides[12]={0,12,23,33,42,50,57,63,68,72,75,77};
  for(i=0;i<12;i++){
     for(j=0;j<=i;j++){
       idx=strides[i]+j;
	 for(x=0;x<lx;x++){
	   for(y=0;y<ly;y++){
	    for(z=0;z<lz;z++) { 
	      cell=x+lx*y+lx*ly*z;
	      printf("lprep: %d %d %d \n",lx,ly,lz);
	      printf("prep: %d %d %d \n",x,y,z);
	      mcell=(lx-x)%lx+((ly-y)%ly)*lx+((lz-z)%lz)*lx*ly;
	      idx1=12*cell+i;
	      idx2=12*mcell+j;
	      if(idx2>=nsites){
		printf("%d %d %d %d\n",i,j,cell,mcell);
		printf("%d %d %d \n",x,y,z);
		printf("%d %d %d \n",lx,ly,lz);
		printf("%d %d %d \n",(lx-x)%lx,(ly-y)%ly,(lz-z)%lz);
	      }
	      //assert(idx1<nsites);
	      //assert(idx2<nsites);
	      assert(idx<78*ncells);

              sqasqbre[idx*NCELLS+cell]+= (sxqr[idx1]*sxqr[idx2]-sxqi[idx1]*sxqi[idx2]);
              sqasqbre[idx*NCELLS+cell]+= (syqr[idx1]*syqr[idx2]-syqi[idx1]*syqi[idx2]);
              sqasqbre[idx*NCELLS+cell]+= (szqr[idx1]*szqr[idx2]-szqi[idx1]*szqi[idx2]);

              sqasqbim[idx*NCELLS+cell]+= (sxqr[idx1]*sxqi[idx2]+sxqi[idx1]*sxqr[idx2]);
              sqasqbim[idx*NCELLS+cell]+= (syqr[idx1]*syqi[idx2]+syqi[idx1]*syqr[idx2]);
              sqasqbim[idx*NCELLS+cell]+= (szqr[idx1]*szqi[idx2]+szqi[idx1]*szqr[idx2]);
	    }//z
	   }//y
	 }//x
     }//j
  }//i

  if(sq_measure%(sq_samples/5)==0){
    FILE *fp;
    for(i=0;i<12;i++){
       for(j=0;j<=i;j++){
         idx=strides[i]+j;
	 sprintf(sfacname,"%si%dj%d.dat",sfacnamepref,i,j);
	 fp=fopen(sfacname,"w");
	 for(cell=0;cell<NCELLS;cell++)
	   fprintf(fp,"%f %f\n",sqasqbre[idx*NCELLS+cell],sqasqbim[idx*NCELLS+cell]);

       }
       fclose(fp);
    }

  }

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
  static double o1list[TAUMAX],o2list[TAUMAX],o3list[TAUMAX];
  static double ao1[TAUMAX],ao2[TAUMAX],ao3[TAUMAX];  
  static double aeo1[TAUMAX],aeo2[TAUMAX],aeo3[TAUMAX];
  static double eo1,eo2,eo3;
  static double avo1[TAUMAX],avo2[TAUMAX],avo3[TAUMAX];
  static double del_avo1[TAUMAX],del_avo2[TAUMAX],del_avo3[TAUMAX];
  static double autocorr1[TAUMAX],autocorr2[TAUMAX],autocorr3[TAUMAX],autoerr1[TAUMAX],autoerr2[TAUMAX],autoerr3[TAUMAX];
  static int nsamples[TAUMAX];
  int i;


  nautocorr++;

  int pos=nautocorr%(taumax);

  o1list[pos]=obs1;
  o2list[pos]=obs2;
  o3list[pos]=obs3;
  
  for(i=0;i<min(nautocorr,taumax);i++){

    int pos2= (taumax+pos-i)%taumax;

    avo1[i]+=o1list[pos];
    avo2[i]+=o2list[pos];
    avo3[i]+=o3list[pos];

    del_avo1[i]+=o1list[pos2];
    del_avo2[i]+=o2list[pos2];
    del_avo3[i]+=o3list[pos2];

    ao1[i]+=o1list[pos]*o1list[pos2];
    ao2[i]+=o2list[pos]*o2list[pos2];
    ao3[i]+=o3list[pos]*o3list[pos2];
    aeo1[i]+=o1list[pos]*o1list[pos2]*o1list[pos]*o1list[pos2];
    aeo2[i]+=o2list[pos]*o2list[pos2]*o2list[pos]*o2list[pos2];
    aeo3[i]+=o3list[pos]*o3list[pos2]*o3list[pos]*o3list[pos2];
    nsamples[i]++;
  }



  if(nautocorr%(MCSTEPS/5)==0){
    FILE *fp;
    fp=fopen(autofname,"w");
    for(i=0;i<taumax;i++){ 
      autocorr1[i]=(ao1[i]/(1.0*nsamples[i]))-avo1[i]*del_avo1[i]/(1.0*(nsamples[i]*nsamples[i]));
      autocorr2[i]=(ao2[i]/(1.0*nsamples[i]))-avo2[i]*del_avo2[i]/(1.0*(nsamples[i]*nsamples[i]));
      autocorr3[i]=(ao3[i]/(1.0*nsamples[i]))-avo3[i]*del_avo3[i]/(1.0*(nsamples[i]*nsamples[i]));


      autoerr1[i]=(aeo1[i]/(1.0*nsamples[i])-ao1[i]*ao1[i]/(1.0*nsamples[i]*nsamples[i]))/(1.0*autocorr1[0]*autocorr1[0]);
      autoerr2[i]=(aeo2[i]/(1.0*nsamples[i])-ao2[i]*ao2[i]/(1.0*nsamples[i]*nsamples[i]))/(1.0*autocorr2[0]*autocorr2[0]);
      autoerr3[i]=(aeo3[i]/(1.0*nsamples[i])-ao3[i]*ao3[i]/(1.0*nsamples[i]*nsamples[i]))/(1.0*autocorr3[0]*autocorr3[0]);
    }
    for(i=1;i<taumax;i++){ 
      autocorr1[i]=autocorr1[i]/autocorr1[0];
      autocorr2[i]=autocorr2[i]/autocorr2[0];
      autocorr3[i]=autocorr3[i]/autocorr3[0];
    }
    autocorr1[0]=autocorr2[0]=autocorr3[0]=1;
    for(i=0;i<taumax;i++) 
      fprintf(fp,"%d %.16f %.16f %.16f %.16f %.16f %.16f\n",i,autocorr1[i], sqrt(autoerr1[i]/(1.0*nsamples[i])),autocorr2[i],sqrt(autoerr2[i]/(1.0*nsamples[i])), autocorr3[i], sqrt(autoerr3[i]/(1.0*nsamples[i])));
    fclose(fp);

  }
}
