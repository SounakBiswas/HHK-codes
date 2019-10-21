#include "global.h"
//#define next_nbr(x,y,z,i,j,k) (LX+x+i)%LX+ ((LY+y+j)%LY)*LX +((LZ+z+k)%LZ)*LX*LY
#define next_nbr(x,y,z,i,j,k) ((LX+x+i)%LX+ ((LY+y+j)%LY)*LX +((LZ+z+k)%LZ)*LX*LY)
void init_genrand64(unsigned long long);
double genrand64_real2(void);

void initialize(){
  int i,j;
  ncells=NCELLS;
  nsites=NSITES;
  temp=TEMP;
  beta=1.0/temp;
  lx=LX; 
  ly=LY; 
  lz=LZ;
  printf("%d %d %d\n",lx,ly,lz);
  j2=J2;
  int x,y,z;
  double theta,phi;
  sq_samples=(int)(MCSTEPS/(1.0*SFACM));
  sq_measure=0;
  taumax=TAUMAX;
  binsize=BINSIZE;
  binno=0;
  nmeasure=0;
  nautocorr=0;
  int seed=SEED;
  init_genrand64(seed);
  for(i=0;i<nsites;i++){
	  theta=M_PI*genrand64_real2();
          phi=2*M_PI*genrand64_real2();
	  sz[i]=cos(theta);
	  sx[i]=sin(theta)*cos(phi);
	  sy[i]=sin(theta)*sin(phi);
  }
  sprintf(binfname,"./outfiles/bin_L%dT%.4fJ2%.4fseed%d.dat",lx,temp,j2,FILENUM);
  sprintf(autofname,"./outfiles/auto_L%dT%.4fJ2%.4fseed%d.dat",lx,temp,j2,FILENUM);
  sprintf(sfacnamepref,"./sfacs/sfac_L%dT%.4fJ2%.4fseed%d",lx,temp,j2,FILENUM);
  double hx,hy,hz;
  
//  energy=0;
//  FILE *fp=fopen("data2.dat","w");
//  for(i=0;i<nsites;i++){ 
//    hx=hy=hz=0;
//    for(j=0; j<4; j++){
//      hx+= -sx[neigh[6*i+j]];
//      hy+= -sy[neigh[6*i+j]];
//      hz+= -sz[neigh[6*i+j]];
//    }
//
//    for(j=4; j<6; j++){
//      hx+= -j2*sx[neigh[6*i+j]];
//      hy+= -j2*sy[neigh[6*i+j]];
//      hz+= -j2*sz[neigh[6*i+j]];
//    }
//  fprintf(fp,"%d %f %f %f %f %f %f\n",i,sx[i],sy[i],sz[i],hx,hy,hz);
//
//    energy+= (sx[i]*hx+sy[i]*hy+sz[i]*hz)/2.0;
//  }
//  fclose(fp);
//  printf("starting energy=%f\n",energy);
//  energy2=energy*energy;

  for(i=0;i<78*NCELLS;i++)
    sqasqbre[i]=sqasqbim[i]=0.0;

}
void make_lattice(){
  int i;
  int x,y,z;
  for(i=0;i<ncells;i++){
    z=(i/(lx*ly));
    y=(i%(lx*ly))/lx;
    x=(i%(lx*ly))%lx;
    neigh[(i*12+0)*6+0 ] = 12*next_nbr(x,y,z,-1,0,0) + 9;
    neigh[(i*12+0)*6+1 ] = 12*next_nbr(x,y,z,0,0,0) + 1;
    neigh[(i*12+0)*6+2 ] = 12*next_nbr(x,y,z,0,-1,0) + 5;
    neigh[(i*12+0)*6+3 ] = 12*next_nbr(x,y,z,0,0,-1) + 3;
    neigh[(i*12+0)*6+4 ] = 12*next_nbr(x,y,z,-1,-1,0) + 8;
    neigh[(i*12+0)*6+5 ] = 12*next_nbr(x,y,z,-1,-1,-1) + 10;

    neigh[(i*12+1)*6+0 ] = 12*next_nbr(x,y,z,0,0,0) + 0;
    neigh[(i*12+1)*6+1 ] = 12*next_nbr(x,y,z,-1,0,0) + 9;
    neigh[(i*12+1)*6+2 ] = 12*next_nbr(x,y,z,0,0,0) + 2;
    neigh[(i*12+1)*6+3 ] = 12*next_nbr(x,y,z,0,0,0) + 4;
    neigh[(i*12+1)*6+4 ] = 12*next_nbr(x,y,z,0,0,0) + 3;
    neigh[(i*12+1)*6+5 ] = 12*next_nbr(x,y,z,-1,-1,0) + 11;

    neigh[(i*12+2)*6+0 ] = 12*next_nbr(x,y,z,0,0,0) + 1;
    neigh[(i*12+2)*6+1 ] = 12*next_nbr(x,y,z,0,0,0) + 4;
    neigh[(i*12+2)*6+2 ] = 12*next_nbr(x,y,z,-1,0,0) + 10;
    neigh[(i*12+2)*6+3 ] = 12*next_nbr(x,y,z,-1,0,0) + 11;
    neigh[(i*12+2)*6+4 ] = 12*next_nbr(x,y,z,0,0,0) + 7;
    neigh[(i*12+2)*6+5 ] = 12*next_nbr(x,y,z,0,0,1) + 5;


    neigh[(i*12+3)*6+0 ] = 12*next_nbr(x,y,z,0,0,1) + 0;
    neigh[(i*12+3)*6+1 ] = 12*next_nbr(x,y,z,0,-1,1) + 5;
    neigh[(i*12+3)*6+2 ] = 12*next_nbr(x,y,z,0,0,0) + 6;
    neigh[(i*12+3)*6+3 ] = 12*next_nbr(x,y,z,0,0,0) + 7;
    neigh[(i*12+3)*6+4 ] = 12*next_nbr(x,y,z,0,0,0) + 1;
    neigh[(i*12+3)*6+5 ] = 12*next_nbr(x,y,z,-1,-1,0) + 11;


    neigh[(i*12+4)*6+0 ] = 12*next_nbr(x,y,z,0,0,0) + 1;
    neigh[(i*12+4)*6+1 ] = 12*next_nbr(x,y,z,0,0,0) + 2;
    neigh[(i*12+4)*6+2 ] = 12*next_nbr(x,y,z,0,0,0) + 5;
    neigh[(i*12+4)*6+3 ] = 12*next_nbr(x,y,z,0,0,0) + 8;
    neigh[(i*12+4)*6+4 ] = 12*next_nbr(x,y,z,0,0,0) + 6;
    neigh[(i*12+4)*6+5 ] = 12*next_nbr(x,y,z,0,0,0) + 9;


    neigh[(i*12+5)*6+0 ] = 12*next_nbr(x,y,z,0,1,0) + 0;
    neigh[(i*12+5)*6+1 ] = 12*next_nbr(x,y,z,0,1,-1) + 3;
    neigh[(i*12+5)*6+2 ] = 12*next_nbr(x,y,z,0,0,0) + 4;
    neigh[(i*12+5)*6+3 ] = 12*next_nbr(x,y,z,0,0,0) + 8;
    neigh[(i*12+5)*6+4 ] = 12*next_nbr(x,y,z,0,0,-1) + 2;
    neigh[(i*12+5)*6+5 ] = 12*next_nbr(x,y,z,0,0,-1) + 7;


    neigh[(i*12+6)*6+0 ] = 12*next_nbr(x,y,z,0,0,0) + 3;
    neigh[(i*12+6)*6+1 ] = 12*next_nbr(x,y,z,0,0,0) + 7;
    neigh[(i*12+6)*6+2 ] = 12*next_nbr(x,y,z,0,-1,0) + 8;
    neigh[(i*12+6)*6+3 ] = 12*next_nbr(x,y,z,0,-1,0) + 11;
    neigh[(i*12+6)*6+4 ] = 12*next_nbr(x,y,z,0,0,0) + 4;
    neigh[(i*12+6)*6+5 ] = 12*next_nbr(x,y,z,0,0,0) + 9;


    neigh[(i*12+7)*6+0 ] = 12*next_nbr(x,y,z,0,0,0) + 3;
    neigh[(i*12+7)*6+1 ] = 12*next_nbr(x,y,z,0,0,0) + 6;
    neigh[(i*12+7)*6+2 ] = 12*next_nbr(x,y,z,0,0,1) + 9;
    neigh[(i*12+7)*6+3 ] = 12*next_nbr(x,y,z,0,0,0) + 10;
    neigh[(i*12+7)*6+4 ] = 12*next_nbr(x,y,z,0,0,0) + 2;
    neigh[(i*12+7)*6+5 ] = 12*next_nbr(x,y,z,0,0,1) + 5;

    neigh[(i*12+8)*6+0 ] = 12*next_nbr(x,y,z,0,0,0) + 4;
    neigh[(i*12+8)*6+1 ] = 12*next_nbr(x,y,z,0,0,0) + 5;
    neigh[(i*12+8)*6+2 ] = 12*next_nbr(x,y,z,0,1,0) + 6;
    neigh[(i*12+8)*6+3 ] = 12*next_nbr(x,y,z,0,0,0) + 11;
    neigh[(i*12+8)*6+4 ] = 12*next_nbr(x,y,z,1,1,0) + 0;
    neigh[(i*12+8)*6+5 ] = 12*next_nbr(x,y,z,0,0,-1) + 10;


    neigh[(i*12+9)*6+0 ] = 12*next_nbr(x,y,z,1,0,0) + 0;
    neigh[(i*12+9)*6+1 ] = 12*next_nbr(x,y,z,1,0,0) + 1;
    neigh[(i*12+9)*6+2 ] = 12*next_nbr(x,y,z,0,0,-1) + 7;
    neigh[(i*12+9)*6+3 ] = 12*next_nbr(x,y,z,0,0,-1) + 10;
    neigh[(i*12+9)*6+4 ] = 12*next_nbr(x,y,z,0,0,0) + 4;
    neigh[(i*12+9)*6+5 ] = 12*next_nbr(x,y,z,0,0,0) + 6;


    neigh[(i*12+10)*6+0 ] = 12*next_nbr(x,y,z,1,0,0) + 2;
    neigh[(i*12+10)*6+1 ] = 12*next_nbr(x,y,z,0,0,0) + 11;
    neigh[(i*12+10)*6+2 ] = 12*next_nbr(x,y,z,0,0,0) + 7;
    neigh[(i*12+10)*6+3 ] = 12*next_nbr(x,y,z,0,0,1) + 9;
    neigh[(i*12+10)*6+4 ] = 12*next_nbr(x,y,z,1,1,1) + 0;
    neigh[(i*12+10)*6+5 ] = 12*next_nbr(x,y,z,0,0,1) + 8;

    neigh[(i*12+11)*6+0 ] = 12*next_nbr(x,y,z,0,0,0) + 8;
    neigh[(i*12+11)*6+1 ] = 12*next_nbr(x,y,z,0,1,0) + 6;
    neigh[(i*12+11)*6+2 ] = 12*next_nbr(x,y,z,0,0,0) + 10;
    neigh[(i*12+11)*6+3 ] = 12*next_nbr(x,y,z,1,0,0) + 2;
    neigh[(i*12+11)*6+4 ] = 12*next_nbr(x,y,z,1,1,0) + 1;
    neigh[(i*12+11)*6+5 ] = 12*next_nbr(x,y,z,1,1,0) + 3;

  }

}
