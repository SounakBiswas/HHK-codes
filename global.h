#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#define LX 8
#define LY LX
#define LZ LX
#define NCELLS LX*LY*LZ 
#define NSITES 12*NCELLS
#define J2 1.0000
#define TEMP 0.2000
#define WARMUP 50000
#define MCSTEPS 500000
#define RELAX_STEPS 10
#define BINSIZE 200
#define SFACM 1
#define SEED 11
#define TAUMAX 100
#define FILENUM 10


int neigh[NSITES*6];
int sq_samples;
int sq_measure; // the current number of times sq have been measured.
int nsites,ncells;
int nautocorr;
int taumax;
double temp,beta;
double j2;
double sx[NSITES], sy[NSITES], sz[NSITES];
double sxqr[NSITES], syqr[NSITES], szqr[NSITES];
double sxqi[NSITES], syqi[NSITES], szqi[NSITES];
double energy;
double energy2;
double energy_bin, energy2_bin;
int lx,ly,lz;
int binsize;
char binfname[200],autofname[200],sfacnamepref[200],sfacname[200];
double sqasqbre[78*NCELLS];
double sqasqbim[78*NCELLS];
long long int nmeasure;
int binno;

