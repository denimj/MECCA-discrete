// *************************
// Program by Luke J. Harmon
// July 21 2011
// Updated by Daniel S. Caetano
// January 29 2013
// *************************

#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
//#define PI 3.141592654 //This does not need to be defined here, for I included the R.h library.

// definitions for ran_ecu()

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52744
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)

// function prototypes
void pMk2(double *pij, double q12, double q21, double t);
double ran_ecu(long *seed);
float bnldev(float pp, int n, long *seed);
float gammln(float xx);
void noTreeSim(int *ntimes, double *times, double *q12, double *q21, int *startState, long *seed, int *nState1, int *nState2, int *vezes, int *loops);

// end function prototypes

// functions

/*The function 'pMk2' is calculating values for 'pij' based on the values for q12 and q21 */
void pMk2(double *pij, double q12, double q21, double t) {
	double x;
	
	x = exp(-(q12+q21)*t);
	
	if( (q12==0) & (q21==0) )
    {
       	pij[1]= 1.0;
        pij[2]= 0.0;
        pij[3]= 0.0;
        pij[4]= 1.0; 
    } else {
        pij[1]= (x*q12 + q21) / (q12 + q21);
        pij[2]= ((1 - x)*q21) / (q12 + q21);
        pij[3]= ((1 - x)*q12) / (q12 + q21);
        pij[4]= (x*q21 + q12) / (q12 + q21);
    }
}

/**************************************************************************************************/
/* Ran_ecu	*/
/* Random number generator of L'Ecuyer with Bays-Durham shuffle and added safeguards from nr in c */
/* Returns good pseudorandoms in range 0-1 */
/* Seed should be positive long int */
/* Can use time(NULL) as initial seed */
/* End ran_ecu */
double ran_ecu(long *seed)
{
	int j;
	long k;
	static long seed2=123456789; /* Default seed if you forget to initialize */
	static long iy=0;
	static long iv[NTAB];
	double result;
	if(*seed <= 0) /* Initialize */
	{
		if ((*seed)==0) *seed = 1;
		seed2=(*seed);
		for(j=NTAB+7; j>=0; j--)
		{
			k=(*seed)/IQ1;
			*seed=IA1*(*seed-k*IQ1)-k*IR1;
			if (*seed<0) *seed += IM1;
			if (j < NTAB) iv[j]=*seed;
		}
		iy=iv[0];
	}
	
	k=(*seed)/IQ1;
	*seed=IA1*(*seed-k*IQ1)-k*IR1;
	if (*seed<0) *seed += IM1;
	k=seed2/IQ2;
	seed2=IA2*(seed2-k*IQ2)-k*IR2;
	if (seed2<0) seed2 += IM2;
	
	j=iy/NDIV;

	iy=iv[j]-seed2;

	iv[j]=*seed;

	if (iy<1) iy += IMM1;
    

	result=AM*iy;
	
	return(result);
}

/**************************************************************************************************/

/* Modified from Numerical Recipes in C */
 
float bnldev(float pp, int n, long *seed) {
	float gammln(float xx);
	int j;
	static int nold=(-1);
	float am,em,g,angle,p,bnl,sq,t,y;
	static float pold=(-1.0),pc,plog,pclog,en,oldg;
	double rr;
	float gam1, gam2;
	p=(pp <= 0.5 ? pp : 1.0-pp);

	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=1;j<=n;j++) {
			rr=ran_ecu(seed);
			if (rr < p) ++bnl;
		}
	} else if (am < 1.0) {

		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran_ecu(seed);
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {

		if (n != nold) { 
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) { 
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc); 
		
		while(1) {
			while(1) {
				
				angle=PI*ran_ecu(seed);
				y=tan(angle);
				em=sq*y+am;
			
				if(!(em < 0.0 || em >= (en+1.0))) break;
			}
			em=floor(em); 
			gam1=gammln(em+1.0);
			gam2=gammln(en-em+1.0);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gam1-gam2+em*plog+(en-em)*pclog);
			rr=ran_ecu(seed);
			if(rr<=t) break;
		}
		
		bnl=em;
	}
   
	if (p != pp) bnl=n-bnl; 
	return bnl;
}

float gammln(float xx) {
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
		return -tmp+log(2.5066282746310005*ser/x);
}

// end functions

// main program

//void noTreeSim(double *times, double *q12, double *q21, int *startState, long *seed, int *nState1, int *nState2){
void noTreeSim(int *ntimes, double *times, double *q12, double *q21, int *startState, long *seed, int *nState1, int *nState2, int *vezes, int *loops){
  int ch12, ch21, i, ns1, ns2, start,loop;
	double wt, pij[4], prT1, rr;
        
        loop = 0;
        start = *startState;

        //ntimes = sizeof(times) / sizeof(double); // This could be wrong. Better to do this in R.
	
	srand(*seed);
	
	if(start==1) {
		ns1=1;
		ns2=0;
	}
	else {
		ns1=0;
		ns2=1;
	}
	
	for(i=1; i< *ntimes; i++) { // Need to verify if this loop will break when i = ntimes or when i = (ntimes-1). And also if this make some difference for the function.
                loop++;
		wt=times[i-1]-times[i];
		// get probs	
		pMk2(pij, *q12, *q21, wt); //Made args 2 and 3 pointers in pMk2 function.
		// how many change?	
		ch12=bnldev(pij[3], ns1, seed);
		ch21=bnldev(pij[2], ns2, seed);
		// add these to the totals
		ns1=ns1-ch12+ch21;
		ns2=ns2-ch21+ch12;
		// which type is affected by the speciation event?
		prT1=(double)ns1/((double)ns1+(double)ns2);
		rr=ran_ecu(seed);
		
		if(rr<prT1) {
			ns1++;
		} else {
			ns2++;
		}
	}
	nState1[0]=ns1;
	nState2[0]=ns2;
	vezes[0]=*ntimes;
	loops[0]=loop;
}	
