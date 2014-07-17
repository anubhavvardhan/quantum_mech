/* $Revision: 1.1 $ */
/*

  SCHROEQ.C  *.MEXW32" The function is described in SCHROEQ.M
		 Solves the Schroedinger Eq. for a 1-Dim 
		 normalised wave function Psi defined on an interval
		 of axis X of given width and undergoing the action
		 of a given potential energy.

  The calling syntax is:

		Psi_2 = SCHROEQ(Psi_1, V, Par);
   with:
	 Psi_1 = input wave function
	 V = potential energy
	 Par = [dt, m_r, xWidth, N];
	 dt = incremental time interval in femtosec [1e-15 sec]
    m_r = mass relative to electron mass
	 xWidth = width of the wave-propagation interval in Angstrom [1e-10 m]
	 N = number of cycles the incremental solution is computed.

  	Author: Renato Nobili - Dept. of Physics of Padova University, November, 2010
*/
static char mc_version[]= "MATLAB Compiler 1.0 infun";

#include <math.h>
#include "mex.h"


/* Input Arguments */

#define PSI_IN prhs[0]
#define POT prhs[1]
#define PAR prhs[2]

/* Output Arguments */

#define PSI_OUT plhs[0]

#if !defined(max)
#define	max(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(min)
#define	min(A, B)	((A) < (B) ? (A) : (B))
#endif

#define pi 3.14159265

/* ------- Complex invertion function ------------- */

typedef struct
{
	double r;
	double i;
} Complex;

Complex _complex(double *x, double *y);

Complex _complex(double *x, double *y)
{
	Complex z;
	z.r=*x;
	z.i=*y;
	return z;
}

#define complex(x, y) _complex(&x, &y)

Complex _Sum(Complex *x, Complex *y);
Complex _Sum(Complex *x, Complex *y)
{
	Complex z;
	z.r = (*x).r + (*y).r;
	z.i= (*x).i + (*y).i;
	return z;	
};
#define Sum(x, y) _Sum(&x, &y)


Complex _Minus(Complex *x, Complex *y);
Complex _Minus(Complex *x, Complex *y)
{
	Complex z;
	z.r = (*x).r - (*y).r;
	z.i= (*x).i - (*y).i;
	return z;	
};
#define Minus(x, y) _Minus(&x, &y)

Complex _CCProd(Complex *x, Complex *y);
Complex _CCProd(Complex *x, Complex *y)
{
	Complex z;
	z.r = ((*x).r)*((*y).r)-((*x).i)*((*y).i);
	z.i= ((*x).i)*((*y).r)+((*x).r)*((*y).i);
	return z;	
};
#define CCProd(x, y) _CCProd(&x, &y)


Complex _Conj(Complex *c);
Complex _Conj(Complex *c)
{
	Complex z;
	z.r=(*c).r;
	z.i=-((*c).i);
	return z;
};
#define Conj(x) _Conj(&x)


Complex _CInvert(Complex *c);
Complex _CInvert(Complex *c)
{
	Complex z;
	double c2;	
 	c2 = ((*c).r)*((*c).r)+((*c).i)*((*c).i);
	z.r=((*c).r)/c2;
	z.i=-((*c).i)/c2;		
 	return z;
};
#define CInvert(x) _CInvert(&x)

Complex _CCDiv(Complex *x, Complex *y);
Complex _CCDiv(Complex *x, Complex *y)
{
	Complex z;
	double y2;	
 	y2 = ((*y).r)*((*y).r)+((*y).i)*((*y).i);
	z.r=(((*x).r)*((*y).r)+((*x).i)*((*y).i))/y2;
	z.i=(((*x).i)*((*y).r)-((*x).r)*((*y).i))/y2;
	return z;
};
#define CCDiv(x, y) _CCDiv(&x, &y)

/*
	The Schroedinger equation:

	i*d Psi(x,t)/dt = -(alpha/m_r)*d^2Psi(x,t)/dx^2 + beta*V(x)*Psi(x,t)

	with alpha=5.7884, beta=1.5193 governs the motion of a particle
	of mass m_r, relative to electron mass, when the time unit is 1 fsec
	(femtosecond), the space unit is 1 Angstrom and the energy unit is
	1 electronVolt. 
	The solution is computed using the Caley's form for the infinitesimal
	time translation U(dt):
	
	 Psi(t+dt) = U(dt)*Psi(t),
		
	 with U(dt)=[1-i*H*(dt/2)]/[1+i*H*(dt/2)] and H=Hamiltonian/h_bar; 		
*/

static double alpha = 5.7884; /* = h_bar*c^2/(2*mu) [Angstrom^2/femtosec] */
				 /*   where mu is the electron mass in eV */
static double beta = 1.5193;	 /* 1/h_bar [1/eV*femtosec] */ 	

void SCHROEQ(Complex *Psi, double *V, double dt, double dx,
	double m_r, unsigned int N, int Tnum);

void SCHROEQ(Complex *Psi, double *V, double dt, double dx,
	double m_r, unsigned int N, int Tnum)
{
	int i, j;
	double c1, c2, c3;

	Complex B, Bstar, B2, tmp;
	Complex *Phi, *A, *Astar;
	Complex tmp1, tmp2;	

	c1 = dt*alpha/(m_r*dx*dx);
	c2 = 0.5*beta*dt;
	c3 = 0.5*c1;
	B.r = 0;
	B.i = -c3;

	Bstar=Conj(B);
	B2 = CCProd(B,B);	
	
	Phi =(Complex *)mxCalloc(N,sizeof(Complex));
	A =(Complex *) mxCalloc(N, sizeof(Complex));
	Astar =(Complex *) mxCalloc(N, sizeof(Complex));
	/* mxCalloc frees memory automatically */

	for(j=0; j<N; j++)
	{	
		Astar[j].r = 1;
		Astar[j].i= -c1-c2*V[j];
	};
	
	for(i=0; i<Tnum; i++)
	{	
		tmp1=CCProd(Astar[0], Psi[0]);
		tmp2=CCProd(Bstar, Psi[1]);
		Phi[0]= Sum(tmp1, tmp2);
		A[0]=Conj(Astar[0]);

		for(j=1; j<N-1; j++) //calculate phi and initialize the diagonal and off diagonal elements
		{
			tmp1=CCProd(Astar[j],Psi[j]); 
			tmp= Sum(Psi[j+1], Psi[j-1]);
			tmp2=CCProd(Bstar,tmp); 
			Phi[j]= Sum(tmp1, tmp2);
			A[j]=Conj(Astar[j]);
		};
		tmp1=CCProd(Astar[N-1],Psi[N-1]);
		tmp2=CCProd(Bstar, Psi[N-2]);
		Phi[N-1]=Sum(tmp1, tmp2);
		A[N-1]=Conj(Astar[N-1]);
		
		for(j=1; j<N; j++) /*sweep down the index to do the LU decomposition */
		{	
			tmp = CInvert(A[j-1]);
			tmp2 =CCProd(B,tmp);
			tmp1 = CCProd(tmp2, Phi[j-1]);
			Phi[j] = Minus(Phi[j], tmp1);
			tmp2=CCProd(B2,tmp);
 			A[j] = Minus(A[j], tmp2);
		};

		for(j=N-2; j>0; j--) /*sweep up the index to do the back substitution*/
		{
			tmp1=CCProd(B,Psi[j+1]);
			tmp = Minus(Phi[j],tmp1);
			Psi[j] = CCDiv(tmp, A[j]);
		};
	};
	
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  int i;		
  unsigned int m, n, m2, n2, N;
  int cycleNum;
  double dx;	
  double *Pntr, *Pnti, *Par, *V; 
  Complex *Psi;
	
  /* Check for proper number of arguments */
  
 
  if (nrhs != 3)
  {	
	mexPrintf("\n\n USAGE:    PSI_2  =  SCHROEQ(PSI_1, V, PAR);\n\n");
  	mexPrintf("	PSI_1, PSI_2 = input and output wave functions (N-dim complex vectors);\n");
  	mexPrintf("	V = potential energy profile in eV (N-dim real vector);\n");
	mexPrintf("	PAR = [dt, xWidth, m_r, cycleNum];\n");
  	mexPrintf("	dt = elementary time interval in femtosec [e.g. 0.001];\n");
  	mexPrintf("	xWidth = x-axis interval in Angstrom [e.g. 20];\n");
  	mexPrintf("	m_r = particle mass relative to electron mass [e.g. 1];\n");
  	mexPrintf("	cycleNum = number of times the Eq. is solved at each call [e.g. 50].\n"); 
	mexPrintf("	Call 'help SCHROEQ' from MATLAB command-line for other instructions.\n\n");
	mexPrintf("	Author: Renato Nobili - Physics Dept. of Padova University - Italy,  November 1999\n");
	 	mexErrMsgTxt("\n	(MexFile SCHROEQ.MEXW32 requires 3 input arguments)");
  }
  else if(nlhs > 1)
  {
		mexErrMsgTxt("MexFile SCHROEQ.MEXW32 accepts at most one output argument.");
  };
  
  m = mxGetM(PSI_IN);
  n = mxGetN(PSI_IN);

 /* Assign pointers to the various parameters */
  
  Pntr = mxGetPr(PSI_IN); 

  if(mxIsComplex(PSI_IN)) Pnti = mxGetPi(PSI_IN);

  N=max(n,m);
 
  if(mxIsSparse(POT)||!mxIsNumeric(PSI_IN)||!mxIsDouble(POT) || min(m,n) != 1)
  {
	   mexErrMsgTxt("MexFile SCHROEQ.MEXW32 requires that the first argument be a complex vector.");
  };

  /* Check the dimensions of V.  V can be N X 1 or 1 X N. */
  
  m2 = mxGetM(POT);
  n2 = mxGetN(POT);

 /* Assign pointer to second argument data */
  V = mxGetPr(POT);
	
  if(mxIsSparse(POT)||!mxIsNumeric(POT)|| !mxIsDouble(POT)|| min(m2,n2)!= 1||max(m2, n2)!= N)
  {
	   mexErrMsgTxt("MexFile SCHROEQ.MEXW32 requires that first and second arguments\n have same dimension.");
  };


	
  m2 = mxGetM(PAR);
  n2 = mxGetN(PAR);
  Par = mxGetPr(PAR);
	
  if(mxIsSparse(POT)||!mxIsNumeric(PAR)||!mxIsDouble(PAR)||min(m2,n2)!=1 || max(m2, n2)!=4)
  {
	    mexErrMsgTxt("MexFile SCHROEQ.MEXW32 requires that third argument be a 4-dim real vector.");
  };

  /* Create a matrix for the return argument */
       
  Psi = (Complex *) mxCalloc(N, sizeof(Complex)); /* allocates memory for PSI */
  			/* mxCalloc frees memory automatically */

  /* Assign pointers to the various parameters */ 
  Pntr = mxGetPr(PSI_IN); 
  for(i=0; i<N; i++) /* copy first argument as complex vector to complex PSI */
  {
		Psi[i].r = Pntr[i];
  };

  if(mxIsComplex(PSI_IN))
  {
		Pnti = mxGetPi(PSI_IN);
	
  		for(i=0; i<N; i++) /* copy first argument as complex vector to complex PSI */
  		{
			Psi[i].i = Pnti[i];	
  		};
   }; 

  	dx= Par[1]/((double) (N-1));  

  	cycleNum = (int) Par[3];
 
 /* Do the actual computations in a subroutine */
  
 	SCHROEQ(Psi, V, Par[0], dx, Par[2], N, cycleNum);  
 
/* Psi : pointer to wave function (complex array)
   V : pointer to potential (real array)
   Par[0]: dT (time interval)
   Par[1] = dx-width (X-width/interval-number)
   Par[2] = particle mass
   N :	number of points at which Psi and V are computed
   Par[3]: number of times the wave function is incremented
*/

  PSI_OUT = mxCreateDoubleMatrix(m, n, mxCOMPLEX);

  Pntr = mxGetPr(PSI_OUT); 
  Pnti = mxGetPi(PSI_OUT);


  for(i=0; i<N; i++) /* copy PSI components to output matlab vector */
  {
	Pntr[i]=Psi[i].r;
	Pnti[i]=Psi[i].i;	
  }

};


