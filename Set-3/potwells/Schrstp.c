/*

 SCHRSTP.C  *.MEXW£" file action described in SCHRSTP.M
		 
 Calling syntax:

			[Psi_2, t_2, V] = SCHRSTP(Psi_1, t_1, Vo, DV, Par);
   with:

INPUT QUANTITIES:
	 Psi_1 = input wave function (N-dim)
	 t_1 = starting time	val in femtoec [1e-15 sec]
	 Vo = unperturbed potential energy (real N-dim)
	 DV = potential perturbation	(real N-dim)
	 Par = [dt, omega, cycleNum, xWidth, m_r];  (5 parameters)
		
	 dt = incremental time 
    omega = radian frequency of the perturbation in radian/femtosec;  
 	 cycleNum = number of cycles the incremental solution is computed.
 	 xWidth = width of the wave-propagation interval on X-axis in Angstrom [1e-10 m]
	 m_r = mass relative to electron mass

	OUTPUT QUANTITIES:	 
	 Psi_2 = output wave function
	 t_2 = t_1 + cycleNum*dt, the potential oscillation phase after computation;
	 V = Vo + DV*f(t_2), the potential profile after computation. 

	SCHRSTP.MEXW32 solves the Schroedinger Eq. for a 1-Dim wave function Psi
	defined on a X-axis interval of given width, whose potential
	undergoes a smooth perturbation accoding to the law 
	V(t) = Vo + f(t)*DV, with f(t)=0 for t<0, f(t)=[1-cos(omega*t)]/2
	for 0 < t < pi/omega, f(t)=1 for t > pi/omega. Omega, passed by
	Par[2], is the radian frequency of the perturbation in radian/femtosec;
	The computation is carried out from time t_1 to t_2 = t_1+cycleNum*dt,
   where t_1 is the initial time; dt is elementary the time interval passed by Par(1);  omega, cycleNum,
	passed by Par(3); cycleNum is the number of time the computation is
	performed at each routine call.

  	Author: Renato Nobili - Dept. of Physics of Padova University, November 2010
*/
static char mc_version[]= "MATLAB Compiler 1.0 infun";

#include <math.h>
#include "mex.h"



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

	i*d Psi(x,t)/dt = -(alpha/m_r)*d^2Psi(x,t)/dx^2 + beta*V(x,t)*Psi(x,t)

	with alpha=5.7884, beta=1.5193 governs the motion of a particle
	of mass m_r, relative to electron mass, when the time unit is 1 fsec
	(femtosecond), the space unit is 1 Angstrom and the energy unit is
	1 electronVolt. 
	V(x,t) = Vo(x) + DV2(x)* (1-cos(omega*(t-to)))/2, where t=Tnum*dt 
	The solution is computed using the Caley's form for the infinitesimal
	time translation U(dt):
	
	 Psi(t+dt) = U(dt)*Psi(t),
		
	 with U(dt)=[1-i*H*(dt/2)]/[1+i*H*(dt/2)] and H=Hamiltonian/h_bar; 		
*/

/* ------------------- GLOBAL VARIABLES --------------------------*/
static double alpha = 5.7884; /* = h_bar*c^2/(2*mu) [Angstrom^2/femtosec] */
				 /*   where mu is the electron mass in eV */
static double beta = 1.5193;	 /* 1/h_bar [1/eV*femtosec] */ 	
double c1, c2, c3, coeff;
Complex B, Bstar, B2, tmp;
Complex *Phi, *A, *Astar;
Complex tmp1, tmp2;	
double tMax;
/*-----------------------------------------------------------------*/

void SCHRSTP(Complex *Psi, double *to, double *Vo, double *DV,
	double dt, double omega, int cycleNum, double dx, double m_r, 
   int N, double *V);

void SCHRSTP(Complex *Psi, double *to, double *Vo, double *DV,
	double dt, double omega, int cycleNum, double dx, double m_r, 
   int N, double *V)
{
	register i, j;
	
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

   tMax=pi/omega;
	
	for(j=0; j<N; j++)
	{	
		Astar[j].r = 1;
	};
	
	for(i=0; i<cycleNum; i++)
	{
		*to += dt;
	/*	coeff = 0.5*(1+tanh((*to)*invDt)); */

	
		if(*to < 0)
		{
			coeff=0;
		}
 		else if(*to<tMax)
		{
			coeff = 0.5*(1-cos(omega*(*to))); 
		}
		else 	coeff=1;

		for(j=0; j<N; j++)
		{	
			V[j]= Vo[j] + coeff*DV[j];
			Astar[j].i = -c1-c2*V[j];
		};

		tmp1=CCProd(Astar[0], Psi[0]); /* left hand side setting */ 
		tmp2=CCProd(Bstar, Psi[1]);
		Phi[0]= Sum(tmp1, tmp2);
		A[0]=Conj(Astar[0]);

		for(j=1; j<N-1; j++) /*calculate phi and initialize the diagonal and off diagonal elements */
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

/*---------------------------------------------------------------*/
/* Input Arguments */

#define PSI_IN prhs[0]
#define T_IN prhs[1]
#define POT_IN prhs[2]
#define DPOT prhs[3]
#define PAR prhs[4]

/* Output Arguments */

#define PSI_OUT plhs[0]
#define T_OUT plhs[1]
#define POT_OUT plhs[2]


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  register i;		
  int m, n, mt, nt, mvo, nvo, mdv, ndv, mpar, npar, N, cycleNum;
  double dx, to;	
  double *Pntr, *Pnti, *PntT_1, *PntT_2, *PntV, *Par, *Vo, *DV, *V; 
  Complex *Psi;
	
  /* Check for proper number of arguments */
  
 
  if (nrhs != 5)
  {	
	mexPrintf("\n\n USAGE:[PSI_2, t_2, ,V] = SCHRSTP(PSI_1, t_1, Vo, DV, PAR);\n\n");
  	mexPrintf("	PSI_1 = input wave functions (N-dim complex vectors);\n"); 
	mexPrintf("	t_1 = initial perturbation time femtosec [e.g. -10];\n");
 	mexPrintf("	Vo = unperturbed potential energy in eV (N-dim real vector);\n");
	mexPrintf("	DV = perturbation potential energy in eV (N-dim real vector);\n");
	mexPrintf("	PAR = [dt, omega, cycleNum, xWidth, m_r];\n");
  	mexPrintf("	dt = elementary time interval in femtosec [e.g. 0.001];\n");
  	mexPrintf("	omega = radian frequency of perturbation in radian/fsec [e.g. 10];\n");
	mexPrintf("	cycleNum = number of times the Eq. is solved at each call [e.g. 50].\n"); 
   mexPrintf("	xWidth = x-axis interval in Angstrom [e.g. 20];\n");
  	mexPrintf("	m_r = particle mass relative to electron mass [e.g. 1];\n");
  	mexPrintf("	Call 'help SCHRSTP' from MATLAB command-line for further comments.\n\n");
	mexPrintf("	Author: Renato Nobili - Physics Dept. of Padova University - Italy,  November 2010\n");
	mexPrintf("	PSI_2 = output wave functions (N-dim complex vectors);\n");
	mexPrintf("	t_2 = time after computation (useful for continuation);\n");
	mexPrintf("	V = final potential profile after computation (useful for graphics).\n");
	mexErrMsgTxt("\n	(Function SCHRSTP.MEXW32 requires 4 input arguments)");  
  }
  else if(nlhs > 3)
  {
		mexErrMsgTxt("MexFile SCHRSTP.MEXW32 accepts at most three output argument.");
  };
  
  m = mxGetM(PSI_IN);
  n = mxGetN(PSI_IN);

 /* Assign pointers to the various parameters */
  
  Pntr = mxGetPr(PSI_IN); 

  if(mxIsComplex(PSI_IN)) Pnti = mxGetPi(PSI_IN);

  N = max(n,m);
 
  if(mxIsSparse(POT_IN)||!mxIsDouble(POT_IN) ||
	 mxIsSparse(DPOT)||!mxIsDouble(DPOT) ||
	!mxIsNumeric(PSI_IN) || min(m,n) != 1)
  {
	   mexErrMsgTxt("MexFile SCHRSTP.MEXW32 requires that the first argument be a complex vector.");
  };


/* Check the dimensions of T_IN.  T_IN is a real number. */
  
  mt = mxGetM(T_IN);
  nt = mxGetN(T_IN);
  PntT_1 = mxGetPr(T_IN);
	
  if(mxIsSparse(T_IN)||!mxIsNumeric(T_IN)|| !mxIsDouble(T_IN)|| mt*nt != 1)
  {
	   mexErrMsgTxt("MexFile SCHRSTP.MEXW32 requires that second argument\n is a real number.");
  };

		

  /* Check the dimensions of Vo and DV.  Vo, DV can be N X 1 or 1 X N. */
  
  mvo = mxGetM(POT_IN);
  nvo = mxGetN(POT_IN);
  Vo = mxGetPr(POT_IN);
		
  mdv = mxGetM(DPOT);
  ndv = mxGetN(DPOT);
  DV = mxGetPr(DPOT);
		
  		

 /* Assign pointer to second argument data */
  		
  if(mxIsSparse(POT_IN)||!mxIsNumeric(POT_IN)|| 
		!mxIsDouble(POT_IN)|| min(mvo,nvo)!= 1||max(mvo, nvo)!= N ||
		mxIsSparse(DPOT)||!mxIsNumeric(DPOT)|| 
		!mxIsDouble(DPOT)|| min(mdv,ndv)!= 1||max(mdv, ndv)!= N)
  {
	   mexErrMsgTxt("MexFile SCHRSTP.MEXW32 requires that first,third and fourth argument\n have same dimension.");
  };
	
  mpar = mxGetM(PAR);
  npar = mxGetN(PAR);
  Par = mxGetPr(PAR);
	
  if(!mxIsNumeric(PAR)||!mxIsDouble(PAR)||min(mpar,npar)!=1 || 
		mxIsComplex(PAR) || max(mpar, npar)!=5)
  {
	    mexErrMsgTxt("MexFile SCHRSTP.MEXW32 requires that fifth argument be a 5-dim real vector.");
  };

  if(Par[1] <=0)
  {
	    mexErrMsgTxt("MexFile SCHRSTP.MEXW32 requires that Par[2] is real positive");
  };

  /* Allocate arrays to be passed to the computational routine  */
       
  Psi = (Complex *) mxCalloc(N, sizeof(Complex)); /* allocates memory for PSI */
  			/* mxCalloc frees memory automatically */

  V = (double *) mxCalloc(N, sizeof(double)); /* allocates memory for V */
  			/* mxCalloc frees memory automatically */
	
  /* Assign pointers to the various parameters */ 

  Pntr = mxGetPr(PSI_IN); 

  for(i=0; i<N; i++) /* copy first argument as real vector to complex PSI */
  {
		Psi[i].r = Pntr[i];
  };
  
		
  if(mxIsComplex(PSI_IN))
  {
		Pnti = mxGetPi(PSI_IN);
	
  		for(i=0; i<N; i++) /* copy the first argument to the complex vector PSI */
  		{
			Psi[i].i = Pnti[i];	
  		};
   }
   else
   {
		for(i=0; i<N; i++) /* set the complex vector PSI to zero */
  		{
			Psi[i].i = 0.0;	
  		};
   };
			

/* [dt, Dt,  cycleNum, xWidth, m_r] */

   to=PntT_1[0];
	cycleNum = (int) Par[2];

  	dx= Par[3]/((double) (N-1));  
   	
  
 /* Do the actual computations by the subroutine:
   SCHRSTP(Complex *Psi, double *t_1, double *Vo, double *DV,
	double dt, double Dt, int cycleNum, double dx, double m_r, 
   int N, double *V) */ 

 	SCHRSTP(Psi, &to, Vo, DV, Par[0], Par[1], cycleNum, dx, Par[4], N, V);  
 
/* Psi : pointer to wave function (complex array)
   Vo, DV2  : pointers to unperturbed and perturbation potential (real arrays)
   Par[0] = dt (time interval)
	Par[1] = omega  (radian frequency of perturbation)
	Par[2] = cycleNum
	Par[3] = dx-width (X-width/interval-number)
	Par[4] = particle-mass/electron-mass ratio
   N =	number of points at which Psi and V are computed
	V = returned potential 
 */

/*---------------- COPYING PSI TO OUPUT ----------------*/

  PSI_OUT = mxCreateDoubleMatrix(m, n, mxCOMPLEX);	
  Pntr = mxGetPr(PSI_OUT); 
  Pnti = mxGetPi(PSI_OUT);
  for(i=0; i<N; i++) /* copy PSI components to output matlab vector */
  {
		Pntr[i]=Psi[i].r;
		Pnti[i]=Psi[i].i;	
  };
/*---------------- COPYING to TO OUPUT ----------------*/
  
  T_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
  PntT_2 = mxGetPr(T_OUT); 		
  PntT_2[0] = to;

/*---------------- COPYING V TO OUPUT ----------------*/
 
  POT_OUT = mxCreateDoubleMatrix(mvo, nvo, mxREAL); 
  PntV = mxGetPr(POT_OUT);
  for(i=0; i<N; i++)
  {
		PntV[i]= V[i]; 
  }; 
	
};
