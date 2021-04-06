#include <stdio.h>
#include <math.h>

#define  NMAX_200M  7500    //Maximum number of input data ordinates
#define  NP_200M      350    //Maximum number of filter coefficients
#define  MMAX_200M     4    //Maximum order of smoothing polynomial

//Note: index 0 not used here.
typedef  float MAT_200M[MMAX_200M+2][MMAX_200M+2];

float signal_200m[NMAX_200M+1], ysave_200m[NMAX_200M+1];
float c_200m[NP_200M+1];
int   index_pt_200m[NP_200M+1];

int   i,j,m,ndata,nl,nr;
float dt,t,tbegin,temp,tend;

FILE  *fp_in, *fp_out;

int IMin_200M(int ia, int ib) 
{
  if (ia<=ib) return ia;
  else return ib;
}


  void LUDCMP_200M(MAT_200M A, int N, int np, int *INDX, int *D, int *CODE);
  void LUBKSB_200M(MAT_200M A, int N, int np, int *INDX, float *B);


  void savgol_200m(float *c_200m, int np, int nl, int nr, int ld, int m)  
  {
    int d,icode,imj,ipj,j,k,kk,mm,i;
    int indx[MMAX_200M+2];
    float fac, sum;
    MAT_200M   a;
    float b[MMAX_200M+2];
    
    if (np<nl+nr+1 || nl<0 || nr<0 || ld>m || m>MMAX_200M || nl+nr<m) 
      {
	printf("\n Bad args in savgol.\n");
	return;
      }
    
    for (i=1; i<=MMAX_200M+1; i++) 
      {
	for (j=1; j<=MMAX_200M+1; j++) a[i][j]=0.0;
	b[i]=0.0;
	indx[i]=0;
      }
    
    for (ipj=0; ipj<=2*m; ipj++) 
      { //Set up the normal equations of the desired leastsquares fit.
	sum=0.0;
	if (ipj==0) sum=1.0;
	for (k=1; k<=nr; k++) sum += (float) pow(k,ipj);
	for (k=1; k<=nl; k++) sum += (float) pow(-k,ipj);
	mm=IMin_200M(ipj,2*m-ipj);
	imj=-mm;
	do 
	  {
	    a[1+(ipj+imj)/2][1+(ipj-imj)/2]=sum;
	    imj += 2;
	  } while (imj<=mm); 
      }
    
    LUDCMP_200M(a,m+1,MMAX_200M+1,indx,&d,&icode);    //Solve them: LU decomposition
    
    for (j=1; j<=m+1; j++) b[j]=0.0;
    b[ld+1]=1.0;    //Right-hand side vector is unit vector, depending on which derivative we want.
    
    LUBKSB_200M(a,m+1,MMAX_200M+1,indx,b);      //Backsubstitute, giving one row of the inverse matrix.
    
    for (kk=1; kk<=np; kk++)          //Zero the output array (it may be bigger than the number
      c_200m[kk]=0.0;                      //of coefficients.
    
    for (k=-nl; k<=nr; k++) 
      {         //Each Savitzky-Golay coefficient is the dot product
	sum=b[1];                       //of powers of an integer with the inverse matrix row.
	fac=1.0;
	for (mm=1; mm<=m; mm++) 
	  {
	    fac *= k;
	    sum += b[mm+1]*fac;
	  } 
	kk=((np-k) % np) + 1;           //Store in wrap-around order}
	c_200m[kk]=sum;
      }
  }


void LUDCMP_200M(MAT_200M A, int N, int np, int *INDX, int *D, int *CODE) 
{
  
#define NMX  100
  
  float AMAX,DUM,SUM,TINY;
  float VV[NMX];
  int   I,IMAX,J,K;
  IMAX=0;
  TINY = (float)1e-12;
  
  *D=1; *CODE=0;
  
  for (I=1; I<=N; I++) 
    {
      AMAX=0.0;
      for (J=1; J<=N; J++)
	if (fabs(A[I][J]) > AMAX) AMAX=(float) fabs(A[I][J]);
      if (AMAX < TINY) 
	{
	  *CODE = 1;
	  return;
	}
      VV[I] = (float)1.0 / AMAX;
    }
  
  for (J=1; J<=N; J++) 
    {
      for (I=1; I<J; I++) 
	{
	  SUM = A[I][J];
	  for (K=1; K<I; K++)
	    SUM -= A[I][K]*A[K][J];
	  A[I][J] = SUM;
	}
      AMAX = 0.0;
      for (I=J; I<=N; I++) 
	{
	  SUM = A[I][J];
	  for (K=1; K<J; K++)
	    SUM -= A[I][K]*A[K][J];
	  A[I][J] = SUM;
	  DUM = VV[I]*(float)fabs(SUM);
	  if (DUM >= AMAX) 
	    {
	      IMAX = I;
	      AMAX = DUM;
	    }
	}  
      
      if (J != IMAX) 
	{
	  for (K=1; K<=N; K++) 
	    {
	      DUM = A[IMAX][K];
	      A[IMAX][K] = A[J][K];
	      A[J][K] = DUM;
	    }
	  *D = -(*D);
	  VV[IMAX] = VV[J];
	}
      
      INDX[J] = IMAX;
      if ((float)fabs(A[J][J]) < TINY)  A[J][J] = TINY;
      
      if (J != N) 
	{
	  DUM = (float)1.0 / A[J][J];
	  for (I=J+1; I<=N; I++)  A[I][J] *= DUM;
	} 
    } //j loop
  
} //LUDCMP_200M()



void LUBKSB_200M(MAT_200M A, int N, int np, int *INDX, float *B)  
{
  float SUM;
  int I,II,J,LL;
  
  II = 0;
  
  for (I=1; I<=N; I++) 
    {
      LL = INDX[I];
      SUM = B[LL];
      B[LL] = B[I];
      if (II != 0)
	for (J=II; J<I; J++)
	  SUM -= A[I][J]*B[J];
      else if (SUM != 0.0)
	II = I;
      B[I] = SUM;
    }
  
  for (I=N; I>0; I--) 
    {
      SUM = B[I];
      if (I < N)
	for (J=I+1; J<=N; J++)
	  SUM -= A[I][J]*B[J];
      B[I] = SUM / A[I][I];
    }
  
}
