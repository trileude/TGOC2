/*
  A simple 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
*/

#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>

#ifndef HMPP
#include "parametres.h"
#include "constoprim.h"
#include "utils.h"

#define CFLOPS(c)		/* {flops+=c;} */
#define IDX(i,j,k)    ( (i*Hstep*Hnxyt) + (j*Hnxyt) + k )
#define IDXE(i,j)     ( (i*Hnxyt) + j )


void
constoprim (const int n,
	    const int Hnxyt,
	    const int Hnvar,
	    const double Hsmallr,
	    const int slices, const int Hstep,
	    double *u, double *q, double *e)
{
  //double u[Hnvar][Hstep][Hnxyt], double q[Hnvar][Hstep][Hnxyt], double e[Hstep][Hnxyt]) {
  //int ijmin, ijmax, IN, i, s;
  //double eken;
  // const int nxyt = Hnxyt;
  WHERE ("constoprim");
  //ijmin = 0;
  //ijmax = n;

  

  int ijmin=0, ijmax=n;
  double eken;
  for (int s = 0; s < slices; s++)
    {
      for (int i = ijmin; i < ijmax; i++)
	    {
	      double qid = MAX (u[IDX (ID, s, i)], Hsmallr);
	      q[IDX (ID, s, i)] = qid;

	      double qiu = u[IDX (IU, s, i)] / qid;
	      double qiv = u[IDX (IV, s, i)] / qid;
	      q[IDX (IU, s, i)] = qiu;
	      q[IDX (IV, s, i)] = qiv;

	      eken = half * (Square (qiu) + Square (qiv));

	      double qip = u[IDX (IP, s, i)] / qid - eken;
	      q[IDX (IP, s, i)] = qip;
	      e[IDXE (s, i)] = qip;

	      CFLOPS (9);
	    }
    }

  if (Hnvar > IP)
  {
      int ijmin=0, ijmax=n;
      for (int IN = IP + 1; IN < Hnvar; IN++)
	    {
	      for (int s = 0; s < slices; s++)
	      {
	        for (int i = ijmin; i < ijmax; i++)
		      {
		        q[IDX (IN, s, i)] = u[IDX (IN, s, i)] / q[IDX (IN, s, i)];
		        CFLOPS (1);
		      }
	      }
	    }
  }
}				// constoprim


#undef IHVW
#undef IDX
#undef IDXE
#endif
//EOF
