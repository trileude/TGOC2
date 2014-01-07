/*
  A simple 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
*/

// #include <stdlib.h>
// #include <unistd.h>
#include <math.h>
#include <stdio.h>

#ifndef HMPP
#include "equation_of_state.h"
#include "parametres.h"
#include "utils.h"

#define CFLOPS(c)		/* {flops+=c;} */
#define IDX(i,j,k)    ( (i*Hstep*Hnxyt) + (j*Hnxyt) + k )
#define IDXE(i,j)     ( (i*Hnxyt) + j )

void
equation_of_state (int imin,
		   int imax,
		   const int Hnxyt,
		   const int Hnvar,
		   const double Hsmallc,
		   const double Hgamma,
		   const int slices, const int Hstep,
		   double *eint, double *q, double *c)
{
  //double eint[Hstep][Hnxyt], double q[Hnvar][Hstep][Hnxyt], double c[Hstep][Hnxyt]) {
  //int k, s;
  //double smallp;

  WHERE ("equation_of_state"); 
  //smallp = Square (Hsmallc) / Hgamma;
  //CFLOPS (1);

    double smallp = Square (Hsmallc) / Hgamma;
    CFLOPS (1);
    for (int s = 0; s < slices; s++)
    {
      for (int k = imin; k < imax; k++)
	    {
	      double rhok = q[IDX (ID, s, k)];
	      double base = (Hgamma - one) * rhok * eint[IDXE (s, k)];
	      base = MAX (base, (double) (rhok * smallp));

	      q[IDX (IP, s, k)] = base;
	      c[IDXE (s, k)] = sqrt (Hgamma * base / rhok);

	      CFLOPS (7);
	    }
    }
}				// equation_of_state

#undef IDX
#undef IDXE

#endif
// EOF
