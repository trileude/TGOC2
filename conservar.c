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
#include "utils.h"
#include "conservar.h"

#define CFLOPS(c)		/* {flops+=c;} */

void
gatherConservativeVars (const int idim,
			const int rowcol,
			const int Himin,
			const int Himax,
			const int Hjmin,
			const int Hjmax,
			const int Hnvar,
			const int Hnxt,
			const int Hnyt,
			const int Hnxyt,
			const int slices, const int Hstep,
			double uold[Hnvar * Hnxt * Hnyt], double *u
			//double uold[Hnvar * Hnxt * Hnyt], double u[Hnvar][Hstep][Hnxyt]
  )
{
  //int i, j, ivar, s;

#define IHU(i, j, v)  ((i) + Hnxt * ((j) + Hnyt * (v)))
#define IDX(i,j,k)    ( (i*Hstep*Hnxyt) + (j*Hnxyt) + k )

  WHERE ("gatherConservativeVars");
  if (idim == 1)
    {
      // Gather conservative variables

    for (int s = 0; s < slices; s++)
	  {
	    for (int i = Himin; i < Himax; i++)
	      {
	        int idxuoID = IHU (i, rowcol + s, ID);
	        u[IDX (ID, s, i)] = uold[idxuoID];

	        int idxuoIU = IHU (i, rowcol + s, IU);
	        u[IDX (IU, s, i)] = uold[idxuoIU];

	        int idxuoIV = IHU (i, rowcol + s, IV);
	        u[IDX (IV, s, i)] = uold[idxuoIV];

	        int idxuoIP = IHU (i, rowcol + s, IP);
	        u[IDX (IP, s, i)] = uold[idxuoIP];
	      }
	  }

  if (Hnvar > IP)
  {   
      for (int ivar = IP + 1; ivar < Hnvar; ivar++)
	    {
	      for (int s = 0; s < slices; s++)
		    {
		      for (int i = Himin; i < Himax; i++)
		      {
		        u[IDX (ivar, s, i)] = uold[IHU (i, rowcol + s, ivar)];
		      }
		    }
	    }
	}
      //
    }
  else
    {
      // Gather conservative variables
	
  for (int s = 0; s < slices; s++)
	{
	  for (int j = Hjmin; j < Hjmax; j++)
	    {
	      u[IDX (ID, s, j)] = uold[IHU (rowcol + s, j, ID)];
	      u[IDX (IU, s, j)] = uold[IHU (rowcol + s, j, IV)];
	      u[IDX (IV, s, j)] = uold[IHU (rowcol + s, j, IU)];
	      u[IDX (IP, s, j)] = uold[IHU (rowcol + s, j, IP)];
	    }
	}
  
  if (Hnvar > IP)
	{
	  for (int ivar = IP + 1; ivar < Hnvar; ivar++)
	  {
      for (int s = 0; s < slices; s++)
	    {
	      for (int j = Hjmin; j < Hjmax; j++)
	      {
	       u[IDX (ivar, s, j)] = uold[IHU (rowcol + s, j, ivar)];
	      }
	    }
    }
	}
    }
}

#undef IHU
#undef IDX

void
updateConservativeVars (const int idim,
			const int rowcol,
			const double dtdx,
			const int Himin,
			const int Himax,
			const int Hjmin,
			const int Hjmax,
			const int Hnvar,
			const int Hnxt,
			const int Hnyt,
			const int Hnxyt,
			const int slices, const int Hstep,
			double *uold, double *u,
			double *flux
			//double uold[Hnvar * Hnxt * Hnyt], double u[Hnvar][Hstep][Hnxyt], double flux[Hnvar][Hstep][Hnxyt]
  )
{
  //int i, j, ivar, s;
  WHERE ("updateConservativeVars");

#define IHU(i, j, v)  ((i) + Hnxt * ((j) + Hnyt * (v)))
#define IDX(i,j,k)    ( (i*Hstep*Hnxyt) + (j*Hnxyt) + k )

  if (idim == 1)
    {

      // Update conservative variables
    for (int ivar = 0; ivar <= IP; ivar++)
	  {
	    for (int s = 0; s < slices; s++)
	    {
	      for (int i = Himin + ExtraLayer; i < Himax - ExtraLayer; i++)
		    {
		      uold[IHU (i, rowcol + s, ivar)] =
		            u[IDX (ivar, s, i)] + (flux[IDX (ivar, s, i - 2)] -
					      flux[IDX (ivar, s, i - 1)]) * dtdx;
		      //CFLOPS (3);
		    }
	    }
	  }
    if (Hnvar > IP){
	    for (int ivar = IP + 1; ivar < Hnvar; ivar++){
	      for (int s = 0; s < slices; s++){
		      for (int i = Himin + ExtraLayer; i < Himax - ExtraLayer; i++){
		        uold[IHU (i, rowcol + s, ivar)] =
			            u[IDX (ivar, s, i)] + (flux[IDX (ivar, s, i - 2)] -
					        flux[IDX (ivar, s, i - 1)]) *dtdx;
		        CFLOPS (3);
		      }
		    }
	    }
	  }
  }else{
      // Update conservative variables
    for (int s = 0; s < slices; s++){
      for (int j = Hjmin + ExtraLayer; j < Hjmax - ExtraLayer; j++){
        uold[IHU (rowcol + s, j, ID)] =
	            u[IDX (ID, s, j)] + (flux[IDX (ID, s, j - 2)] -
			       flux[IDX (ID, s, j - 1)]) * dtdx;
        CFLOPS (3);

        uold[IHU (rowcol + s, j, IV)] =
	            u[IDX (IU, s, j)] + (flux[IDX (IU, s, j - 2)] -
			        flux[IDX (IU, s, j - 1)]) * dtdx;
        CFLOPS (3);

        uold[IHU (rowcol + s, j, IU)] =
	            u[IDX (IV, s, j)] + (flux[IDX (IV, s, j - 2)] -
			        flux[IDX (IV, s, j - 1)]) * dtdx;
        CFLOPS (3);

        uold[IHU (rowcol + s, j, IP)] =
	            u[IDX (IP, s, j)] + (flux[IDX (IP, s, j - 2)] -
			        flux[IDX (IP, s, j - 1)]) * dtdx;
        CFLOPS (3);
	    }
	  }

    if (Hnvar > IP){
    
	    for (int ivar = IP + 1; ivar < Hnvar; ivar++){
        for (int s = 0; s < slices; s++){
	        for (int j = Hjmin + ExtraLayer; j < Hjmax - ExtraLayer; j++){
	            uold[IHU (rowcol + s, j, ivar)] =
		            u[IDX (ivar, s, j)] + (flux[IDX (ivar, s, j - 2)] -
				        flux[IDX (ivar, s, j - 1)]) *dtdx;
	            CFLOPS (3);
	        }
	      }
	    }
    }
  }
}

#undef IHU
#undef IDX
#endif
//EOF
