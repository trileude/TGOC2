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
#include <strings.h>
#include <string.h>

#include "hmpp.h"
#include "parametres.h"
#include "hydro_godunov.h"
#include "hydro_funcs.h"
#include "utils.h"
#include "make_boundary.h"

#include "riemann.h"
#include "qleftright.h"
#include "trace.h"
#include "slope.h"
#include "equation_of_state.h"
#include "constoprim.h"

#include "cmpflx.h"
#include "conservar.h"

// variables auxiliaires pour mettre en place le mode resident de HMPP
void
hydro_godunov (int idimStart, double dt, const hydroparam_t H,
	       hydrovar_t * Hv, hydrowork_t * Hw, hydrovarwork_t * Hvw)
{
  // Local variables
  int j;
  double dtdx;
  int clear = 0;

  double (*e)[H.nxyt];
  double (*flux)[H.nxystep][H.nxyt];
  double (*qleft)[H.nxystep][H.nxyt];
  double (*qright)[H.nxystep][H.nxyt];
  double (*c)[H.nxyt];
  double *uold;
  int (*sgnm)[H.nxyt];
  double (*qgdnv)[H.nxystep][H.nxyt];
  double (*u)[H.nxystep][H.nxyt];
  double (*qxm)[H.nxystep][H.nxyt];
  double (*qxp)[H.nxystep][H.nxyt];
  double (*q)[H.nxystep][H.nxyt];
  double (*dq)[H.nxystep][H.nxyt];

  static FILE *fic = NULL;

  if (fic == NULL && H.prt == 1)
    {
      char logname[256];
      sprintf (logname, "TRACE.%04d_%04d.txt", H.nproc, H.mype);
      fic = fopen (logname, "w");
    }

  WHERE ("hydro_godunov");

  // int hmppGuard = 1;
  int idimIndex = 0;
  // Allocate work space for 1D sweeps
  allocate_work_space (H.nxyt, H, Hw, Hvw);

  for (idimIndex = 0; idimIndex < 2; idimIndex++)
    {
      int idim = (idimStart - 1 + idimIndex) % 2 + 1;
      // constant
      dtdx = dt / H.dx;

      // Update boundary conditions
      if (H.prt)
	{
	  fprintf (fic, "godunov %d\n", idim);
	  PRINTUOLD (fic, H, Hv);
	}
      // if (H.mype == 1) fprintf(fic, "Hydro makes boundary.\n");
      make_boundary (idim, H, Hv);
      if (H.prt)
	{
	  fprintf (fic, "MakeBoundary\n");
	}
      PRINTUOLD (fic, H, Hv);

      uold = Hv->uold;
      qgdnv = (double (*)[H.nxystep][H.nxyt]) Hvw->qgdnv;
      flux = (double (*)[H.nxystep][H.nxyt]) Hvw->flux;
      c = (double (*)[H.nxyt]) Hw->c;
      e = (double (*)[H.nxyt]) Hw->e;
      qleft = (double (*)[H.nxystep][H.nxyt]) Hvw->qleft;
      qright = (double (*)[H.nxystep][H.nxyt]) Hvw->qright;
      sgnm = (int (*)[H.nxyt]) Hw->sgnm;
      q = (double (*)[H.nxystep][H.nxyt]) Hvw->q;
      dq = (double (*)[H.nxystep][H.nxyt]) Hvw->dq;
      u = (double (*)[H.nxystep][H.nxyt]) Hvw->u;
      qxm = (double (*)[H.nxystep][H.nxyt]) Hvw->qxm;
      qxp = (double (*)[H.nxystep][H.nxyt]) Hvw->qxp;

      int Hmin, Hmax, Hstep;
      int Hdimsize;
      int Hndim_1;

      if (idim == 1)
	{
	  Hmin = H.jmin + ExtraLayer;
	  Hmax = H.jmax - ExtraLayer;
	  Hdimsize = H.nxt;
	  Hndim_1 = H.nx + 1;
	  Hstep = H.nxystep;
	}
      else
	{
	  Hmin = H.imin + ExtraLayer;
	  Hmax = H.imax - ExtraLayer;
	  Hdimsize = H.nyt;
	  Hndim_1 = H.ny + 1;
	  Hstep = H.nxystep;
	}

      if (!H.nstep && idim == 1)
	{
	  /* LM -- HERE a more secure implementation should be used: a new parameter ? */
	}
      // if (H.mype == 1) fprintf(fic, "Hydro computes slices.\n");


  


   
  for (j = Hmin; j < Hmax; j += Hstep)
	{
	  // we try to compute many slices each pass
	  int jend = j + Hstep;
	  if (jend >= Hmax)
	    jend = Hmax;
	  int slices = jend - j;	// numbre of slices to compute
	  // fprintf(stderr, "Godunov idim=%d, j=%d %d \n", idim, j, slices);

	  if (clear)
	    Dmemset ((H.nxyt) * H.nxystep * H.nvar, (double *) dq, 0);
	  gatherConservativeVars (idim, j, H.imin, H.imax, H.jmin, H.jmax,
				  H.nvar, H.nxt, H.nyt, H.nxyt, slices, Hstep,
				  uold, (double*)u);
	  if (H.prt)
	    {
	      fprintf (fic, "ConservativeVars %d %d %d %d %d %d\n", H.nvar,
		       H.nxt, H.nyt, H.nxyt, slices, Hstep);
	    }
	  PRINTARRAYV2 (fic, u, Hdimsize, "u", H);

	  if (clear)
	    Dmemset ((H.nxyt) * H.nxystep * H.nvar, (double *) dq, 0);

	  // Convert to primitive variables
	  constoprim (Hdimsize, H.nxyt, H.nvar, H.smallr, slices, Hstep, (double*)u, (double*)q,(double*) e);
	  PRINTARRAY (fic, e, Hdimsize, "e", H);
	  PRINTARRAYV2 (fic, q, Hdimsize, "q", H);

	  equation_of_state (0, Hdimsize, H.nxyt, H.nvar, H.smallc, H.gamma,
			     slices, Hstep, (double*)e, (double*)q, (double*)c);
	  PRINTARRAY (fic, c, Hdimsize, "c", H);
	  PRINTARRAYV2 (fic, q, Hdimsize, "q", H);

	  // Characteristic tracing
	  if (H.iorder != 1)
	    {
	      slope (Hdimsize, H.nvar, H.nxyt, H.slope_type, slices, Hstep, (double*)q,		     (double*)dq);
	      PRINTARRAYV2 (fic, dq, Hdimsize, "dq", H);
	    }

	  if (clear)
	    Dmemset ((H.nxyt + 2) * H.nxystep * H.nvar, (double *) qxm, 0);
	  if (clear)
	    Dmemset ((H.nxyt + 2) * H.nxystep * H.nvar, (double *) qxp, 0);
	  if (clear)
	    Dmemset ((H.nxyt + 2) * H.nxystep * H.nvar, (double *) qleft, 0);
	  if (clear)
	    Dmemset ((H.nxyt + 2) * H.nxystep * H.nvar, (double *) qright, 0);
	  if (clear)
	    Dmemset ((H.nxyt + 2) * H.nxystep * H.nvar, (double *) flux, 0);
	  if (clear)
	    Dmemset ((H.nxyt + 2) * H.nxystep * H.nvar, (double *) qgdnv, 0);
	  trace (dtdx, Hdimsize, H.scheme, H.nvar, H.nxyt, slices, Hstep, (double*)q,
		 (double*)dq, (double*)c, (double*)qxm, (double*)qxp);
	  PRINTARRAYV2 (fic, qxm, Hdimsize, "qxm", H);
	  PRINTARRAYV2 (fic, qxp, Hdimsize, "qxp", H);

	  qleftright (idim, H.nx, H.ny, H.nxyt, H.nvar, slices, Hstep, (double*)qxm,
		      (double*)qxp, (double*)qleft, (double*)qright);
	  PRINTARRAYV2 (fic, qleft, Hdimsize, "qleft", H);
	  PRINTARRAYV2 (fic, qright, Hdimsize, "qright", H);

	  riemann (Hndim_1, H.smallr, H.smallc, H.gamma, H.niter_riemann,
		   H.nvar, H.nxyt, slices, Hstep, (double*)qleft, (double*)qright, (double*)qgdnv, (int*)sgnm);
	  PRINTARRAYV2 (fic, qgdnv, Hdimsize, "qgdnv", H);

	  cmpflx (Hdimsize, H.nxyt, H.nvar, H.gamma, slices, Hstep, (double*)qgdnv,
		  (double*)flux);
	  PRINTARRAYV2 (fic, flux, Hdimsize, "flux", H);
	  PRINTARRAYV2 (fic, u, Hdimsize, "u", H);

	  updateConservativeVars (idim, j, dtdx, H.imin, H.imax, H.jmin,
				  H.jmax, H.nvar, H.nxt, H.nyt, H.nxyt,
				  slices, Hstep, (double*)uold, (double*)u, (double*)flux);


	  PRINTUOLD (fic, H, Hv);
	  }			// for j

      

      if (H.prt)
	{
	  // printf("[%d] After pass %d\n", H.mype, idim);
	  PRINTUOLD (fic, H, Hv);
	}
    }
  // Deallocate work space
  deallocate_work_space (H, Hw, Hvw);
  if ((H.t + dt >= H.tend) || (H.nstep + 1 >= H.nstepmax))
    {
      /* LM -- HERE a more secure implementation should be used: a new parameter ? */
    }

}				// hydro_godunov


// EOF
