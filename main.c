/*
  A simple 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
*/
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#include "parametres.h"
#include "hydro_funcs.h"
#include "vtkfile.h"
#include "compute_deltat.h"
#include "hydro_godunov.h"
#include "utils.h"
hydroparam_t H;
hydrovar_t Hv;			// nvar
hydrovarwork_t Hvw;		// nvar
hydrowork_t Hw;
unsigned long flops = 0;
int
main (int argc, char **argv)
{
  double dt = 0;
  int nvtk = 0;
  char outnum[80];
  int time_output = 0;
  int nthreads = -1;

  // double output_time = 0.0;
  double next_output_time = 0;
  double start_time = 0, end_time = 0;
  double start_iter = 0, end_iter = 0;
  double elaps = 0;

  omp_set_num_threads(12);
  #pragma omp parallel
  {
    #pragma omp master
    {
      nthreads = omp_get_num_threads();
      printf("nthreads : %d\n", nthreads);
    }
  }

  MPI_Init (&argc, &argv);

  start_time = cclock ();
  if (H.mype == 1)
    fprintf (stdout, "Hydro starts.\n");
  process_args (argc, argv, &H);
  hydro_init (&H, &Hv);
  // PRINTUOLD(H, &Hv);
  if (H.nproc > 1)
    MPI_Barrier (MPI_COMM_WORLD);

  if (H.dtoutput > 0)
    {

      // outputs are in physical time not in time steps
      time_output = 1;
      next_output_time = next_output_time + H.dtoutput;
    }
  if (H.dtoutput || H.noutput)
    vtkfile (++nvtk, H, &Hv);
  if (H.mype == 1)
    fprintf (stdout, "Hydro starts main loop.\n");
  while ((H.t < H.tend) && (H.nstep < H.nstepmax))
    {
      start_iter = cclock ();
      outnum[0] = 0;
      flops = 0;
      if ((H.nstep % 2) == 0)
	{
	  // if (H.mype == 1) fprintf(stdout, "Hydro computes deltat.\n");
	  compute_deltat (&dt, H, &Hw, &Hv, &Hvw);
	  if (H.nstep == 0)
	    {

	      dt = dt / 2.0;
	    }
	  if (H.nproc > 1)
	    {
	      double dtmin;
	      MPI_Allreduce (&dt, &dtmin, 1, MPI_DOUBLE, MPI_MIN,
			     MPI_COMM_WORLD);
	      dt = dtmin;
	    }
	}
      // if (H.mype == 1) fprintf(stdout, "Hydro starts godunov.\n");
      if ((H.nstep % 2) == 0)
	    {
	      hydro_godunov (1, dt, H, &Hv, &Hw, &Hvw);
    //            hydro_godunov(2, dt, H, &Hv, &Hw, &Hvw);
	    }
      else
	    {
	      hydro_godunov (2, dt, H, &Hv, &Hw, &Hvw);
    //            hydro_godunov(1, dt, H, &Hv, &Hw, &Hvw);
	    }
      end_iter = cclock ();
      H.nstep++;
      H.t += dt;
      {
	double iter_time = (double) (end_iter - start_iter);
	if (flops > 0)
	  {
	    if (iter_time > 1.e-9)
	      {
		double mflops = (double) flops / (double) 1.e+6 / iter_time;
		sprintf (outnum, "%s {%.3f Mflops %lu} (%.3fs)", outnum,
			 mflops, flops, iter_time);
	      }
	  }
	else
	  {
	    if (H.nx == 400 && H.ny == 400)
	      {			/* LM -- Got from input !! REMOVE !!  */
		flops = 31458268;
		double mflops = (double) flops / (double) 1.e+6 / iter_time;
		sprintf (outnum, "%s {~%.3f Mflops} (%.3fs)", outnum, mflops,
			 iter_time);
	      }
	    else
	      sprintf (outnum, "%s (%.3fs)", outnum, iter_time);
	  }
      }
      if (time_output == 0)
	{
	  if ((H.nstep % H.noutput) == 0)
	    {
	      vtkfile (++nvtk, H, &Hv);
	      sprintf (outnum, "%s [%04d]", outnum, nvtk);
	    }
	}
      else
	{
	  if (H.t >= next_output_time)
	    {
	      vtkfile (++nvtk, H, &Hv);
	      next_output_time = next_output_time + H.dtoutput;
	      sprintf (outnum, "%s [%04d]", outnum, nvtk);
	    }
	}
      if (H.mype == 0)
	{
	  fprintf (stdout, "--> Step=%4d, %12.5e, %10.5e %s\n", H.nstep, H.t,
		   dt, outnum);
	  fflush (stdout);
	}
    }
  hydro_finish (H, &Hv);
  end_time = cclock ();
  elaps = (double) (end_time - start_time);
  timeToString (outnum, elaps);
  if (H.mype == 0)
    fprintf (stdout, "Hydro ends in %ss (%.3lf).\n", outnum, elaps);
  MPI_Finalize ();
  return 0;
}
