/*
 * CitcomCU is a Finite Element Code that solves for thermochemical
 * convection within a three dimensional domain appropriate for convection
 * within the Earth's mantle. Cartesian and regional-spherical geometries
 * are implemented. See the file README contained with this distribution
 * for further details.
 * 
 * Copyright (C) 1994-2005 California Institute of Technology
 * Copyright (C) 2000-2005 The University of Colorado
 *
 * Authors: Louis Moresi, Shijie Zhong, and Michael Gurnis
 *
 * For questions or comments regarding this software, you may contact
 *
 *     Luis Armendariz <luis@geodynamics.org>
 *     http://geodynamics.org
 *     Computational Infrastructure for Geodynamics (CIG)
 *     California Institute of Technology
 *     2750 East Washington Blvd, Suite 210
 *     Pasadena, CA 91007
 *
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 2 of the License, or any
 * later version.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program; if not, write to the Free Software 
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

void general_stokes_solver(struct All_variables *E)
{
	//float vmag;
	double *delta_U;
	//double *force, Udot_mag, dUdot_mag;
	double Udot_mag, dUdot_mag;
	double time;
	//int count, i, j, k;
	int count, i;

	static double alpha,alpha1;
	static double *oldU;
	static int damp=0,visits = 0;

	//const int nno = E->lmesh.nno;
	const int neq = E->lmesh.neq;
	//const int dims = E->mesh.nsd;
	//const int vpts = vpoints[E->mesh.nsd];

	static int powerlaw; 
	
	if(visits == 0)
	{
	  
	  powerlaw = (E->viscosity.SDEPV || E->viscosity.BDEPV)?(1):(0);
	  if(powerlaw){
	    /* damping factors */
	    alpha = E->viscosity.sdepv_iter_damp;
	    alpha1 = 1 - alpha;
	    if(fabs(alpha-1) > 1e-7){
	      if(E->parallel.me == 0)
		fprintf(stderr,"damping stress dependent iteration velocities by %g\n",alpha);
	      damp = 1;
	    } else{
	      damp = 0;
	    }
	  } /* end powerlaw */
		oldU = (double *)malloc((neq + 2) * sizeof(double));
		for(i = 1; i <= neq; i++)
			oldU[i] = 0.0;
		visits++;
	}

	dUdot_mag = 0.0;

	delta_U = (double *)malloc((neq + 2) * sizeof(double));

	/* FIRST store the old velocity field */

	E->monitor.elapsed_time_vsoln1 = E->monitor.elapsed_time_vsoln;
	E->monitor.elapsed_time_vsoln = E->monitor.elapsed_time;

	if(E->parallel.me == 0)
		time = CPU_time0();

	velocities_conform_bcs(E, E->U);

	assemble_forces(E, 0);

/*	
	if(E->parallel.me==0)
	{
		fprintf(stderr,"time1= %g seconds\n",CPU_time0()-time);
		time=CPU_time0();
	}
*/

	count = 1;

	do
	{

		if(E->viscosity.update_allowed)
			get_system_viscosity(E, 1, E->EVI[E->mesh.levmax], E->VI[E->mesh.levmax]);

		construct_stiffness_B_matrix(E);

		solve_constrained_flow_iterative(E);

		if(powerlaw)
		{

		  if(damp){
		    /* add some of the old solution */
		    for(i = 1; i <= neq; i++)
		      E->U[i] = alpha * E->U[i] + alpha1 * oldU[i];
		  }
			for(i = 1; i <= neq; i++)
			{
				delta_U[i] = E->U[i] - oldU[i];
				oldU[i] = E->U[i];
			}
			Udot_mag = sqrt(global_vdot(E, E->U, E->U, E->mesh.levmax));
			dUdot_mag = sqrt(global_vdot(E, delta_U, delta_U, E->mesh.levmax));

			if(Udot_mag != 0.0)
				dUdot_mag /= Udot_mag;

			if(E->control.sdepv_print_convergence  && (E->parallel.me == 0))
			{
				fprintf(stderr, "Stress dependent viscosity: DUdot = %.4e (%.4e) for iteration %d\n", dUdot_mag, Udot_mag, count);
				fprintf(E->fp, "Stress dependent viscosity: DUdot = %.4e (%.4e) for iteration %d\n", dUdot_mag, Udot_mag, count);
				fflush(E->fp);
			}
			count++;
		}						/* end for SDEPV / BDEPV  */
	} while((count < 50) && (dUdot_mag > E->viscosity.sdepv_misfit) && powerlaw);

	free((void *)delta_U);

	return;
}
