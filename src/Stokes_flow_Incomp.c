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

/*   Functions which solve for the velocity and pressure fields using Uzawa-type iteration loop.  */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

extern int Emergency_stop;

/* Master loop for pressure and (hence) velocity field */


void solve_constrained_flow_iterative(struct All_variables *E)
{
	//double *D1;
	//double *u;
	//double *R, *Bp;
	double residual_ddash;
	//double vmag;

	static int been_here = 0;

	//int steps, cycles;
	int cycles;
	//int i, j, k, doff, vel_cycles_previous, vel_calls_previous;

	double time;

	//const int npno = E->lmesh.npno;
	//const int gnpno = E->mesh.npno;
	//const int nno = E->lmesh.nno;
	//const int dims = E->mesh.nsd;
	//const int neq = E->lmesh.neq;
	//const int gneq = E->mesh.neq;

	time = CPU_time0();

	cycles = E->control.p_iterations;

	/* Solve for velocity and pressure, correct for bc's */

	residual_ddash = solve_Ahat_p_fhat(E, E->U, E->P, E->F, E->control.accuracy, &cycles);

	been_here = 1;

	v_from_vector(E, E->V, E->U);
#ifdef USE_GZDIR
	if(E->control.gzdir)
	  p_to_nodes(E,E->P,E->NP,E->mesh.levmax); 
#endif

	return;
}



/*  ==========================================================================  

adjusted CitcomS style (doesn't quite work yet)


*/

float solve_Ahat_p_fhat_new(struct All_variables *E, double *V, double *P, 
			    double *FF, double imp, int *steps_max)
{
  int i, j, count, convergent, valid, problems, lev, npno, neq;
  int gnpno, gneq;
  
  static int been_here = 0;
  static double *r1, *r2, *z1, *s1, *s2, *u1;
  double *F;
  double *shuffle;
  double alpha, delta, r0dotr0, r1dotz1, r0dotz0;
  double residual, initial_residual, res_magnitude, v_res;
  char message[500];
  double time0, time;
  static double timea;
  float dpressure, dvelocity;

  npno = E->lmesh.npno;
  neq = E->lmesh.neq;
  lev = E->mesh.levmax;  
  gnpno = E->mesh.npno;
  gneq = E->mesh.neq;
  
  if(been_here == 0)  {
    r1 = (double *)malloc((npno + 1) * sizeof(double));
    r2 = (double *)malloc((npno + 1) * sizeof(double));
    z1 = (double *)malloc((npno + 1) * sizeof(double));
    s1 = (double *)malloc((npno + 1) * sizeof(double));
    s2 = (double *)malloc((npno + 1) * sizeof(double));
    u1 = (double *)malloc((neq + 2) * sizeof(double));
    F = (double *)malloc((neq+2)*sizeof(double));
   
    timea = CPU_time0();
  }
  been_here++;
  
  problems = 0;
  time0 = time = CPU_time0();

  /* copy the F vector */
  for(j=0;j<neq;j++)
    F[j] = FF[j];

  v_res = sqrt(global_vdot(E, E->F, E->F, lev));


  initial_vel_residual(E, V, P, F, u1,v_res);

  assemble_div_u(E, V, r1, lev);

  E->monitor.incompressibility = sqrt(global_div_norm2(E, r1)
				      / (1e-32 + E->monitor.vdotv));
  
  dvelocity = 1.0;
  dpressure = 1.0;
  convergent = 0;
  r0dotz0 = 0;

  if(E->parallel.me == 0)
    fprintf(stderr, "initial residue of momentum equation %g %d inc %g\n", 
	    v_res, gneq,E->monitor.incompressibility );

  count = 0;
  E->monitor.vdotv = global_vdot(E, V, V, lev);
  E->monitor.pdotp = global_pdot(E, P, P, lev);
  E->monitor.vdotv = sqrt(E->monitor.vdotv/gneq);
  E->monitor.pdotp = sqrt(E->monitor.pdotp/gnpno);

  generate_log_message(count,time0,timea,dvelocity, dpressure,E);
  residual = 0;
  while((count < *steps_max) && (dpressure >= imp || dvelocity >= imp)){

    for(j = 1; j <= npno; j++)
      z1[j] = E->BPI[lev][j] * r1[j];
    
    r1dotz1 = global_pdot(E, r1, z1, lev);
    assert(r1dotz1 != 0.0  /* Division by zero in head of incompressibility iteration */);

    if((count == 0))
      for(j = 1; j <= npno; j++)
	s2[j] = z1[j];
    else{
      delta = r1dotz1 / r0dotr0;
      for(j = 1; j <= npno; j++)
	s2[j] = z1[j] + delta * s1[j];
    }
    
    assemble_grad_p(E, s2, F, lev);

    valid = solve_del2_u(E, u1, F, imp * v_res, lev);
    strip_bcs_from_residual(E, u1, lev);
    
    assemble_div_u(E, u1, F, lev);
    
    alpha = r1dotz1 / global_pdot(E, s2, F, lev);

    /* r2 = r1 - alpha * div(u1) */
    for(j=1; j<=npno; j++)
      r2[j] = r1[j] - alpha * F[j];

    /* P = P + alpha * s2 */
    for(j=1; j<=npno; j++)
      P[j] += alpha * s2[j];

    /* V = V - alpha * u1 */
    for(j=0; j<neq; j++)
      V[j] -= alpha * u1[j];

    E->monitor.vdotv = global_vdot(E, V, V, lev);
    E->monitor.pdotp = global_pdot(E, P, P, lev);

    assemble_div_u(E, V, z1, lev);
    E->monitor.incompressibility = sqrt(global_div_norm2(E, z1)
					/ (1e-32 + E->monitor.vdotv));

    dpressure = alpha * sqrt(global_pdot(E, s2, s2, lev) / (1.0e-32 + E->monitor.pdotp));
    dvelocity = alpha * sqrt(global_vdot(E, u1, u1, lev) / (1.0e-32 + E->monitor.vdotv));
    /* keep the normalized versions for the message */
    E->monitor.vdotv = sqrt(E->monitor.vdotv/gneq);
    E->monitor.pdotp = sqrt(E->monitor.pdotp/gnpno);

    
    count++;
    
    generate_log_message(count,time0,timea,dvelocity, dpressure,E);

    shuffle = s1;
    s1 = s2;
    s2 = shuffle;

    shuffle = r1;
    r1 = r2;
    r2 = shuffle;

    /* shift <r0, z0> = <r1, z1> */
    r0dotz0 = r1dotz1;
    
    
  }							/* end loop for conjugate gradient   */
  assemble_div_u(E, V, z1, lev);  
  if(problems){
    fprintf(E->fp, "Convergence of velocity solver may affect continuity\n");
    fprintf(E->fp, "Consider running with the `see_convergence=on' option\n");
    fprintf(E->fp, "To evaluate the performance of the current relaxation parameters\n");
    fflush(E->fp);
  }
  
  if(E->control.print_convergence && E->parallel.me == 0){
    fprintf(E->fp, "after (%03d) pressure loops and %g sec for step %d\n", count, CPU_time0() - timea, E->monitor.solution_cycles);
    fprintf(stderr, "after (%03d) pressure loops and %g sec for step %d\n", count, CPU_time0() - timea, E->monitor.solution_cycles);
    fflush(E->fp);
  }

  
  *steps_max = count;
  
  return (residual);
}

void initial_vel_residual(struct All_variables *E,
			  double *V, double *P, double *F,double *u1,
			  double acc)
{
  int neq = E->lmesh.neq;
  int lev = E->mesh.levmax;
  int i, valid;
  
  /* F = F - grad(P) - K*V */
  assemble_grad_p(E, P, u1, lev);
  for(i=0; i<neq; i++)
      F[i] -= u1[i];
  
  assemble_del2_u(E, V, u1, lev, 1);
  for(i=0; i<neq; i++)
    F[i] -= u1[i];

  strip_bcs_from_residual(E, F, lev);

  /* solve K*u1 = F for u1 */
  valid = solve_del2_u(E, u1, F, acc, lev);

  if(!valid && (E->parallel.me==0)) {
    fputs("Warning: solver not converging! 0\n", stderr);
    fputs("Warning: solver not converging! 0\n", E->fp);
  }
  strip_bcs_from_residual(E, u1, lev);

  /* V = V + u1 */
  for(i=0; i < neq; i++)
    V[i] += u1[i];
  
  return;
}

float solve_Ahat_p_fhat(struct All_variables *E, double *V, double *P, 
			double *F, double imp, int *steps_max)
{
	//int i, j, k, ii, count, convergent, valid, problems, lev, lev_low, npno, neq, steps;
	int i, j, count, convergent, valid, problems, lev, npno, neq;
	int gnpno, gneq;

	static int been_here = 0;
	//double *u;
	//static double *p1, *r1, *r0, *r2, *z0, *z1, *s1, *s2, *Ah, *u1;
	static double *r1, *r0, *r2, *z0, *z1, *s1, *s2, *Ah, *u1;
	//double *shuffle, *R;
	double *shuffle;
	double alpha, delta, s2dotAhat, r0dotr0, r1dotz1;
	//double residual, initial_residual, last_residual, res_magnitude, v_res;
	double residual, initial_residual, res_magnitude, v_res;
	char message[500];
	double time0, time;
	static double timea;
	float dpressure, dvelocity;

	//const int dims = E->mesh.nsd;
	//const int n = loc_mat_size[E->mesh.nsd];

	npno = E->lmesh.npno;
	neq = E->lmesh.neq;

	gnpno = E->mesh.npno;
	gneq = E->mesh.neq;

	if(been_here == 0)	/* been here gets incremented further down */
	{
		r0 = (double *)malloc((npno + 1) * sizeof(double));
		r1 = (double *)malloc((npno + 1) * sizeof(double));
		r2 = (double *)malloc((npno + 1) * sizeof(double));
		z0 = (double *)malloc((npno + 1) * sizeof(double));
		z1 = (double *)malloc((npno + 1) * sizeof(double));
		s1 = (double *)malloc((npno + 1) * sizeof(double));
		s2 = (double *)malloc((npno + 1) * sizeof(double));
		Ah = (double *)malloc((neq + 1) * sizeof(double));
		u1 = (double *)malloc((neq + 2) * sizeof(double));

		timea = CPU_time0();
	}

	problems = 0;
	time0 = time = CPU_time0();

	been_here++;

	/* calculate the velocity residual, note there are tricks involved here */

	lev = E->mesh.levmax;

	v_res = sqrt(global_vdot(E, F, F, lev) / gneq);

	if(E->parallel.me == 0)
		fprintf(stderr, "initial residue of momentum equation %g %d\n", v_res, gneq);

	assemble_grad_p(E, P, Ah, lev);
	assemble_del2_u(E, V, u1, lev, 1);

	for(i = 0; i < neq; i++)
		Ah[i] = F[i] - Ah[i] - u1[i];

	strip_bcs_from_residual(E, Ah, lev);

	valid = solve_del2_u(E, u1, Ah, imp * v_res, E->mesh.levmax);
	strip_bcs_from_residual(E, u1, lev);

	for(i = 0; i < neq; i++)
		V[i] += u1[i];

	assemble_div_u(E, V, r1, lev);

	E->monitor.pdotp = sqrt(global_pdot(E, r1, r1, lev) / gnpno);
	residual = initial_residual = 	E->monitor.pdotp;

	E->monitor.vdotv = sqrt(global_vdot(E, V, V, lev) / gneq);

	E->monitor.incompressibility = residual / E->monitor.vdotv;

	if(E->control.print_convergence && E->parallel.me == 0)
		fprintf(stderr, "Loop to reduce pressure residual %g\n", residual);

	count = 0;
	convergent = 0;

	dpressure = 1.0;
	dvelocity = 1.0;

	res_magnitude = residual;

	generate_log_message(count,time0,timea,dvelocity, dpressure,E);


/*   while( (count < *steps_max) && (E->monitor.incompressibility >= E->control.tole_comp || dvelocity >= imp) )  {     
*/ 
	while((count < *steps_max) && (dpressure >= imp || dvelocity >= imp))
	{

		for(j = 1; j <= npno; j++)
			z1[j] = E->BPI[lev][j] * r1[j];

		r1dotz1 = global_pdot(E, r1, z1, lev);

		if((count == 0))
			for(j = 1; j <= npno; j++)
				s2[j] = z1[j];
		else
		{
			r0dotr0 = global_pdot(E, r0, z0, lev);
			assert(r0dotr0 != 0.0 /* Division by zero in head of incompressibility iteration */ );
			delta = r1dotz1 / r0dotr0;
			for(j = 1; j <= npno; j++)
				s2[j] = z1[j] + delta * s1[j];
		}

		assemble_grad_p(E, s2, Ah, lev);

		valid = solve_del2_u(E, u1, Ah, imp * v_res, lev);
		strip_bcs_from_residual(E, u1, lev);

		assemble_div_u(E, u1, Ah, lev);

		s2dotAhat = global_pdot(E, s2, Ah, lev);

		if(valid)
			alpha = r1dotz1 / s2dotAhat;
		else
			alpha = 0.0;

		for(j = 1; j <= npno; j++)
		{
			r2[j] = r1[j] - alpha * Ah[j];
			P[j] += alpha * s2[j];
		}

		for(j = 0; j < neq; j++)
			V[j] -= alpha * u1[j];

		assemble_div_u(E, V, Ah, lev);
		/* this is how it was computed before */
		E->monitor.vdotv = global_vdot(E, V, V, lev);
		E->monitor.pdotp = global_pdot(E, P, P, lev);
		
		E->monitor.incompressibility = sqrt((gneq / gnpno) * (1.0e-32 + global_pdot(E, Ah, Ah, lev) / (1.0e-32 + E->monitor.vdotv)));
		dpressure = alpha * sqrt(global_pdot(E, s2, s2, lev) / (1.0e-32 + E->monitor.pdotp));
		dvelocity = alpha * sqrt(global_vdot(E, u1, u1, lev) / (1.0e-32 + E->monitor.vdotv));
		/* keep the normalized versions for the message */
		E->monitor.vdotv = sqrt(E->monitor.vdotv/gneq);
		E->monitor.pdotp = sqrt(E->monitor.pdotp/gnpno);

		count++;

		generate_log_message(count,time0,timea,dvelocity, dpressure,E);


		shuffle = s1;
		s1 = s2;
		s2 = shuffle;
		shuffle = r0;
		r0 = r1;
		r1 = r2;
		r2 = shuffle;
		shuffle = z0;
		z0 = z1;
		z1 = shuffle;

	}							/* end loop for conjugate gradient   */

	if(problems)
	{
		fprintf(E->fp, "Convergence of velocity solver may affect continuity\n");
		fprintf(E->fp, "Consider running with the `see_convergence=on' option\n");
		fprintf(E->fp, "To evaluate the performance of the current relaxation parameters\n");
		fflush(E->fp);
	}

	if(E->control.print_convergence && E->parallel.me == 0)
	{
		fprintf(E->fp, "after (%03d) pressure loops and %g sec for step %d\n", count, CPU_time0() - timea, E->monitor.solution_cycles);
		/**/ fprintf(stderr, "after (%03d) pressure loops and %g sec for step %d\n", count, CPU_time0() - timea, E->monitor.solution_cycles);
		/**/ fflush(E->fp);
	}


/*
  free((void *) r0);
  free((void *) r1);       
  free((void *) r2);
  free((void *) z0);
  free((void *) z1);
  free((void *) s1);
  free((void *) s2);
  free((void *) Ah);
  free((void *) u1);
*/

	*steps_max = count;

	return (residual);
}



void generate_log_message(int count,double time0,double timea,double dvelocity,double dpressure, struct All_variables *E){

  const int old_version=0;
  char message[500];
  if(E->control.print_convergence && E->parallel.me == 0){
    if(old_version){		/* old version for backward compat */
      if(count == 0)
	/* first message was like this */
	sprintf(message, "AhatP (%03d) after %g sec %g sec with div/v=%.3e for step %d\n", 
		count, CPU_time0() - time0, CPU_time0() - timea, E->monitor.incompressibility, E->monitor.solution_cycles);
      else
	/* other messages like that */
      	sprintf(message,  "AhatP (%03d) after %g sec with div/v=%.3e, dv/v=%.3e & dp/p=%.3e for step %d\n", 
		count, CPU_time0() - time0, E->monitor.incompressibility, dvelocity, dpressure, E->monitor.solution_cycles);
    }else{
      sprintf(message,  "AhatP (%03d) after %8.2f sec with v %12.6e dv/v %9.3e div/v %9.3e  p %12.6e dp/p %9.3e for step %d\n", 
	      count, CPU_time0() - time0, 
	      E->monitor.vdotv,E->monitor.incompressibility, dvelocity, 
	      E->monitor.pdotp,dpressure, E->monitor.solution_cycles);
      

    }
    fprintf(E->fp,"%s",message);//fflush(E->fp);
    fprintf(stderr,"%s",message);
  }
}

/*  ==========================================================================  */

void v_from_vector(struct All_variables *E, float **V, double *F)
{
	int node;
	//int node, d;
	//unsigned int type;

	const int nno = E->lmesh.nno;
	//const int dofs = E->mesh.dof;

	for(node = 1; node <= nno; node++)
	{
	  V[1][node] = F[E->id[node].doff[1]];
	  V[2][node] = F[E->id[node].doff[2]];
	  V[3][node] = F[E->id[node].doff[3]];
	  if(E->node[node] & VBX)
	    V[1][node] = E->VB[1][node];
	  if(E->node[node] & VBY)
	    V[2][node] = E->VB[2][node];
	  if(E->node[node] & VBZ)
	    V[3][node] = E->VB[3][node];
	}
	return;
}
/* for a restart, we might have a velocity solution, so assign this to
   the F vector */
void vector_from_v(struct All_variables *E, double *F, float **V)
{
	int node;

	const int nno = E->lmesh.nno;

	for(node = 1; node <= nno; node++)
	{
	  F[E->id[node].doff[1]] = V[1][node];
	  F[E->id[node].doff[2]] = V[2][node];
	  F[E->id[node].doff[3]] = V[3][node];
	  
	}
	return;
}
