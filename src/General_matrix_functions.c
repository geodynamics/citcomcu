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

/*
 *   Shijie Zhong's version of full multigrid with consistent projection
 *   scheme that speeds up the solution procedure often by X3, compared with
 *   the original one in the code. This full mutligrid solver was implemented
 *   with parallel computing in 1998.
 *   */


#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

#ifdef _UNICOS
#include <fortran.h>
#endif

/* *INDENT-OFF* */
int epsilon[4][4] = {  {0,  0,  0,  0}, /* Levi-Cita epsilon */
                       {0,  1, -1,  1},
                       {0, -1,  1, -1},
                       {0,  1, -1,  1}  };
/* *INDENT-ON* */

static float cost_per_level[MAX_LEVELS];	/* this will accumulate data over the run */
static int total_cycles[MAX_LEVELS];


/*=====================================================================
  Variable dimension matrix allocation  function from numerical recipes
  Note: ANSII consistency requires some additional features !
  =====================================================================  */

double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	double **m;

	/* allocate pointer to rows  */
	m = (double **)malloc((nrow + 1) * sizeof(double *));
	m += 1;
	m -= nrl;

	/*  allocate rows and set the pointers accordingly   */
	m[nrl] = (double *)malloc((nrow * ncol + 1) * sizeof(double));
	m[nrl] += 1;
	m[nrl] -= ncl;

	for(i = nrl + 1; i <= nrh; i++)
		m[i] = m[i - 1] + ncol;

	return (m);
}


float **fmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	float **m;

	/* allocate pointer to rows  */
	m = (float **)malloc((unsigned)((nrow + 1) * sizeof(float *)));
	m += 1;
	m -= nrl;

	/*  allocate rows and set the pointers accordingly   */
	m[nrl] = (float *)malloc((unsigned)((nrow * ncol + 1) * sizeof(float)));
	m[nrl] += 1;
	m[nrl] -= ncl;

	for(i = nrl + 1; i <= nrh; i++)
		m[i] = m[i - 1] + ncol;

	return (m);
}


void dfree_matrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;
	for(i = nrh; i >= nrl; i--)
		free((void *)(m[i] + ncl));
	free((void *)(m + nrl));
	return;
}

void ffree_matrix(float **m, int nrl, int nrh, int ncl, int nch)
{
	int i;
	for(i = nrh; i >= nrl; i--)
		free((void *)(m[i] + ncl));
	free((void *)(m + nrl));
	return;
}

/*=============================================================
  Functions to allocate/remove space for variable sized vector.
  =============================================================  */

double *dvector(int nl, int nh)
{
	double *v;
	v = (double *)malloc((unsigned)(nh - nl + 1) * sizeof(double));
	return (v - nl);
}

float *fvector(int nl, int nh)
{
	float *v;
	v = (float *)malloc((unsigned)(nh - nl + 1) * sizeof(float));
	return (v - nl);
}

void dfree_vector(double *v, int nl, int nh)
{
	free((char *)(v + nl));
}

void ffree_vector(float *v, int nl, int nh)
{
	free((char *)(v + nl));
}

int *sivector(int nl, int nh)
{
	int *v;
	v = (int *)malloc((unsigned)(nh - nl + 1) * sizeof(int));
	return (v - nl);
}

void sifree_vector(int *v, int nl, int nh)
{
	free((char *)(v + nl));
}


double pdot(struct All_variables *E, double *A, double *B, int lev)

{
	double prod;
	int e;
	const int n = E->mesh.NPNO[lev];

	prod = 0.0;
	for(e = 1; e <= n; e++)
		prod += A[e] * B[e];
	return (prod);
}

double pselfdot(struct All_variables *E, double *A)

{
	double prod;
	int e;
	const int n = E->mesh.npno;

	prod = 0.0;
	for(e = 1; e <= n; e++)
		prod += A[e] * A[e];
	return (prod);
}

/* Beware of the alias if A=B on vector machines, use vselfdot instead  */

double vdot(struct All_variables *E, double *A, double *B, int level)
{
	double prod, mprod[1];
	int i, incx = 1;

	char trans = 'N';
	double alpha = 1.0;

	const int n = E->mesh.NEQ[level];

#ifdef blas


#else
	prod = 0.0;
	for(i = 0; i < n; i++)
		prod += A[i] * B[i];
#endif

	return (prod);
}


double vselfdot(struct All_variables *E, double *A, int level)
{
	double prod;
	int i, n;

	n = E->mesh.NEQ[level];

	prod = 0.0;
	for(i = 0; i < n; i++)
		prod += A[i] * A[i];

	return (prod);
}


double vfselfdot(struct All_variables *E, float *A, int level)
{
	float prod;
	int i, n;


	n = E->mesh.NEQ[level];

	prod = 0.0;
	for(i = 0; i < n; i++)
		prod += A[i] * A[i];


	return (prod);
}

float fdot(float *A, float *B, int n1, int n2)
{
	float prod;
	int i;

	prod = 0.0;
	for(i = n1; i <= n2; i++)
		prod += A[i] * B[i];


	return (prod);
}

float fselfdot(float *A, int n1, int n2)
{
	float prod;
	int i;

	prod = 0.0;
	for(i = n1; i <= n2; i++)
		prod += A[i] * A[i];

	return (prod);
}


float dot(struct All_variables *E, float *A, float *B)
{
	float prod = 0.0;
	float domega;
	int e, i, j;

	for(e = 1; e <= E->mesh.nel; e++)
		for(i = 1; i <= vpoints[E->mesh.nsd]; i++)
		{
			j = E->ien[e].node[i];
			domega = E->ECO[E->mesh.levmax][e].area;
			prod += A[j] * B[j] * domega;
		}

	return (prod);
}

float selfdot(struct All_variables *E, float *A)
{
	double prod = 0.0;
	double domega;
	int e, i, j;

	for(e = 1; e <= E->mesh.nel; e++)
		for(i = 1; i <= vpoints[E->mesh.nsd]; i++)
		{
			j = E->ien[e].node[i];
			domega = E->ECO[E->mesh.levmax][e].area;
			prod += A[j] * A[j] * domega;
		}

	return (prod);
}

void dvcopy(double *A, double *B, int a, int b)
{
	int i;

	for(i = a; i <= b; i++)
		A[i] = B[i];

	return;
}

void vcopy(float *A, float *B, int a, int b)
{
	int i;

	for(i = a; i <= b; i++)
		A[i] = B[i];

	return;
}


void vprod(double *R, double *A, double *B, int a, int b)
{
	int i;

	for(i = a; i <= b; i++)
		R[i] = A[i] * B[i];

	return;
}

float fnmax(struct All_variables *E, float *A, int a, int b)
{
	float maxm = -1.0e32;
	int i;

	for(i = a; i <= b; i++)
		if(A[i] > maxm)
			maxm = A[i];

	return (maxm);
}




/*  ===========================================================
    Iterative solver also using multigrid  ........
    ===========================================================  */

int solve_del2_u(struct All_variables *E, double *d0, double *F, double acc, int high_lev, int ic)
{
	static int been_here = 0;
	static int up_heavy, down_heavy, v_steps_high;
	int valid, count, cycles, convergent;
	int i, neq, gneq;

	double initial_time, time;
	double residual, prior_residual, r0;
	static double *D1, *r, *Au;

	neq = E->lmesh.NEQ[high_lev];
	gneq = E->mesh.NEQ[high_lev];

	if(been_here == 0)
	{
		r = (double *)malloc((neq + 10) * sizeof(double));
		D1 = (double *)malloc((neq + 10) * sizeof(double));
		Au = (double *)malloc((neq + 10) * sizeof(double));
		E->control.total_iteration_cycles = 0;
		E->control.total_v_solver_calls = 0;
		for(i = E->mesh.levmin; i <= E->mesh.levmax; i++)
		{
			cost_per_level[i] = 0.0;
			total_cycles[i] = 0;
		}
		v_steps_high = E->control.v_steps_high;
		up_heavy = E->control.up_heavy;
		down_heavy = E->control.down_heavy;
		been_here++;
	}

	for(i = 0; i < neq; i++)
	{
		r[i] = F[i];
		Au[i] = d0[i] = D1[i] = 0.0;
	}

	r0 = residual = sqrt(global_vdot(E, r, r, high_lev) / gneq);

	acc = max(acc, r0 * E->control.accuracy);

	prior_residual = 2 * residual;
	count = 0;
	initial_time = CPU_time0();

	if(!(E->control.NMULTIGRID || E->control.EMULTIGRID))
	{
		cycles = E->control.v_steps_low;
		time = CPU_time0();
		residual = conj_grad(E, D1, r, Au, acc, &cycles, high_lev);
		for(i = 0; i < neq; i++)
		{
			d0[i] += D1[i];
			D1[i] = 0.0;
		}
/*
        count =0;
        do {
            gauss_seidel(E,D1,r,Au,acc,&cycles,high_lev,0);
            for(i=0;i<neq;i++) {
                r[i] = r[i]-Au[i];
                d0[i] += D1[i];
                }
            count ++;
            residual=sqrt(global_vdot(E,r,r,high_lev)/gneq);
            if (E->parallel.me==0) {
                   fprintf(stderr,"resi = %.12e for iter %d acc %.6e\n",residual,count,acc);
                   fprintf(E->fp,"resi = %.12e for iter %d acc %.6e\n",residual,count,acc);
                   }
            } while (residual > acc);
*/


		cost_per_level[high_lev] += CPU_time0() - time;
		total_cycles[high_lev] += cycles;

	}

	else
	{							/* multigrid  */

		cycles = 0;
		count = 0;
		if(E->monitor.solution_cycles % 5 == 0 && E->parallel.me == 0)
		{
			fprintf(E->fp, "resi = %.6e for iter %d acc %.6e\n", residual, count, acc);
		}
		/*    E->control.v_steps_high = v_steps_high;
		 * E->control.down_heavy = down_heavy;
		 * E->control.up_heavy = up_heavy;
		 * if (residual < acc)  {
		 * E->control.v_steps_high = v_steps_high/2;
		 * E->control.down_heavy = down_heavy/2;
		 * E->control.up_heavy = up_heavy/2;
		 * }
		 */
		valid = (residual < acc) ? 0 : 1;
		while(residual > acc)
		{
			residual = multi_grid(E, D1, r, Au, acc, high_lev);

			for(i = 0; i < neq; i++)
			{
				d0[i] += D1[i];
				D1[i] = 0.0;
			}
			count++;
			if(E->monitor.solution_cycles % 5 == 0 && E->parallel.me == 0)
			{
				fprintf(E->fp, "resi = %.6e for iteri %d acc %.6e\n", residual, count, acc);
			}
		}
	}

	if((count > 0) && (residual > r0 * 2.0) || (fabs(residual - prior_residual) < acc * 0.1 && (residual > acc * 10.0)))
		convergent = 0;
	else
	{
		convergent = 1;
		prior_residual = residual;
	}

	if(E->control.print_convergence && E->parallel.me == 0)
	{
		fprintf(E->fp, "%s residual (%03d)(%03d) = %.3e from %.3e to %.3e in %5.2f secs \n", (convergent ? " * " : "!!!"), count, cycles, residual, r0, acc, CPU_time0() - initial_time);
		fflush(E->fp);
	}

	count++;


	if(E->control.verbose)
	{
		printf("Total time for %d loops = %g \n", count, CPU_time0() - initial_time);
		for(i = E->mesh.levmax; i >= E->mesh.levmin; i--)
		{
			printf("Level %d, total time = %g\n", i, cost_per_level[i]);
			printf("     total cycles = %d (%g)\n", total_cycles[i], cost_per_level[i] / (1.0e-5 + (float)total_cycles[i]));
		}
		printf("projection time = %g\n", E->monitor.cpu_time_on_mg_maps);
		printf("Del sq u solved to accuracy of  %g \n", residual);
	}

	E->control.total_iteration_cycles += count;
	E->control.total_v_solver_calls += 1;

	return (valid);
}

double multi_grid(
	struct All_variables *E,
	double *d1,
	double *F,
	double *Au,
	double acc,
	int hl						/* higher level of two */
)
{
	double residual, AudotAu;
	int lev, dlev, ulev, i, j, Vn, Vnmax, ic, cycles;
	double residuaa, alpha, beta;

	FILE *fp;
	char filename[1000];

	const int levmin = E->mesh.levmin;
	const int levmax = E->mesh.levmax;

	double time;

	static int been_here = 0;
	static double *res[MAX_LEVELS], *rhs[MAX_LEVELS], *AU[MAX_LEVELS];
	static double *vel[MAX_LEVELS], *fl[MAX_LEVELS], *del_vel[MAX_LEVELS];
	if(0 == been_here)
	{
		for(i = E->mesh.levmin; i <= E->mesh.levmax; i++)
		{
			vel[i] = (double *)malloc((E->lmesh.NEQ[i] + 2) * sizeof(double));
			res[i] = (double *)malloc((E->lmesh.NEQ[i] + 2) * sizeof(double));
			rhs[i] = (double *)malloc((E->lmesh.NEQ[i] + 2) * sizeof(double));
			fl[i] = (double *)malloc((E->lmesh.NEQ[i] + 2) * sizeof(double));
			del_vel[i] = (double *)malloc((E->lmesh.NEQ[i] + 2) * sizeof(double));
			AU[i] = (double *)malloc((E->lmesh.NEQ[i] + 2) * sizeof(double));
		}
	}
	been_here = 1;

	Vnmax = E->control.mg_cycle;

	for(j = 0; j < E->lmesh.NEQ[levmax]; j++)
		fl[levmax][j] = F[j];

	/* Project residual onto all the lower levels */

	for(lev = levmax; lev > levmin; lev--)
	{
		project_vector(E, lev, fl[lev], fl[lev - 1], 1);
		strip_bcs_from_residual(E, fl[lev - 1], lev - 1);
	}

	/* Solve for the lowest level */

	cycles = E->control.v_steps_low;
/*    (void) conj_grad(E,vel[levmin],fl[levmin],AU[levmin],acc*0.001,&cycles,levmin); 
 */ gauss_seidel(E, vel[levmin], fl[levmin], AU[levmin], acc * 0.01, &cycles, levmin, 0);

	for(lev = levmin + 1; lev <= levmax; lev++)
	{
		time = CPU_time0();

		/* Utilize coarse solution and smooth at this level */
		interp_vector(E, lev - 1, vel[lev - 1], vel[lev]);
		strip_bcs_from_residual(E, vel[lev], lev);

		for(j = 0; j < E->lmesh.NEQ[lev]; j++)
			rhs[lev][j] = fl[lev][j];

		for(Vn = 1; Vn <= Vnmax; Vn++)
		{
			/*    Downward stoke of the V    */
			for(dlev = lev; dlev >= levmin + 1; dlev--)
			{


				/* Pre-smoothing  */
				cycles = ((dlev == levmax) ? E->control.v_steps_high : E->control.down_heavy);
				ic = ((dlev == lev) ? 1 : 0);
				gauss_seidel(E, vel[dlev], rhs[dlev], AU[dlev], 0.01, &cycles, dlev, ic);

				/* Update residual  */
				for(i = 0; i < E->lmesh.NEQ[dlev]; i++)
					res[dlev][i] = rhs[dlev][i] - AU[dlev][i];

				/* Project residual to the lower levels */
				project_vector(E, dlev, res[dlev], rhs[dlev - 1], 1);
				strip_bcs_from_residual(E, rhs[dlev - 1], dlev - 1);

			}


			/*    Bottom of the V    */
			cycles = E->control.v_steps_low;
/*      (void) conj_grad(E,vel[levmin],rhs[levmin],AU[levmin],acc*0.001,&cycles,levmin); 
*/ gauss_seidel(E, vel[levmin], rhs[levmin], AU[levmin], acc * 0.01, &cycles, levmin, 0);


			/*    Upward stoke of the V    */

			for(ulev = levmin + 1; ulev <= lev; ulev++)
			{
				cycles = ((ulev == levmax) ? E->control.v_steps_high : E->control.up_heavy);

				interp_vector(E, ulev - 1, vel[ulev - 1], del_vel[ulev]);
				strip_bcs_from_residual(E, del_vel[ulev], ulev);
				gauss_seidel(E, del_vel[ulev], res[ulev], AU[ulev], 0.01, &cycles, ulev, 1);

				AudotAu = global_vdot(E, AU[ulev], AU[ulev], ulev);
				alpha = global_vdot(E, AU[ulev], res[ulev], ulev) / AudotAu;

				for(i = 0; i < E->lmesh.NEQ[ulev]; i++)
					vel[ulev][i] += alpha * del_vel[ulev][i];

				if(ulev == levmax)
					for(i = 0; i < E->lmesh.NEQ[ulev]; i++)
					{
						res[ulev][i] -= alpha * AU[ulev][i];
					}


			}
		}
	}

	for(j = 0; j < E->lmesh.NEQ[levmax]; j++)
	{
		F[j] = res[levmax][j];
		d1[j] = vel[levmax][j];
	}

	residual = sqrt(global_vdot(E, F, F, levmax) / E->mesh.NEQ[levmax]);


	return (residual);
}

/*  ===========================================================
    Conjugate gradient relaxation for the matrix equation Kd = f
    Returns the residual reduction after itn iterations ... 
    ===========================================================  */


double conj_grad(struct All_variables *E, double *d0, double *F, double *Au, double acc, int *cycles, int level)
{
	static double *r0, *r1, *r2;
	static double *z0, *z1, *z2;
	static double *p1, *p2;
	static double *Ap;
	static double *BI;
	static double *shuffle;
	static int been_here = 0;

	int count, i, steps;
	double residual;
	double alpha, beta, dotprod, dotr1z1, dotr0z0;

	double time;

	const int mem_lev = E->mesh.levmax;
	const int high_neq = E->lmesh.NEQ[level];

	steps = *cycles;

	if(0 == been_here)			/* only used at low level (even if low=high) */
	{
		r0 = (double *)malloc((1 + E->lmesh.NEQ[mem_lev]) * sizeof(double));
		r1 = (double *)malloc((1 + E->lmesh.NEQ[mem_lev]) * sizeof(double));
		r2 = (double *)malloc((1 + E->lmesh.NEQ[mem_lev]) * sizeof(double));
		z0 = (double *)malloc((1 + E->lmesh.NEQ[mem_lev]) * sizeof(double));
		z1 = (double *)malloc((1 + E->lmesh.NEQ[mem_lev]) * sizeof(double));
		p1 = (double *)malloc((1 + E->lmesh.NEQ[mem_lev]) * sizeof(double));
		p2 = (double *)malloc((1 + E->lmesh.NEQ[mem_lev]) * sizeof(double));
		Ap = (double *)malloc((1 + E->lmesh.NEQ[mem_lev]) * sizeof(double));
		been_here++;
	}

	for(i = 0; i < high_neq; i++)
	{
		r1[i] = F[i];
		d0[i] = 0.0;			/* this may get used in multigrid ? SJZ   */
	}

	residual = sqrt(global_vdot(E, r1, r1, level) / E->mesh.NEQ[level]);

	assert(residual != 0.0 /* initial residual for CG = 0.0 */ );
	count = 0;

	while(((residual > acc) && (count < steps)) || count == 0)
	{

		for(i = 0; i < high_neq; i++)
			z1[i] = E->BI[level][i] * r1[i];

		dotr1z1 = global_vdot(E, r1, z1, level);

		if(0 == count)
			for(i = 0; i < high_neq; i++)
				p2[i] = z1[i];
		else
		{
			assert(dotr0z0 != 0.0 /* in head of conj_grad */ );
			beta = dotr1z1 / dotr0z0;
			for(i = 0; i < high_neq; i++)
				p2[i] = z1[i] + beta * p1[i];
		}

		dotr0z0 = dotr1z1;

		assemble_del2_u(E, p2, Ap, level, 1);

		dotprod = global_vdot(E, p2, Ap, level);

		if(0.0 == dotprod)
			alpha = 1.0e-3;
		else
			alpha = dotr1z1 / dotprod;

		for(i = 0; i < high_neq; i++)
		{
			d0[i] += alpha * p2[i];
			r2[i] = r1[i] - alpha * Ap[i];
		}

		residual = sqrt(global_vdot(E, r2, r2, level) / E->mesh.NEQ[level]);

		shuffle = r0;
		r0 = r1;
		r1 = r2;
		r2 = shuffle;
		shuffle = z0;
		z0 = z1;
		z1 = shuffle;
		shuffle = p1;
		p1 = p2;
		p2 = shuffle;

		count++;
		/* end of while-loop */

	}

	*cycles = count;

	strip_bcs_from_residual(E, d0, level);


	return (residual);
}


/* ===========================================================================

 */


void jacobi(struct All_variables *E, double *d0, double *F, double *Ad, double acc, int *cycles, int level, int guess)
{
	static double *r1;
	static int been_here = 0;

	int count, steps;
	int i, j, k, eqn1, eqn2, eqn3;
	int *C;
	double U1, U2, U3;

	higher_precision *B1, *B2, *B3;

	const int dims = E->mesh.nsd, dofs = E->mesh.dof;
	const int neq = E->lmesh.NEQ[level];
	const int max_eqn = max_eqn_interaction[dims];

	if(0 == been_here)
	{
		r1 = (double *)malloc((neq + 1) * sizeof(double));
		been_here++;
	}


	steps = *cycles;
	count = 0;

	if(guess)
	{
		assemble_del2_u(E, d0, Ad, level, 1);
		for(i = 0; i < neq; i++)
			r1[i] = F[i] - Ad[i];
	}
	else
		for(i = 0; i < neq; i++)
		{
			r1[i] = F[i];
			d0[i] = 0.0;
		}

	while(count <= steps)
	{
		for(i = 0; i < neq; i++)
		{
			d0[i] += r1[i] * E->BI[level][i];
			Ad[i] = 0.0;
		}

		for(i = 1; i <= E->lmesh.NNO[level]; i++)
		{
			eqn1 = E->ID[level][i].doff[1];
			eqn2 = E->ID[level][i].doff[2];
			U1 = d0[eqn1];
			U2 = d0[eqn2];

			C = E->Node_map[level] + (i - 1) * max_eqn;
			B1 = E->Eqn_k1[level] + (i - 1) * max_eqn;
			B2 = E->Eqn_k2[level] + (i - 1) * max_eqn;
			if(3 == dims)
			{
				U3 = d0[eqn3];
				eqn3 = E->ID[level][i].doff[3];
				B3 = E->Eqn_k3[level] + (i - 1) * max_eqn;
				for(j = 3; j < max_eqn; j++)
				{
					Ad[eqn1] += B1[j] * d0[C[j]];
					Ad[eqn2] += B2[j] * d0[C[j]];
					Ad[eqn3] += B3[j] * d0[C[j]];
				}

				for(j = 0; j < max_eqn; j++)
					Ad[C[j]] += B1[j] * U1 + B2[j] * U2 + B3[j] * U3;
			}
			else
			{
				for(j = 2; j < max_eqn; j++)
				{
					Ad[eqn1] += B1[j] * d0[C[j]];
					Ad[eqn2] += B2[j] * d0[C[j]];
				}

				for(j = 0; j < max_eqn; j++)
					Ad[C[j]] += B1[j] * U1 + B2[j] * U2;
			}
		}

		exchange_id_d20(E, Ad, level);

		for(i = 0; i < neq; i++)
		{
			r1[i] = F[i] - Ad[i];
			Ad[i] = 0.0;
		}

		count++;

	}

	return;
}

/* ========================================================================================
   An element by element version of the gauss-seidel routine. Initially this is a test 
   platform, we want to know if it handles discontinuities any better than the node/equation
   versions

   Be careful: no-parallel
   =========================================================================================*/

void element_gauss_seidel(struct All_variables *E, double *d0, double *F, double *Ad, double acc, int *cycles, int level, int guess)
{
	int count, i, j, k, l, m, ns, nc, d, steps, loc;
	int p1, p2, p3, q1, q2, q3;
	int e, eq, node, node1;
	int element, eqn1, eqn2, eqn3, eqn11, eqn12, eqn13;

	double U1[24], AD1[24], F1[24];
	double w1, w2, w3;
	double w11, w12, w13;
	double w[24];
	static double *Ad0, *dd, *elt_k;
	static int *vis, been_here = 0;

	const int dims = E->mesh.nsd;
	const int ends = enodes[dims];
	const int n = loc_mat_size[E->mesh.nsd];
	const int neq = E->lmesh.NEQ[level];
	const int nel = E->lmesh.NEL[level];
	const int nno = E->lmesh.NNO[level];


	steps = *cycles;

	if(0 == been_here)
	{
		dd = (double *)malloc((neq + 1) * sizeof(double));
		vis = (int *)malloc((nno + 1) * sizeof(int));
		elt_k = (double *)malloc((24 * 24) * sizeof(double));
		been_here++;
	}

	if(guess)
	{
		e_assemble_del2_u(E, d0, Ad, level, 1);
	}
	else
	{
		for(i = 0; i < neq; i++)
			Ad[i] = d0[i] = 0.0;
	}

	count = 0;
	while(count <= steps)
	{
		for(i = 1; i <= nno; i++)
			vis[i] = 0;

		for(e = 1; e <= nel; e++)
		{

			elt_k = E->elt_k[level][e].k;

			for(i = 1; i <= ends; i++)
			{
				node = E->IEN[level][e].node[i];
				p1 = (i - 1) * dims;
				w[p1] = w[p1 + 1] = 1.0;
				if(dims == 3)
					w[p1 + 2] = 1.0;
				if(E->NODE[level][node] & VBX)
					w[p1] = 0.0;
				if(E->NODE[level][node] & VBY)
					w[p1 + 1] = 0.0;
				if(E->NODE[level][node] & VBZ)
					w[p1 + 2] = 0.0;

			}


			for(i = 1; i <= ends; i++)
			{
				node = E->IEN[level][e].node[i];
				if(!vis[node])
					continue;

				eqn1 = E->ID[level][node].doff[1];
				eqn2 = E->ID[level][node].doff[2];
				eqn3 = E->ID[level][node].doff[3];
				p1 = (i - 1) * dims * n;
				p2 = p1 + n;
				p3 = p2 + n;


				/* update Au */
				for(j = 1; j <= ends; j++)
				{
					node1 = E->IEN[level][e].node[j];

					if(3 == dims)
					{
						eqn11 = E->ID[level][node1].doff[1];
						eqn12 = E->ID[level][node1].doff[2];
						eqn13 = E->ID[level][node1].doff[3];
						q1 = (j - 1) * 3;

						Ad[eqn11] += w[q1] * (elt_k[p1 + q1] * dd[eqn1] + elt_k[p2 + q1] * dd[eqn2] + elt_k[p3 + q1] * dd[eqn3]);
						Ad[eqn12] += w[q1 + 1] * (elt_k[p1 + q1 + 1] * dd[eqn1] + elt_k[p2 + q1 + 1] * dd[eqn2] + elt_k[p3 + q1 + 1] * dd[eqn3]);
						Ad[eqn13] += w[q1 + 2] * (elt_k[p1 + q1 + 2] * dd[eqn1] + elt_k[p2 + q1 + 2] * dd[eqn2] + elt_k[p3 + q1 + 2] * dd[eqn3]);
					}

					else
					{
						eqn11 = E->ID[level][node1].doff[1];
						eqn12 = E->ID[level][node1].doff[2];
						q1 = (j - 1) * 2;

						Ad[eqn11] += w[q1] * (elt_k[p1 + q1] * dd[eqn1] + elt_k[p2 + q1] * dd[eqn2]);
						Ad[eqn12] += w[q1 + 1] * (elt_k[p1 + q1 + 1] * dd[eqn1] + elt_k[p2 + q1 + 1] * dd[eqn2]);
					}
				}
			}


			for(i = 1; i <= ends; i++)
			{
				node = E->IEN[level][e].node[i];
				if(vis[node])
					continue;

				eqn1 = E->ID[level][node].doff[1];
				eqn2 = E->ID[level][node].doff[2];
				eqn3 = E->ID[level][node].doff[3];
				p1 = (i - 1) * dims * n;
				p2 = p1 + n;
				p3 = p2 + n;

				/* update dd, d0 */
				d0[eqn1] += (dd[eqn1] = w[(i - 1) * dims] * (F[eqn1] - Ad[eqn1]) * E->BI[level][eqn1]);
				d0[eqn2] += (dd[eqn2] = w[(i - 1) * dims + 1] * (F[eqn2] - Ad[eqn2]) * E->BI[level][eqn2]);
				if(3 == dims)
				{
					d0[eqn3] += (dd[eqn3] = w[(i - 1) * dims + 2] * (F[eqn3] - Ad[eqn3]) * E->BI[level][eqn3]);
				}

				vis[node] = 1;

				/* update Au */
				for(j = 1; j <= ends; j++)
				{
					node1 = E->IEN[level][e].node[j];

					if(3 == dims)
					{
						eqn11 = E->ID[level][node1].doff[1];
						eqn12 = E->ID[level][node1].doff[2];
						eqn13 = E->ID[level][node1].doff[3];
						q1 = (j - 1) * 3;
						q2 = q1 + 1;
						q3 = q1 + 2;

						Ad[eqn11] += w[q1] * (elt_k[p1 + q1] * dd[eqn1] + elt_k[p2 + q1] * dd[eqn2] + elt_k[p3 + q1] * dd[eqn3]);
						Ad[eqn12] += w[q2] * (elt_k[p1 + q2] * dd[eqn1] + elt_k[p2 + q2] * dd[eqn2] + elt_k[p3 + q2] * dd[eqn3]);
						Ad[eqn13] += w[q3] * (elt_k[p1 + q3] * dd[eqn1] + elt_k[p2 + q3] * dd[eqn2] + elt_k[p3 + q3] * dd[eqn3]);
					}
					else
					{
						eqn11 = E->ID[level][node1].doff[1];
						eqn12 = E->ID[level][node1].doff[2];
						q1 = (j - 1) * 2;
						q2 = q1 + 1;
						Ad[eqn11] += w[q1] * (elt_k[p1 + q1] * dd[eqn1] + elt_k[p2 + q1] * dd[eqn2]);
						Ad[eqn12] += w[q2] * (elt_k[p1 + q2] * dd[eqn1] + elt_k[p2 + q2] * dd[eqn2]);
					}

				}
			}

		}
		/* completed cycle */

		count++;

	}

	return;
}

void gauss_seidel1(struct All_variables *E, double *d0, double *F, double *Ad, double acc, int *cycles, int level, int guess)
{

	int count, i, j, k, l, m, ns, steps;
	int *C;
	int eqn1, eqn2, eqn3;

	double UU, U1, U2, U3;
	static double zeroo = 0.0;

	higher_precision node_k[4][81];
	higher_precision *temp1, *temp, *B1, *B2, *B3;


	const int dims = E->mesh.nsd;
	const int ends = enodes[dims];
	const int n = loc_mat_size[E->mesh.nsd];
	const int neq = E->lmesh.NEQ[level];
	const int nno = E->lmesh.NNO[level];
	const int nox = E->lmesh.NOX[level];
	const int noz = E->lmesh.NOY[level];
	const int noy = E->lmesh.NOZ[level];
	const int max_eqn = max_eqn_interaction[dims];

	steps = *cycles;

	if(guess)
	{
		d0[neq] = 0.0;
		n_assemble_del2_u(E, d0, Ad, level, 1);
	}
	else
		for(i = 0; i < neq; i++)
		{
			d0[i] = Ad[i] = zeroo;
		}


	temp = (higher_precision *) malloc((neq + 1) * sizeof(higher_precision));

	count = 0;

	while(count < steps)
	{
		for(i = 0; i <= neq; i++)
			temp[i] = zeroo;
		for(i = 1; i <= nno; i++)
		{
			C = E->Node_map[level] + (i - 1) * max_eqn;
			B1 = E->Eqn_k1[level] + (i - 1) * max_eqn;
			B2 = E->Eqn_k2[level] + (i - 1) * max_eqn;

			eqn1 = E->ID[level][i].doff[1];
			temp[eqn1] = (F[eqn1] - Ad[eqn1]) * E->BI[level][eqn1];
			d0[eqn1] += temp[eqn1];

			eqn2 = E->ID[level][i].doff[2];
			temp[eqn2] = (F[eqn2] - Ad[eqn2]) * E->BI[level][eqn2];
			d0[eqn2] += temp[eqn2];

			for(j = 2; j < max_eqn; j++)
			{
				UU = temp[C[j]];
				Ad[eqn1] += B1[j] * UU;
				Ad[eqn2] += B2[j] * UU;
			}

			for(j = 0; j < max_eqn; j++)
				Ad[C[j]] += B1[j] * temp[eqn1] + B2[j] * temp[eqn2];

		}
		count++;
	}


	*cycles = count;
	return;
}

/* ============================================================================
   Multigrid Gauss-Seidel relaxation scheme which requires the storage of local
   information, otherwise some other method is required. NOTE this is a bit worse
   than real gauss-seidel because it relaxes all the equations for a node at one
   time (Jacobi at a node). It does the job though.
   ============================================================================ */

void gauss_seidel(struct All_variables *E, double *d0, double *F, double *Ad, double acc, int *cycles, int level, int guess)
{

	int count, i, j, k, l, m, ns, steps;
	int *C;
	int eqn1, eqn2, eqn3;

	double residual, *r, UU, U1, U2, U3;
	static double zeroo = 0.0;
	static int been_here = 0;
	static higher_precision *temp1, *temp;

	higher_precision node_k[4][81];
	higher_precision *B1, *B2, *B3;


	const int dims = E->mesh.nsd;
	const int ends = enodes[dims];
	const int n = loc_mat_size[E->mesh.nsd];
	const int neq = E->lmesh.NEQ[level];
	const int nno = E->lmesh.NNO[level];
	const int nox = E->lmesh.NOX[level];
	const int noz = E->lmesh.NOY[level];
	const int noy = E->lmesh.NOZ[level];
	const int max_eqn = max_eqn_interaction[dims];

	steps = *cycles;

	if(guess)
	{
		d0[neq] = 0.0;
		n_assemble_del2_u(E, d0, Ad, level, 1);
	}
	else
		for(i = 0; i < neq; i++)
		{
			d0[i] = Ad[i] = zeroo;
		}

	count = 0;

	if(been_here == 0)
	{
		temp = (higher_precision *) malloc((E->lmesh.NEQ[E->mesh.levmax] + 1) * sizeof(higher_precision));
		temp1 = (higher_precision *) malloc((E->lmesh.NEQ[E->mesh.levmax] + 1) * sizeof(higher_precision));
		been_here++;
	}

	while(count < steps)
	{
		for(i = 0; i <= neq; i++)
			temp[i] = zeroo;

		if(dims == 3)
		{

			for(i = 1; i <= nno; i++)
				if(E->NODE[level][i] & OFFSIDE)
				{
					eqn1 = E->ID[level][i].doff[1];
					eqn2 = E->ID[level][i].doff[2];
					eqn3 = E->ID[level][i].doff[3];
					temp[eqn1] = (F[eqn1] - Ad[eqn1]) * E->BI[level][eqn1];
					temp[eqn2] = (F[eqn2] - Ad[eqn2]) * E->BI[level][eqn2];
					temp[eqn3] = (F[eqn3] - Ad[eqn3]) * E->BI[level][eqn3];
					temp1[eqn1] = Ad[eqn1];
					temp1[eqn2] = Ad[eqn2];
					temp1[eqn3] = Ad[eqn3];
				}
			for(i = 1; i <= nno; i++)
			{
				eqn1 = E->ID[level][i].doff[1];
				eqn2 = E->ID[level][i].doff[2];
				eqn3 = E->ID[level][i].doff[3];
				C = E->Node_map[level] + (i - 1) * max_eqn;
				B1 = E->Eqn_k1[level] + (i - 1) * max_eqn;
				B2 = E->Eqn_k2[level] + (i - 1) * max_eqn;
				B3 = E->Eqn_k3[level] + (i - 1) * max_eqn;

				for(j = 3; j < max_eqn; j++)
				{
					UU = temp[C[j]];
					Ad[eqn1] += B1[j] * UU;
					Ad[eqn2] += B2[j] * UU;
					Ad[eqn3] += B3[j] * UU;
				}
				if(!(E->NODE[level][i] & OFFSIDE))
				{
					temp[eqn1] = (F[eqn1] - Ad[eqn1]) * E->BI[level][eqn1];
					temp[eqn2] = (F[eqn2] - Ad[eqn2]) * E->BI[level][eqn2];
					temp[eqn3] = (F[eqn3] - Ad[eqn3]) * E->BI[level][eqn3];
				}
				for(j = 0; j < max_eqn; j++)
					Ad[C[j]] += B1[j] * temp[eqn1] + B2[j] * temp[eqn2] + B3[j] * temp[eqn3];

				d0[eqn1] += temp[eqn1];
				d0[eqn2] += temp[eqn2];
				d0[eqn3] += temp[eqn3];
			}

			for(i = 1; i <= nno; i++)
				if(E->NODE[level][i] & OFFSIDE)
				{
					eqn1 = E->ID[level][i].doff[1];
					eqn2 = E->ID[level][i].doff[2];
					eqn3 = E->ID[level][i].doff[3];
					Ad[eqn1] -= temp1[eqn1];
					Ad[eqn2] -= temp1[eqn2];
					Ad[eqn3] -= temp1[eqn3];
				}

			exchange_id_d20(E, Ad, level);

			for(i = 1; i <= nno; i++)
				if(E->NODE[level][i] & OFFSIDE)
				{
					eqn1 = E->ID[level][i].doff[1];
					eqn2 = E->ID[level][i].doff[2];
					eqn3 = E->ID[level][i].doff[3];
					Ad[eqn1] += temp1[eqn1];
					Ad[eqn2] += temp1[eqn2];
					Ad[eqn3] += temp1[eqn3];
				}

		}

		else if(2 == dims)
		{

			for(i = 1; i <= nno; i++)
				if(E->NODE[level][i] & OFFSIDE)
				{
					eqn1 = E->ID[level][i].doff[1];
					eqn2 = E->ID[level][i].doff[2];
					temp[eqn1] = (F[eqn1] - Ad[eqn1]) * E->BI[level][eqn1];
					temp[eqn2] = (F[eqn2] - Ad[eqn2]) * E->BI[level][eqn2];
					temp1[eqn1] = Ad[eqn1];
					temp1[eqn2] = Ad[eqn2];
				}
			for(i = 1; i <= nno; i++)
			{
				eqn1 = E->ID[level][i].doff[1];
				eqn2 = E->ID[level][i].doff[2];
				C = E->Node_map[level] + (i - 1) * max_eqn;
				B1 = E->Eqn_k1[level] + (i - 1) * max_eqn;
				B2 = E->Eqn_k2[level] + (i - 1) * max_eqn;

				for(j = 2; j < max_eqn; j++)
				{
					UU = temp[C[j]];
					Ad[eqn1] += B1[j] * UU;
					Ad[eqn2] += B2[j] * UU;
				}
				if(!(E->NODE[level][i] & OFFSIDE))
				{
					temp[eqn1] = (F[eqn1] - Ad[eqn1]) * E->BI[level][eqn1];
					temp[eqn2] = (F[eqn2] - Ad[eqn2]) * E->BI[level][eqn2];
				}
				for(j = 0; j < max_eqn; j++)
					Ad[C[j]] += B1[j] * temp[eqn1] + B2[j] * temp[eqn2];

				d0[eqn1] += temp[eqn1];
				d0[eqn2] += temp[eqn2];
			}

			for(i = 1; i <= nno; i++)
				if(E->NODE[level][i] & OFFSIDE)
				{
					eqn1 = E->ID[level][i].doff[1];
					eqn2 = E->ID[level][i].doff[2];
					Ad[eqn1] -= temp1[eqn1];
					Ad[eqn2] -= temp1[eqn2];
				}

			exchange_id_d20(E, Ad, level);

			for(i = 1; i <= nno; i++)
				if(E->NODE[level][i] & OFFSIDE)
				{
					eqn1 = E->ID[level][i].doff[1];
					eqn2 = E->ID[level][i].doff[2];
					Ad[eqn1] += temp1[eqn1];
					Ad[eqn2] += temp1[eqn2];
				}
		}

		count++;
	}



	*cycles = count;

/*
    free((void *)temp);
    free((void *)temp1);
*/
	return;

}

void print_elt_k(struct All_variables *E, double a[24 * 24])
{
	int l, ll, n;

	printf("elt k is ...\n");


	n = loc_mat_size[E->mesh.nsd];

	for(l = 0; l < n; l++)
	{
		fprintf(stderr, "\n");
		fflush(stderr);
		for(ll = 0; ll < n; ll++)
		{
			fprintf(stderr, "%s%.3e ", a[ll * n + l] >= 0.0 ? "+" : "", a[ll * n + l]);
			fflush(stderr);
		}
	}
	fprintf(stderr, "\n");
	fflush(stderr);

	return;
}


double cofactor(double A[4][4], int i, int j, int n)
{
	int k, l, p, q;
	static int been_here = 0;
	double B[4][4];				/* because of recursive behaviour of det/cofac, need to use
								 * new copy of B at each 'n' level of this routine */

	if(n > 3)
		printf("Error, no cofactors for matrix more than 3x3\n");

	p = q = 1;

	for(k = 1; k <= n; k++)
	{
		if(k == i)
			continue;
		for(l = 1; l <= n; l++)
		{
			if(l == j)
				continue;
			B[p][q] = A[k][l];
			q++;
		}
		q = 1;
		p++;
	}


	return (epsilon[i][j] * determinant(B, n - 1));


}


/* Fast (conditional) determinant for 3x3 or 2x2 ... otherwise calls general routine */

double determinant(double A[4][4], int n)
{
	switch (n)
	{
	case 1:
		return (A[1][1]);
		break;
	case 2:
		return (A[1][1] * A[2][2] - A[1][2] * A[2][1]);
		break;
	case 3:
		return (A[1][1] * (A[2][2] * A[3][3] - A[2][3] * A[3][2]) - A[1][2] * (A[2][1] * A[3][3] - A[2][3] * A[3][1]) + A[1][3] * (A[2][1] * A[3][2] - A[2][2] * A[3][1]));
		break;
	default:
		return (1);
/*		return(gen_determinant(A,n)); */
	}
}


/* recursive function to determine matrix determinant */
/* TODO: need to modify cofactor() in order to use it.
 * Note that, currently, this function is never called */
/*double gen_determinant(A, n)
	double **A;
	int n;
{
	double det;
	double cofactor();

	int i;

	if(n == 1)
		return (A[1][1]);		// need a way to break the recursion

	det = 0.0;
	for(i = 1; i <= n; i++)
		det += A[1][i] * cofactor(A, 1, i, n);

	return (det);
}
*/


float area_of_4node(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4)
{
	float area;

	area = fabs(0.5 * (x1 * (y2 - y4) + x2 * (y4 - y1) + x4 * (y1 - y2))) + fabs(0.5 * (x2 * (y3 - y4) + x3 * (y4 - y2) + x4 * (y2 - y3)));

	return area;
}

 /*=====================================*/
double modified_plgndr_a(int l, int m, double t)
{
	int i, ll;
	double x, fact1, fact2, fact, pll, pmm, pmmp1, somx2, plgndr;
	const double three = 3.0;
	const double two = 2.0;
	const double one = 1.0;

	x = cos(t);
	pmm = one;
	if(m > 0)
	{
		somx2 = sqrt((one - x) * (one + x));
		fact1 = three;
		fact2 = two;
		for(i = 1; i <= m; i++)
		{
			fact = sqrt(fact1 / fact2);
			pmm = -pmm * fact * somx2;
			fact1 += two;
			fact2 += two;
		}
	}

	if(l == m)
		plgndr = pmm;
	else
	{
		pmmp1 = x * sqrt(two * m + three) * pmm;
		if(l == m + 1)
			plgndr = pmmp1;
		else
		{
			for(ll = m + 2; ll <= l; ll++)
			{
				fact1 = sqrt((4.0 * ll * ll - one) * (double)(ll - m) / (double)(ll + m));
				fact2 = sqrt((2.0 * ll + one) * (ll - m) * (ll + m - one) * (ll - m - one) / (double)((two * ll - three) * (ll + m)));
				pll = (x * fact1 * pmmp1 - fact2 * pmm) / (ll - m);
				pmm = pmmp1;
				pmmp1 = pll;
			}
			plgndr = pll;
		}
	}

	plgndr /= sqrt(4.0 * M_PI);

	if(m != 0)
		plgndr *= sqrt(two);

	return plgndr;
}

 /* ===================================  */
double sqrt_multis(int jj, int ii)
{
	int i;
	double sqrt_multisa;

	sqrt_multisa = 1.0;
	if(jj > ii)
		for(i = jj; i > ii; i--)
			sqrt_multisa *= 1.0 / sqrt((double)i);

	return sqrt_multisa;
}



 /* ===================================  */
double multis(int ii)
{
	int i;
	double multisa;

	multisa = 1.0;
	if(ii)
		for(i = 2; i <= ii; i++)
			multisa *= (double)i;

	return multisa;
}


 /* ===================================  */
int int_multis(int ii)
{
	int i, multisa;

	multisa = 1;
	if(ii)
		for(i = 2; i <= ii; i++)
			multisa *= i;

	return multisa;
}

/* =====================================*/
double plgndr_a(int l, int m, double t)
{
	int i, ll;
	double x, fact, pll, pmm, pmmp1, somx2, plgndr;
	const double two = 2.0;
	const double one = 1.0;

	x = cos(t);
	pmm = one;
	if(m > 0)
	{
		somx2 = sqrt((one - x) * (one + x));
		fact = one;
		for(i = 1; i <= m; i++)
		{
			pmm = -pmm * fact * somx2;
			fact = fact + two;
		}
	}

	if(l == m)
		plgndr = pmm;
	else
	{
		pmmp1 = x * (2 * m + 1) * pmm;
		if(l == m + 1)
			plgndr = pmmp1;
		else
		{
			for(ll = m + 2; ll <= l; ll++)
			{
				pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
				pmm = pmmp1;
				pmmp1 = pll;
			}
			plgndr = pll;
		}
	}

	return plgndr;
}

/* =====================================*/
double sphere_h(int l, int m, double t, double f, int ic)
{

	double sphere_hamonics;

	sphere_hamonics = 0.0;
	if(ic == 0)
		sphere_hamonics = cos(m * f) * plgndr_a(l, m, t);
	else if(m)
		sphere_hamonics = sin(m * f) * plgndr_a(l, m, t);

	return sphere_hamonics;
}
