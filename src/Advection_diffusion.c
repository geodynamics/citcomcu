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

/*   Functions which solve the heat transport equations using Petrov-Galerkin
     streamline-upwind methods. The process is basically as described in Alex
     Brooks PhD thesis (Caltech) which refers back to Hughes, Liu and Brooks.  */

/* Extended Boussinesq approximation is implemented by SZ in 2003 */
/* Regional spherical geometry is implemented by SZ in 2003 */

#include <mpi.h>

#include <malloc.h>
#include <sys/types.h>
#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"

extern int Emergency_stop;

struct el
{
	double gpt[9];
};

/* ============================================
   Generic adv-diffusion for temperature field.
   ============================================ */


void advection_diffusion_parameters(struct All_variables *E)
{
	/* Set intial values, defaults & read parameters */
  int m = E->parallel.me;

	E->advection.fixed_timestep = 0.0;
	E->advection.temp_iterations = 2;	/* petrov-galerkin iterations: minimum value. */
	E->advection.total_timesteps = 1;
	E->advection.sub_iterations = 1;
	E->advection.last_sub_iterations = 1;
	E->advection.gamma = 0.5;
	E->advection.ADVECTION = 1;
	
	input_boolean("ADV", &(E->advection.ADVECTION), "on",m);

	input_int("minstep", &(E->advection.min_timesteps), "1",m);
	input_int("maxstep", &(E->advection.max_timesteps), "1000",m);
	input_int("maxtotstep", &(E->advection.max_total_timesteps), "1000000",m);
	input_float("finetunedt", &(E->advection.fine_tune_dt), "0.9",m);
	input_float("fixed_timestep", &(E->advection.fixed_timestep), "0.0",m);
	input_int("adv_sub_iterations", &(E->advection.temp_iterations), "2,2,nomax",m);
	input_float("maxadvtime", &(E->advection.max_dimensionless_time), "10.0",m);

	input_float("sub_tolerance", &(E->advection.vel_substep_aggression), "0.005",m);
	input_int("maxsub", &(E->advection.max_substeps), "25",m);

	input_float("liddefvel", &(E->advection.lid_defining_velocity), "0.01",m);
	input_float("sublayerfrac", &(E->advection.sub_layer_sample_level), "0.5",m);
	E->control.adi_heating = 0;
	E->control.visc_heating = 0;
	input_int("adi_heating", &(E->control.adi_heating), "1",m);
	input_int("visc_heating", &(E->control.visc_heating), "1",m);

	E->viscosity.zcomp = E->control.Q0ER = 0;

	input_int("markers_per_ele", &(E->advection.markers_per_ele), "0",m);
	input_float("comp_depth", &(E->viscosity.zcomp), "0.0",m);
	input_float("Q0_enriched", &(E->control.Q0ER), "0.0",m);

	E->viscosity.zcomp = 1.0 - E->viscosity.zcomp;


	/* allocate memory */

	return;
}

void advection_diffusion_allocate_memory(struct All_variables *E)
{
	int i;

	E->Tdot = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));

	for(i = 1; i <= E->lmesh.nno; i++)
		E->Tdot[i] = 0.0;

	return;
}

void PG_timestep_particle(struct All_variables *E)
{
	float T_interior1;

	//int iredo, i, j, psc_pass, count, steps;
	int iredo, i, psc_pass, count;

	float *DTdot, *Tdot1, *T1, T_maxvaried;

	//static int loops_since_new_eta = 0;
	static int been_here = 0;
	static int on_off = 0;

	DTdot = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	Tdot1 = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	T1 = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));

	if(been_here++ == 0)
	{
	}


	if(on_off == 0)
	{
		E->advection.timesteps++;
		std_timestep(E);
		E->advection.total_timesteps++;
	}

	if(on_off == 1)
	{
		Runge_Kutta(E, E->C, E->V, on_off);
	}
	else if(on_off == 0)
	{
		for(i = 1; i <= E->lmesh.nno; i++)
		{
			T1[i] = E->T[i];
			Tdot1[i] = E->Tdot[i];
		}

		T_maxvaried = 1.01;

		T_interior1 = Tmax(E, E->T);

		E->advection.dt_reduced = 1.0;
		E->advection.last_sub_iterations = 1;

		count = 0;

		do
		{

			E->advection.timestep *= E->advection.dt_reduced;

			iredo = 0;

			if(E->advection.ADVECTION)
			{
				predictor(E, E->T, E->Tdot);

				for(psc_pass = 0; psc_pass < E->advection.temp_iterations; psc_pass++)
				{
					pg_solver(E, E->T, E->Tdot, DTdot, E->V, E->convection.heat_sources, 1.0, 1, E->TB, E->node);
					corrector(E, E->T, E->Tdot, DTdot);
				}
			}

			/* get the max temperature for new T */
			E->monitor.T_interior = Tmax(E, E->T);

			if(E->monitor.T_interior / T_interior1 > T_maxvaried)
			{
				for(i = 1; i <= E->lmesh.nno; i++)
				{
					E->T[i] = T1[i];
					E->Tdot[i] = Tdot1[i];
				}
				iredo = 1;
				E->advection.dt_reduced *= 0.5;
				E->advection.last_sub_iterations++;
			}

		} while(iredo == 1 && E->advection.last_sub_iterations <= 5);


		count++;

		temperatures_conform_bcs(E);

		E->advection.last_sub_iterations = count;


		Euler(E, E->C, E->V, on_off);

		E->monitor.elapsed_time += E->advection.timestep;
	}							/* end for on_off==0  */


	thermal_buoyancy(E);

	if(E->monitor.solution_cycles < E->advection.max_timesteps)
		E->control.keep_going = 1;
	else
		E->control.keep_going = 0;

	on_off = (on_off == 0) ? 1 : 0;

	free((void *)DTdot);		/* free memory for vel solver */
	free((void *)Tdot1);		/* free memory for vel solver */
	free((void *)T1);			/* free memory for vel solver */

	return;
}



void PG_timestep(struct All_variables *E)
{
	float T_interior1;

	//int iredo, i, j, psc_pass, count, steps;
	int iredo, i, psc_pass, count;

	float *DTdot, *Tdot1, *T1, T_maxvaried;

	//static int loops_since_new_eta = 0;
	static int been_here = 0;

	DTdot = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	Tdot1 = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	T1 = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));

	if(been_here++ == 0)
	{
	}

	E->advection.timesteps++;

	std_timestep(E);

	for(i = 1; i <= E->lmesh.nno; i++)
	{
		T1[i] = E->T[i];
		Tdot1[i] = E->Tdot[i];
	}

	T_maxvaried = 1.01;

	T_interior1 = Tmax(E, E->T);

	E->advection.dt_reduced = 1.0;
	E->advection.last_sub_iterations = 1;

	count = 0;

	do
	{

		E->advection.timestep *= E->advection.dt_reduced;

		iredo = 0;

		if(E->advection.ADVECTION)
		{
			predictor(E, E->T, E->Tdot);

			for(psc_pass = 0; psc_pass < E->advection.temp_iterations; psc_pass++)
			{
				pg_solver(E, E->T, E->Tdot, DTdot, E->V, E->convection.heat_sources, 1.0, 1, E->TB, E->node);
				corrector(E, E->T, E->Tdot, DTdot);
			}
		}

		/* get the max temperature for new T */
		E->monitor.T_interior = Tmax(E, E->T);

		if(E->monitor.T_interior / T_interior1 > T_maxvaried)
		{
			for(i = 1; i <= E->lmesh.nno; i++)
			{
				E->T[i] = T1[i];
				E->Tdot[i] = Tdot1[i];
			}
			iredo = 1;
			E->advection.dt_reduced *= 0.5;
			E->advection.last_sub_iterations++;
		}

	} while(iredo == 1 && E->advection.last_sub_iterations <= 5);

	E->advection.total_timesteps++;
	E->monitor.elapsed_time += E->advection.timestep;
	count++;


	temperatures_conform_bcs(E);
	thermal_buoyancy(E);

	E->advection.last_sub_iterations = count;

	if(E->monitor.solution_cycles < E->advection.max_timesteps)
		E->control.keep_going = 1;
	else
		E->control.keep_going = 0;

	free((void *)DTdot);		/* free memory for vel solver */
	free((void *)Tdot1);		/* free memory for vel solver */
	free((void *)T1);			/* free memory for vel solver */

	return;
}


/* ==============================
   predictor and corrector steps.
   ============================== */

void predictor(struct All_variables *E, float *field, float *fielddot)
{
	int node;
	float multiplier;

	multiplier = (1.0 - E->advection.gamma) * E->advection.timestep;

	for(node = 1; node <= E->lmesh.nno; node++)
	{
		if(!(E->node[node] & (TBX | TBZ | TBY)))
			field[node] += multiplier * fielddot[node];
		fielddot[node] = 0.0;
	}

	return;
}

void corrector(struct All_variables *E, float *field, float *fielddot, float *Dfielddot)

{
	int node;
	float multiplier;

	multiplier = E->advection.gamma * E->advection.timestep;

	for(node = 1; node <= E->lmesh.nno; node++)
	{
		if(!(E->node[node] & (TBX | TBZ | TBY)))
			field[node] += multiplier * Dfielddot[node];
		fielddot[node] += Dfielddot[node];


	}

	return;
}

/* ===================================================
   The solution step -- determine residual vector from
   advective-diffusive terms and solve for delta Tdot
   Two versions are available -- one for Cray-style 
   vector optimizations etc and one optimized for 
   workstations.
   =================================================== */


void pg_solver(struct All_variables *E, float *T, float *Tdot, float *DTdot, float **V, struct SOURCES Q0, float diff, int bc, float **TBC, unsigned int *FLAGS)
{
	int el, e, a, i, a1;
	double rtf[4][9], Eres[9];	/* correction to the (scalar) Tdot field */

	struct Shape_function PG;

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int ends = enodes[dims];

	for(i = 1; i <= E->lmesh.nno; i++)
		DTdot[i] = 0.0;

	for(el = 1; el <= E->lmesh.nel; el++)
	{

		if(E->control.Rsphere)
			get_rtf(E, el, 0, rtf, E->mesh.levmax);

		e = (el - 1) % E->lmesh.elz + 1;
		diff = (E->diffusivity[e] + E->diffusivity[e + 1]) * 0.5;

		pg_shape_fn(E, el, &PG, V, rtf, diff);
		element_residual(E, el, PG, V, T, Tdot, Q0, Eres, rtf, diff, TBC, FLAGS);

		for(a = 1; a <= ends; a++)
		{
			a1 = E->ien[el].node[a];
			DTdot[a1] += Eres[a];
		}

	}							/* next element */

	exchange_node_f20(E, DTdot, E->mesh.levmax);

	for(i = 1; i <= E->lmesh.nno; i++)
	{
		DTdot[i] *= E->Mass[i];	/* lumped mass matrix */
	}

	return;
}



/* ===================================================
   Petrov-Galerkin shape functions for a given element
   =================================================== */

void pg_shape_fn(struct All_variables *E, int el, struct Shape_function *PG, float **V, double rtf[4][9], float diffusion)
{
	int i, j;
	int *ienmatrix;

	//double xsi1, xsi2, xsi3;
	double uc1, uc2, uc3;
	//double aa1, aa2, aa3;
	//double ah1, ah2, ah3;
	double size1, size2, size3;
	//double neg1, neg2, neg3;
	double u1, u2, u3;
	//double adiff, dx1, dx2, dx3, uxse, ueta, ufai, xse, eta, fai;
	double adiff, uxse, ueta, ufai, xse, eta, fai;

	//double twodiff, prod, prod1;
	double twodiff, prod1;

	//double K, unorm;
	double unorm;

	//const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	//const int ends = enodes[dims];
	//const int vpts = vpoints[dims];

	ienmatrix = E->ien[el].node;

	twodiff = 2.0 * diffusion;


	size1 = (double)E->eco[el].size[1];
	size2 = (double)E->eco[el].size[2];
	size3 = (double)E->eco[el].size[3];

	uc1 = uc2 = uc3 = 0.0;

	for(i = 1; i <= ENODES3D; i++)
	{
		uc1 += E->N.ppt[GNPINDEX(i, 1)] * V[1][ienmatrix[i]];
		uc2 += E->N.ppt[GNPINDEX(i, 1)] * V[2][ienmatrix[i]];
		uc3 += E->N.ppt[GNPINDEX(i, 1)] * V[3][ienmatrix[i]];
	}

	uxse = fabs(uc1 * size1);
	ueta = fabs(uc2 * size2);
	ufai = fabs(uc3 * size3);

	xse = (uxse > twodiff) ? (1.0 - twodiff / uxse) : 0.0;
	eta = (ueta > twodiff) ? (1.0 - twodiff / ueta) : 0.0;
	fai = (ufai > twodiff) ? (1.0 - twodiff / ufai) : 0.0;

	unorm = uc1 * uc1 + uc2 * uc2 + uc3 * uc3;

	adiff = (unorm > 0.000001) ? ((uxse * xse + ueta * eta + ufai * fai) / (2.0 * unorm)) : 0.0;

	if(E->control.Rsphere)
	{
		for(i = 1; i <= VPOINTS3D; i++)
		{
			u1 = u2 = u3 = 0.0;

			rtf[2][i] = rtf[3][i] / sin(rtf[1][i]);	/* temporary */

			for(j = 1; j <= ENODES3D; j++)	/* this line heavily used */
			{
				u1 += V[1][ienmatrix[j]] * E->N.vpt[GNVINDEX(j, i)];
				u2 += V[2][ienmatrix[j]] * E->N.vpt[GNVINDEX(j, i)];
				u3 += V[3][ienmatrix[j]] * E->N.vpt[GNVINDEX(j, i)];
			}

			for(j = 1; j <= ENODES3D; j++)
			{
				prod1 = (u1 * E->gNX[el].vpt[GNVXINDEX(0, j, i)] * rtf[3][i] + u2 * E->gNX[el].vpt[GNVXINDEX(1, j, i)] * rtf[2][i] + u3 * E->gNX[el].vpt[GNVXINDEX(2, j, i)]);

				PG->vpt[GNVINDEX(j, i)] = E->N.vpt[GNVINDEX(j, i)] + adiff * prod1;
			}
		}
	}
	else if(E->control.CART3D)
	{
		for(i = 1; i <= VPOINTS3D; i++)
		{
			u1 = u2 = u3 = 0.0;

			for(j = 1; j <= ENODES3D; j++)	/* this line heavily used */
			{
				u1 += V[1][ienmatrix[j]] * E->N.vpt[GNVINDEX(j, i)];
				u2 += V[2][ienmatrix[j]] * E->N.vpt[GNVINDEX(j, i)];
				u3 += V[3][ienmatrix[j]] * E->N.vpt[GNVINDEX(j, i)];
			}

			for(j = 1; j <= ENODES3D; j++)
			{
				prod1 = (u1 * E->gNX[el].vpt[GNVXINDEX(0, j, i)] + u2 * E->gNX[el].vpt[GNVXINDEX(1, j, i)] + u3 * E->gNX[el].vpt[GNVXINDEX(2, j, i)]);

				PG->vpt[GNVINDEX(j, i)] = E->N.vpt[GNVINDEX(j, i)] + adiff * prod1;
			}
		}
	}


	return;
}



/* ==========================================
   Residual force vector from heat-transport.
   Used to correct the Tdot term.
   =========================================  */

void element_residual(struct All_variables *E, int el, struct Shape_function PG, float **vel, float *field, float *fielddot, struct SOURCES Q0, double Eres[9], double rtf[4][9], float diff, float **BC, unsigned int *FLAGS)
{
	//int i, j, a, k, node, nodes[5], d, aid, back_front, onedfns;
	int i, j, a, k, node, nodes[5], aid, back_front, onedfns;
	double Q;
	double dT[9];
	double tx1[9], tx2[9], tx3[9];
	double v1[9], v2[9], v3[9];
	//double adv_dT, t2[4];
	double T, DT;

	//register double prod, sfn;
	register double sfn;
	struct Shape_function1 GM;
	struct Shape_function1_dA dGamma;
	//double temp;

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int ends = enodes[dims];
	const int vpts = vpoints[dims];
	const int diffusion = (diff != 0.0);

	for(i = 1; i <= vpts; i++)
	{
		dT[i] = 0.0;
		v1[i] = tx1[i] = 0.0;
		v2[i] = tx2[i] = 0.0;
		v3[i] = tx3[i] = 0.0;
	}


	if(E->control.CART3D)
		for(j = 1; j <= ends; j++)
		{
			node = E->ien[el].node[j];
			T = field[node];
			if(E->node[node] & (TBX | TBY | TBZ))
				DT = 0.0;
			else
				DT = fielddot[node];

			for(i = 1; i <= vpts; i++)
			{
				dT[i] += DT * E->N.vpt[GNVINDEX(j, i)];
				tx1[i] += E->gNX[el].vpt[GNVXINDEX(0, j, i)] * T;
				tx2[i] += E->gNX[el].vpt[GNVXINDEX(1, j, i)] * T;
				tx3[i] += E->gNX[el].vpt[GNVXINDEX(2, j, i)] * T;
				sfn = E->N.vpt[GNVINDEX(j, i)];
				v1[i] += vel[1][node] * sfn;
				v2[i] += vel[2][node] * sfn;
				v3[i] += vel[3][node] * sfn;
			}
		}
	else if(E->control.Rsphere)
		for(j = 1; j <= ends; j++)
		{
			node = E->ien[el].node[j];
			T = field[node];
			if(E->node[node] & (TBX | TBY | TBZ))
				DT = 0.0;
			else
				DT = fielddot[node];

			for(i = 1; i <= vpts; i++)
			{
				dT[i] += DT * E->N.vpt[GNVINDEX(j, i)];
				tx1[i] += E->gNX[el].vpt[GNVXINDEX(0, j, i)] * T * rtf[3][i];
				tx2[i] += E->gNX[el].vpt[GNVXINDEX(1, j, i)] * T * rtf[2][i];
				tx3[i] += E->gNX[el].vpt[GNVXINDEX(2, j, i)] * T;
				sfn = E->N.vpt[GNVINDEX(j, i)];
				v1[i] += vel[1][node] * sfn;
				v2[i] += vel[2][node] * sfn;
				v3[i] += vel[3][node] * sfn;
			}
		}

	Q = E->rad_heat.total;
	if(E->control.composition)
		Q = (1 - E->CE[el]) * E->rad_heat.total + E->CE[el] * E->control.Q0ER;

	Q = (Q - E->heating_adi[el] + E->heating_visc[el]) * E->heating_latent[el];


/*
    Q=0.0;
    for(i=0;i<Q0.number;i++)
	Q += Q0.Q[i] * exp(-Q0.lambda[i] * (E->monitor.elapsed_time+Q0.t_offset));
*/

	/* construct residual from this information */


	if(diffusion)
	{
		if(E->control.CART3D)
			for(j = 1; j <= ends; j++)
			{
				Eres[j] = 0.0;
				for(i = 1; i <= vpts; i++)
					Eres[j] -= PG.vpt[GNVINDEX(j, i)] * E->gDA[el].vpt[i] * (dT[i] - Q + v1[i] * tx1[i] + v2[i] * tx2[i] + v3[i] * tx3[i]) + diff * E->heating_latent[el] * E->gDA[el].vpt[i] * (E->gNX[el].vpt[GNVXINDEX(0, j, i)] * tx1[i] + E->gNX[el].vpt[GNVXINDEX(1, j, i)] * tx2[i] + E->gNX[el].vpt[GNVXINDEX(2, j, i)] * tx3[i]);
			}
		else if(E->control.Rsphere)
			for(j = 1; j <= ends; j++)
			{
				Eres[j] = 0.0;
				for(i = 1; i <= vpts; i++)
					Eres[j] -= PG.vpt[GNVINDEX(j, i)] * E->gDA[el].vpt[i] * (dT[i] - Q + v1[i] * tx1[i] + v2[i] * tx2[i] + v3[i] * tx3[i]) + diff * E->heating_latent[el] * E->gDA[el].vpt[i] * (E->gNX[el].vpt[GNVXINDEX(0, j, i)] * tx1[i] * rtf[3][i] + E->gNX[el].vpt[GNVXINDEX(1, j, i)] * tx2[i] * rtf[2][i] + E->gNX[el].vpt[GNVXINDEX(2, j, i)] * tx3[i]);
			}
	}
	else
	{							/* no diffusion term */
		for(j = 1; j <= ends; j++)
		{
			Eres[j] = 0.0;
			for(i = 1; i <= vpts; i++)
				Eres[j] -= PG.vpt[GNVINDEX(j, i)] * E->gDA[el].vpt[i] * (dT[i] - Q + v1[i] * tx1[i] + v2[i] * tx2[i] + v3[i] * tx3[i]);
		}
	}

	/* See brooks etc: the diffusive term is excused upwinding for 
	 * rectangular elements  */

	/* include BC's for fluxes at (nominally horizontal) edges (X-Y plane) */

	if(FLAGS != NULL)
	{
		onedfns = 0;
		for(a = 1; a <= ends; a++)
			if(FLAGS[E->ien[el].node[a]] & FBZ)
			{
				if(!onedfns++)
					get_global_1d_shape_fn(E, el, &GM, &dGamma, 1);

				nodes[1] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][0];
				nodes[2] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][0];
				nodes[4] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][1];
				nodes[3] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][1];

				for(aid = 0, j = 1; j <= onedvpoints[E->mesh.nsd]; j++)
					if(a == nodes[j])
						aid = j;
				if(aid == 0)
					printf("%d: mixed up in pg-flux int: looking for %d\n", el, a);

				if(loc[a].plus[1] != 0)
					back_front = 0;
				else
					back_front = dims;

				for(j = 1; j <= onedvpoints[dims]; j++)
					for(k = 1; k <= onedvpoints[dims]; k++)
						Eres[a] += dGamma.vpt[GMVGAMMA(1 + back_front, j)] * E->M.vpt[GMVINDEX(aid, j)] * g_1d[j].weight[dims - 1] * BC[2][E->ien[el].node[a]] * E->M.vpt[GMVINDEX(k, j)];
			}
	}

	return;
}




/* =====================================================
   Obtain largest possible timestep (no melt considered)
   =====================================================  */


void std_timestep(struct All_variables *E)
{
	static int been_here = 0;
	//static float diff_timestep, root3, root2;
	static float diff_timestep;
	int i, d, n, nel;

	float adv_timestep;
	//float ts, uc2, uc3, uc1, uc, size, step;
	float ts, uc2, uc3, uc1, uc, step;

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	//const int ends = enodes[dims];

	nel = E->lmesh.nel;

	if(E->advection.fixed_timestep != 0.0)
	{
		E->advection.timestep = E->advection.fixed_timestep;
		return;
	}

	if(been_here == 0)
	{
		diff_timestep = 1.0e8;
		for(i = 1; i <= nel; i++)
		{
			for(d = 1; d <= dims; d++)
			{
				ts = E->eco[i].size[d] * E->eco[i].size[d];
				if(diff_timestep > ts)
					diff_timestep = ts;
			}
		}
		diff_timestep = 0.5 * global_fmin(E, diff_timestep);
	}

	adv_timestep = 1.0e8;
	for(i = 1; i <= nel; i++)
	{
		uc = uc1 = uc2 = uc3 = 0.0;

		if(3 == dims)
		{
			for(n = 1; n <= ENODES3D; n++)
			{
				uc1 += E->N.ppt[GNPINDEX(n, 1)] * E->V[1][E->ien[i].node[n]];
				uc2 += E->N.ppt[GNPINDEX(n, 1)] * E->V[2][E->ien[i].node[n]];
				uc3 += E->N.ppt[GNPINDEX(n, 1)] * E->V[3][E->ien[i].node[n]];
			}
			uc = fabs(uc1) / E->eco[i].size[1] + fabs(uc2) / E->eco[i].size[2] + fabs(uc3) / E->eco[i].size[3];
			step = 0.5 / uc;
			adv_timestep = min(adv_timestep, step);
		}
		else
		{
			for(n = 1; n <= ENODES2D; n++)
			{
				uc1 += E->N.ppt[GNPINDEX(n, 1)] * E->V[1][E->ien[i].node[n]];
				uc2 += E->N.ppt[GNPINDEX(n, 1)] * E->V[2][E->ien[i].node[n]];
			}
			uc = fabs(uc1) / E->eco[i].size[1] + fabs(uc2) / E->eco[i].size[2];
			step = 0.5 / uc;
			adv_timestep = min(adv_timestep, step);
		}
	}

	adv_timestep = 1.0e-32 + min(E->advection.fine_tune_dt * adv_timestep, diff_timestep);

	E->advection.timestep = global_fmin(E, adv_timestep);

	return;
}


void process_heating(struct All_variables *E)
{
	//int m, e, i, j, ee;
	int e, i, j, ee;
	static int been = 0;
	//static float para1;
	double temp1, temp2, temp3, temp4, temp5, temp6;
	FILE *fp;
	char filename[250];

	const int dims = E->mesh.nsd;
	const int ends = enodes[dims];
	//const int lev = E->mesh.levmax;
	//const int nno = E->lmesh.nno;
	const int vpts = vpoints[dims];


	if(been == 0)
	{
//		para1 = 1e6*E->data.radius_km*E->data.radius_km/(E->data.density*E->data.Cp*E->data.ref_temperature*E->data.therm_diff);
//		E->data.disptn_number = 1e3*E->data.therm_exp*E->data.grav_acc*E->data.radius_km/E->data.Cp;
		been++;
	}

	E->rad_heat.total = E->control.Q0;

	temp1 = E->data.disptn_number / E->control.Atemp;
	temp3 = temp4 = 0;

	if(E->control.visc_heating)
	{

		strain_rate_2_inv(E, E->heating_visc, 0);

		for(e = 1; e <= E->lmesh.nel; e++)
		{
			ee = (e - 1) % E->lmesh.elz + 1;
			temp2 = 0.0;
			for(i = 1; i <= vpts; i++)
				temp2 += E->EVi[(e - 1) * vpts + i];
			temp2 = temp2 / vpts;
			E->heating_visc[e] = temp1 * temp2 * E->heating_visc[e];

			temp4 = temp4 + E->heating_visc[e] * E->eco[e].area;

		}
	}

	if(E->control.adi_heating)
	{

		for(e = 1; e <= E->lmesh.nel; e++)
		{
			ee = (e - 1) % E->lmesh.elz + 1;

			temp2 = 0.0;
			for(i = 1; i <= ends; i++)
			{
				j = E->ien[e].node[i];
				temp2 = temp2 + E->V[3][j] * (E->T[j] + E->data.surf_temp) * E->data.disptn_number;
			}
			temp2 = temp2 / ends;
			E->heating_adi[e] = temp2 * (E->expansivity[ee] + E->expansivity[ee + 1]) * 0.5;

			temp3 = temp3 + E->heating_adi[e] * E->eco[e].area;

		}
	}

	MPI_Allreduce(&temp4, &temp6, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&temp3, &temp5, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	if(E->parallel.me == 0)
		fprintf(E->fp, "%g %g\n", temp5, temp6);

	if(E->control.Ra_670 != 0.0 || E->control.Ra_410 != 0)
	{
		for(e = 1; e <= E->lmesh.nel; e++)
			E->heating_latent[e] = 1.0;
	}

	if(E->control.Ra_670 != 0.0)
	{
		temp1 = 2.0 * E->control.width670 * E->control.clapeyron670 * E->control.Ra_670 / E->control.Atemp;

		for(e = 1; e <= E->lmesh.nel; e++)
		{
			temp2 = 0;
			temp3 = 0;
			for(i = 1; i <= ends; i++)
			{
				j = E->ien[e].node[i];
				temp2 = temp2 + temp1 * (1.0 - E->Fas670[j]) * E->Fas670[j] * E->V[3][j] * (E->T[j] + E->data.surf_temp) * E->data.disptn_number;
				temp3 = temp3 + temp1 * E->control.clapeyron670 * (1.0 - E->Fas670[j]) * E->Fas670[j] * (E->T[j] + E->data.surf_temp) * E->data.disptn_number;
			}
			temp2 = temp2 / ends;
			temp3 = temp3 / ends;
			E->heating_adi[e] += temp2;
			E->heating_latent[e] += temp3;
		}
	}


	if(E->control.Ra_410 != 0.0)
	{
		temp1 = 2.0 * E->control.width410 * E->control.clapeyron410 * E->control.Ra_410 / E->control.Atemp;

		for(e = 1; e <= E->lmesh.nel; e++)
		{
			temp2 = 0;
			temp3 = 0;
			for(i = 1; i <= ends; i++)
			{
				j = E->ien[e].node[i];
				temp2 = temp2 + temp1 * (1.0 - E->Fas410[j]) * E->Fas410[j] * E->V[3][j] * (E->T[j] + E->data.surf_temp) * E->data.disptn_number;
				temp3 = temp3 + temp1 * E->control.clapeyron410 * (1.0 - E->Fas410[j]) * E->Fas410[j] * (E->T[j] + E->data.surf_temp) * E->data.disptn_number;
			}
			temp2 = temp2 / ends;
			temp3 = temp3 / ends;
			E->heating_adi[e] += temp2;
			E->heating_latent[e] += temp3;
		}
	}

	if(E->control.Ra_670 != 0.0 || E->control.Ra_410 != 0)
	{
		for(e = 1; e <= E->lmesh.nel; e++)
			E->heating_latent[e] = 1.0 / E->heating_latent[e];
	}


	if(E->monitor.solution_cycles % 1000 == 0)
	{
#ifdef USE_GZDIR
	  if(E->control.gzdir)
		sprintf(filename, "%s/heating.%d", E->control.data_file2, E->parallel.me);
	  else
#endif
		sprintf(filename, "%s.heating.%d", E->control.data_file2, E->parallel.me);
	  
		fp = fopen(filename, "w");
		fprintf(fp, "QQ %g %g %g\n", E->control.Atemp, E->data.disptn_number, E->rad_heat.total);
		for(e = 1; e <= E->lmesh.nel; e++)
			fprintf(fp, "%d %g %g %g\n", e, E->heating_adi[e], E->heating_visc[e], E->heating_latent[e]);
		fclose(fp);
	}
/*
*/


	return;

}
