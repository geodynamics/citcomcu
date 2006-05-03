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
#include "Construct_arrays.h"

/* get_elt_k(), get_aug_k(), build_diagonal_of_K()
 * get_elt_g(), build_diagonal_of_Ahat() */
#include "Element_calculations.h"

/* exchange_id_d20() */
#include "Parallel_related.h" 

/* layers() */
#include "Viscosity_structures.h"

/* project_viscosity() */
#include "Solver_multigrid.h"


/*========================================================
  Function to make the IEN array for a mesh of given 
  dimension. IEN is an externally defined structure array

  NOTE: this is not really general enough for new elements:
  it should be done through a pre-calculated lookup table.
  ======================================================== */

void construct_ien(struct All_variables *E)
{
	int lev, p, q, r, rr, i, a, node, e;
	int element, start, nel, nno;
	int pmax, qmax, rmax;
	int pnmax, qnmax, rnmax;

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int ends = enodes[dims];

	for(lev = E->mesh.levmax; lev >= E->mesh.levmin; lev--)
	{
		pmax = E->lmesh.ELZ[lev];
		qmax = E->lmesh.ELX[lev];
		rmax = E->lmesh.ELY[lev];
		pnmax = E->lmesh.NOZ[lev];
		qnmax = E->lmesh.NOX[lev];
		rnmax = E->lmesh.NOY[lev];
		nel = E->lmesh.NEL[lev];
		nno = E->lmesh.NNO[lev];

		for(r = 1; r <= rmax; r++)
			for(q = 1; q <= qmax; q++)
				for(p = 1; p <= pmax; p++)
				{
					element = (r - 1) * pmax * qmax + (q - 1) * pmax + p;
					start = (r - 1) * pnmax * qnmax + (q - 1) * pnmax + p;
					for(rr = 1; rr <= ends; rr++)
						E->IEN[lev][element].node[rr] = start + offset[rr].vector[0] + offset[rr].vector[1] * pnmax + offset[rr].vector[2] * pnmax * qnmax;
				}


		/*   NEI is used for constructing Node_map for assembling del2 node by node; 
		 * and is used in Solver_multigrid for interp and project et al. */

		for(i = 1; i <= nno; i++)
			E->NEI[lev].nels[i] = 0;

		for(e = 1; e <= nel; e++)
			for(a = 1; a <= ends; a++)
			{
				node = E->IEN[lev][e].node[a];
				E->NEI[lev].nels[node]++;
				E->NEI[lev].element[(node - 1) * ends + E->NEI[lev].nels[node] - 1] = e;
				E->NEI[lev].lnode[(node - 1) * ends + E->NEI[lev].nels[node] - 1] = a;
			}

	}							/* end loop for lev */

	/*  determine surface things */
	e = 0;
	for(element = 1; element <= E->lmesh.NEL[E->mesh.levmax]; element++)
		if(element % E->lmesh.elz == 0)
		{
			e++;
			E->sien[e].node[1] = E->ien[element].node[5] / E->lmesh.noz;
			E->sien[e].node[2] = E->ien[element].node[6] / E->lmesh.noz;
			E->sien[e].node[3] = E->ien[element].node[7] / E->lmesh.noz;
			E->sien[e].node[4] = E->ien[element].node[8] / E->lmesh.noz;
			E->surf_element[e] = element;
		}
	E->lmesh.snel = e;
	for(i = 1; i <= E->lmesh.nsf; i++)
		E->surf_node[i] = i * E->lmesh.noz;


	if(E->control.verbose)
	{
		fprintf(E->fp, "output_IEN_arrays \n");
		if(dims == 2)
			for(i = 1; i <= E->mesh.nel; i++)
				fprintf(E->fp, "%d %d %d %d %d\n", i, E->ien[i].node[1], E->ien[i].node[2], E->ien[i].node[3], E->ien[i].node[4]);
		else if(dims == 3)
			for(i = 1; i <= E->mesh.nel; i++)
				fprintf(E->fp, "%d %d %d %d %d %d %d %d %d\n", i, E->ien[i].node[1], E->ien[i].node[2], E->ien[i].node[3], E->ien[i].node[4], E->ien[i].node[5], E->ien[i].node[6], E->ien[i].node[7], E->ien[i].node[8]);
	}

	return;
}

/*============================================
  Function to make the ID array for above case
  ============================================ */

void construct_id(struct All_variables *E)
{
	int i, j, k, i1, i2, j1, j2, k1, k2;
	int eqn_count, node;
	unsigned int doff;
	int lev;
	int nox, noy, noz;
	int elx, ely, elz;

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	//const int ends = enodes[dims];

	for(lev = E->mesh.levmax; lev >= E->mesh.levmin; lev--)
	{
		eqn_count = 0;
		elz = E->lmesh.ELZ[lev];
		elx = E->lmesh.ELX[lev];
		ely = E->lmesh.ELY[lev];
		nox = E->lmesh.NOX[lev];
		noz = E->lmesh.NOZ[lev];
		noy = E->lmesh.NOY[lev];

		for(node = 1; node <= E->lmesh.NNO[lev]; node++)
		{
			for(doff = 1; doff <= dims; doff++)
			{
				E->ID[lev][node].doff[doff] = eqn_count;
				eqn_count += 1;
			}
		}

		E->lmesh.NEQ[lev] = eqn_count;

		j1 = 1;
		j2 = nox - 1;
		if(E->parallel.me_loc[1] == E->parallel.nprocx - 1)
			j2 = nox;

		i1 = 1;
		i2 = noz - 1;
		if(E->parallel.me_loc[3] == E->parallel.nprocz - 1)
			i2 = noz;

		k1 = 1;
		k2 = noy - 1;
		if(E->parallel.me_loc[2] == E->parallel.nprocy - 1)
			k2 = noy;

		for(k = k1; k <= k2; k++)
			for(j = j1; j <= j2; j++)
				for(i = i1; i <= i2; i++)
				{
					node = i + (j - 1) * noz + (k - 1) * noz * nox;
					for(doff = 1; doff <= dims; doff++)
						E->parallel.IDD[lev][E->ID[lev][node].doff[doff]] = 1;
				}

	}

	E->lmesh.neq = E->lmesh.NEQ[E->mesh.levmax];	/*  Total NUMBER of independent variables  */


	lev = E->mesh.levmax;
	if(E->control.verbose)
	{
		fprintf(E->fp, "output_ID_arrays \n");
		if(dims == 2)
			for(i = 1; i <= E->lmesh.nno; i++)
				fprintf(E->fp, "%d %d %d %d \n", eqn_count, i, E->ID[lev][i].doff[1], E->ID[lev][i].doff[2]);
		else if(dims == 3)
			for(i = 1; i <= E->lmesh.nno; i++)
				fprintf(E->fp, "%d %d %d %d %d\n", eqn_count, i, E->ID[lev][i].doff[1], E->ID[lev][i].doff[2], E->ID[lev][i].doff[3]);
	}

	return;
}

/*==========================================================
  Function to construct  the LM array from the ID and IEN arrays 
  ========================================================== */

void construct_lm(struct All_variables *E)
{
	int lev;
	int a, e;
	//int lev, eqn_no;
	int nel, nel2;

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int ends = enodes[dims];


	for(lev = E->mesh.levmin; lev <= E->mesh.levmax; lev++)
	{
		nel = E->lmesh.NEL[lev];
		nel2 = nel / 2;
		for(e = 1; e <= nel; e++)
		{
			for(a = 1; a <= ends; a++)
			{
				E->LMD[lev][e].node[a].doff[1] = E->ID[lev][E->IEN[lev][e].node[a]].doff[1];
				E->LMD[lev][e].node[a].doff[2] = E->ID[lev][E->IEN[lev][e].node[a]].doff[2];
				E->LMD[lev][e].node[a].doff[3] = E->ID[lev][E->IEN[lev][e].node[a]].doff[3];
			}
		}
	}

	if(E->control.verbose)
	{
		if(dims == 3)
			for(lev = E->mesh.levmin; lev <= E->mesh.levmax; lev++)
			{
				nel = E->lmesh.NEL[lev];
				for(e = 1; e <= nel; e++)
				{
					for(a = 1; a <= ends; a++)
						fprintf(E->fp, "%d %d %d %d %d %d\n", lev, e, a, E->LMD[lev][e].node[a].doff[1], E->LMD[lev][e].node[a].doff[2], E->LMD[lev][e].node[a].doff[3]);
				}
			}
		else if(dims == 2)
			for(lev = E->mesh.levmin; lev <= E->mesh.levmax; lev++)
			{
				nel = E->lmesh.NEL[lev];
				for(e = 1; e <= nel; e++)
					for(a = 1; a <= ends; a++)
						fprintf(E->fp, "%d %d %d %d %d\n", lev, e, a, E->LMD[lev][e].node[a].doff[1], E->LMD[lev][e].node[a].doff[2]);
			}
	}
	return;
}


/* =====================================================
 *    Function to build the local node matrix indexing maps
 *
 *    Shijie Zhong modified it in 1998 to only store half of
 *    the matrix by  taking into account of the symmetry
 *    of stiffness matrix. This also affects how the matrix
 *    operation is done in assembling A*f, thus affecting
 *    a lot of routines.
 * ===================================================== */


void construct_node_maps(struct All_variables *E)
{
	//float initial_time;

	//int el, n, nn, lev, i, j, k, ja, jj, ii, kk, ia, is, ie, js, je, ks, ke, dims2;
	int nn, lev, i, j, k, ja, jj, ii, kk, ia, is, ie, js, je, ks, ke;
	//int doff, nox, noy, noz, noxz, node1, eqn1, loc1, count, found, element;
	int doff, nox, noy, noz, noxz;
	int neq, nno, matrix;
	//int *node_map;

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	//const int ends = enodes[dims];
	const int max_eqn = max_eqn_interaction[dims];

	for(lev = E->mesh.levmax; lev >= E->mesh.levmin; lev--)
	{
		neq = E->lmesh.NEQ[lev];
		nno = E->lmesh.NNO[lev];
		nox = E->lmesh.NOX[lev];
		noy = E->lmesh.NOY[lev];
		noz = E->lmesh.NOZ[lev];
		noxz = E->lmesh.NOZ[lev] * E->lmesh.NOX[lev];

		matrix = (nno + 3) * max_eqn;
		E->Node_map[lev] = (int *)malloc((matrix + 3) * sizeof(int));

		for(i = 0; i <= matrix; i++)
			E->Node_map[lev][i] = neq + 1;	/* DANGER !!! */

		for(ii = 1; ii <= noy; ii++)
			for(jj = 1; jj <= nox; jj++)
				for(kk = 1; kk <= noz; kk++)
				{
					nn = kk + (jj - 1) * noz + (ii - 1) * noxz;
					for(doff = 1; doff <= dims; doff++)
						E->Node_map[lev][(nn - 1) * max_eqn + doff - 1] = E->ID[lev][nn].doff[doff];
					/* the first three d.o.f. in Node_map are for the current node */

					ia = 0;
					is = (ii == 1) ? 2 : 1;
					ie = 2;
					js = (jj == 1) ? 2 : 1;
					je = (jj == nox) ? 2 : 3;
					ks = (kk == 1) ? 2 : 1;
					ke = (kk == noz) ? 2 : 3;

					for(i = is; i <= ie; i++)
						for(j = js; j <= je; j++)
							for(k = ks; k <= ke; k++)
							{
								ja = nn - ((2 - i) * noxz + (2 - j) * noz + 2 - k);
								/* the other d.o.f. in Node_map are for nodes with smaller
								 * node number less the current node */
								if(ja < nn)
								{
									ia++;
									for(doff = 1; doff <= dims; doff++)
										E->Node_map[lev][(nn - 1) * max_eqn + ia * dims + doff - 1] = E->ID[lev][ja].doff[doff];
								}
							}
				}

		E->Eqn_k1[lev] = (higher_precision *) malloc((matrix + 5) * sizeof(higher_precision));
		E->Eqn_k2[lev] = (higher_precision *) malloc((matrix + 5) * sizeof(higher_precision));
		if(dims == 3)
			E->Eqn_k3[lev] = (higher_precision *) malloc((matrix + 5) * sizeof(higher_precision));

		E->mesh.matrix_size[lev] = matrix + 1;
	}							/* end for level and m */

	return;
}


void construct_node_ks(struct All_variables *E)
{
	//int lev, level, i, j, k, e;
	int level, i, j, k;
	//int node, node1, eqn1, eqn2, eqn3, loc0, loc1, loc2, loc3, found, element, index, pp, qq;
	int node, node1, eqn1, eqn2, eqn3, loc0, found, element, index, pp, qq;
	int neq, nno, nel;

	double elt_K[24 * 24];
	//static int been_here = 0;
	double w1, w2, w3, ww1, ww2, ww3;

	//higher_precision *B1, *B2, *B3;
	higher_precision *B1, *B2;
	//int *C;

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int ends = enodes[dims];
	const int lms = loc_mat_size[E->mesh.nsd];
	const int max_eqn = max_eqn_interaction[dims];
	//const int vpts = vpoints[E->mesh.nsd];

	const double zero = 0.0;

	for(level = E->mesh.levmax; level >= E->mesh.levmin; level--)
	{
		neq = E->lmesh.NEQ[level];
		nel = E->lmesh.NEL[level];
		nno = E->lmesh.NNO[level];

		for(i = 0; i <= (neq + 1); i++)
			E->BI[level][i] = zero;
		for(i = 0; i <= E->mesh.matrix_size[level]; i++)
		{
			E->Eqn_k1[level][i] = zero;
			E->Eqn_k2[level][i] = zero;
			E->Eqn_k3[level][i] = zero;
		}

		for(element = 1; element <= nel; element++)
		{

			get_elt_k(E, element, elt_K, level, 0);

			if(E->control.augmented_Lagr)
				get_aug_k(E, element, elt_K, level);

			build_diagonal_of_K(E, element, elt_K, level);

			for(i = 1; i <= ends; i++)
			{					/* i, is the node we are storing to */
				node = E->IEN[level][element].node[i];

				pp = (i - 1) * dims;
				w1 = w2 = w3 = 1.0;

				loc0 = (node - 1) * max_eqn;

				if(E->NODE[level][node] & VBX)
					w1 = 0.0;
				if(E->NODE[level][node] & VBY)
					w2 = 0.0;
				if(E->NODE[level][node] & VBZ)
					w3 = 0.0;

				for(j = 1; j <= ends; j++)
				{				/* j is the node we are receiving from */
					node1 = E->IEN[level][element].node[j];
					if(node1 <= node)
					{
						ww1 = ww2 = ww3 = 1.0;
						qq = (j - 1) * dims;
						eqn1 = E->LMD[level][element].node[j].doff[1];
						eqn2 = E->LMD[level][element].node[j].doff[2];
						eqn3 = E->LMD[level][element].node[j].doff[3];

						if(E->NODE[level][node1] & VBX)
							ww1 = 0.0;
						if(E->NODE[level][node1] & VBY)
							ww2 = 0.0;
						if(E->NODE[level][node1] & VBZ)
							ww3 = 0.0;

						/* search for direction 1 */

						found = 0;
						for(k = 0; k < max_eqn; k++)
							if(E->Node_map[level][loc0 + k] == eqn1)
							{	/* found, index next equation */
								index = k;
								found++;
								break;
							}

						assert(found /* direction 1 */ );

						E->Eqn_k1[level][loc0 + index] += w1 * ww1 * elt_K[pp * lms + qq];	/* direction 1 */
						E->Eqn_k2[level][loc0 + index] += w2 * ww1 * elt_K[(pp + 1) * lms + qq];	/* direction 1 */
						E->Eqn_k3[level][loc0 + index] += w3 * ww1 * elt_K[(pp + 2) * lms + qq];	/* direction 1 */

						/* search for direction 2 */

						found = 0;
						for(k = 0; k < max_eqn; k++)
							if(E->Node_map[level][loc0 + k] == eqn2)
							{	/* found, index next equation */
								index = k;
								found++;
								break;
							}

						assert(found /* direction 2 */ );

						E->Eqn_k1[level][loc0 + index] += w1 * ww2 * elt_K[pp * lms + qq + 1];	/* direction 1 */
						E->Eqn_k2[level][loc0 + index] += w2 * ww2 * elt_K[(pp + 1) * lms + qq + 1];	/* direction 2 */
						E->Eqn_k3[level][loc0 + index] += w3 * ww2 * elt_K[(pp + 2) * lms + qq + 1];	/* direction 3 */

						/* search for direction 3 */

						found = 0;
						for(k = 0; k < max_eqn; k++)
							if(E->Node_map[level][loc0 + k] == eqn3)
							{	/* found, index next equation */
								index = k;
								found++;
								break;
							}


						assert(found /* direction 3 */ );

						E->Eqn_k1[level][loc0 + index] += w1 * ww3 * elt_K[pp * lms + qq + 2];	/* direction 1 */
						E->Eqn_k2[level][loc0 + index] += w2 * ww3 * elt_K[(pp + 1) * lms + qq + 2];	/* direction 2 */
						E->Eqn_k3[level][loc0 + index] += w3 * ww3 * elt_K[(pp + 2) * lms + qq + 2];	/* direction 3 */

					}			/* for node1< node */
				}
			}
		}						/* element */


		exchange_id_d20(E, E->BI[level], level);


		for(j = 0; j < neq; j++)
		{
			if(E->BI[level][j] == 0.0)
				fprintf(stderr, "me= %d level %d, equation %d/%d has zero diagonal term\n", E->parallel.me, level, j, neq);
			assert(E->BI[level][j] != 0 /* diagonal of matrix = 0, not acceptable */ );
			E->BI[level][j] = (float)1.0 / E->BI[level][j];
		}

		E->control.B_is_good[level] = 0;

		if(E->control.verbose)
		{
			fprintf(stderr, "output stiffness matrix!!!\n");
			fprintf(E->fp, "level %d\n", level);
			for(j = 1; j <= nno; j++)
			{
				eqn1 = E->ID[level][j].doff[1];
				eqn2 = E->ID[level][j].doff[2];
				loc0 = (j - 1) * max_eqn;
				B1 = E->Eqn_k1[level] + loc0;
				B2 = E->Eqn_k2[level] + loc0;
				for(i = 0; i < max_eqn; i++)
				{
					fprintf(E->fp, "%d %d %g %g\n", j, i, B1[i], B2[i]);
				}
			}
		}

	}

	return;
}



/* ============================================
   Function to set up the boundary condition
   masks and other indicators.
   ============================================  */

void construct_masks(struct All_variables *E)	/* Add lid/edge masks/nodal weightings */
{
	//int i, j, k, l, node, el, elt;
	int i, j, node;
	int lev, elx, elz, ely, nno, nox, noz, noy;

	for(lev = E->mesh.levmax; lev >= E->mesh.levmin; lev--)
	{
		elx = E->lmesh.ELX[lev];
		elz = E->lmesh.ELZ[lev];
		ely = E->lmesh.ELY[lev];
		nox = E->lmesh.NOX[lev];
		noy = E->lmesh.NOY[lev];
		noz = E->lmesh.NOZ[lev];
		nno = E->lmesh.NNO[lev];

		if(E->parallel.me_loc[3] == 0)
			for(i = 1; i <= nox; i++)	/* Horizontal  */
				for(j = 1; j <= noy; j++)
				{
					node = 1 + (i - 1) * noz + (j - 1) * noz * nox;
					E->NODE[lev][node] = E->NODE[lev][node] | TZEDGE;
				}
		if(E->parallel.me_loc[3] == E->parallel.nprocz - 1)
			for(i = 1; i <= nox; i++)	/* Horizontal  */
				for(j = 1; j <= noy; j++)
				{
					node = noz + (i - 1) * noz + (j - 1) * noz * nox;
					E->NODE[lev][node] = E->NODE[lev][node] | TZEDGE;
				}

		if(E->parallel.me_loc[2] == 0)
			for(i = 1; i <= noz; i++)	/* vertical */
				for(j = 1; j <= nox; j++)
				{
					node = i + (j - 1) * noz;
					E->NODE[lev][node] = E->NODE[lev][node] | TYEDGE;
				}
		if(E->parallel.me_loc[2] == E->parallel.nprocy - 1)
			for(i = 1; i <= noz; i++)
				for(j = 1; j <= nox; j++)
				{
					node = i + (j - 1) * noz + (noy - 1) * nox * noz;
					E->NODE[lev][node] = E->NODE[lev][node] | TYEDGE;
				}

		if(E->parallel.me_loc[1] == 0)	/* vertical */
			for(i = 1; i <= noz; i++)
				for(j = 1; j <= noy; j++)
				{
					node = i + (j - 1) * noz * nox;
					E->NODE[lev][node] = E->NODE[lev][node] | TXEDGE;
				}
		if(E->parallel.me_loc[1] == E->parallel.nprocx - 1)
			for(i = 1; i <= noz; i++)
				for(j = 1; j <= noy; j++)
				{
					node = i + (nox - 1) * noz + (j - 1) * nox * noz;
					E->NODE[lev][node] = E->NODE[lev][node] | TXEDGE;
				}


		for(i = 1; i <= nno; i++)
			E->TW[lev][i] = 0.0;

		for(i = 1; i <= nno; i++)
		{
			E->TW[lev][i] = enodes[E->mesh.nsd];
			if(E->NODE[lev][i] & TXEDGE)
				E->TW[lev][i] /= 2.0;
			if(E->NODE[lev][i] & TZEDGE)
				E->TW[lev][i] /= 2.0;
			if(E->NODE[lev][i] & TYEDGE)
				E->TW[lev][i] = max(1, E->TW[lev][i] / 2);
		}

		for(i = 1; i <= nno; i++)
		{
			assert(E->TW[lev][i] != 0.0 /* setting weightings failed */ );
			E->TW[lev][i] = 1.0 / (E->TW[lev][i]);
		}

	}

	return;
}


/*   ==========================================
     build the sub-element reference matrices
     ==========================================   */

void construct_sub_element(struct All_variables *E)
{
	int i, j, k, l;
	int lev, elx, elz, ely, elzu, elxu, elt, eltu;


	for(lev = E->mesh.levmax - 1; lev >= E->mesh.levmin; lev--)
	{
		elx = E->lmesh.ELX[lev];
		elz = E->lmesh.ELZ[lev];
		ely = E->lmesh.ELY[lev];
		elzu = 2 * elz;
		elxu = 2 * elx;

		for(i = 1; i <= elx; i++)
			for(j = 1; j <= elz; j++)
				for(k = 1; k <= ely; k++)
				{
					elt = j + (i - 1) * elz + (k - 1) * elz * elx;
					eltu = (j * 2 - 1) + elzu * 2 * (i - 1) + elxu * elzu * 2 * (k - 1);

					for(l = 1; l <= enodes[E->mesh.nsd]; l++)
					{
						E->EL[lev][elt].sub[l] = eltu + offset[l].vector[0] + offset[l].vector[1] * elzu + offset[l].vector[2] * elzu * elxu;
					}
				}
	}

	return;
}


void construct_elt_ks(struct All_variables *E)
{
	//int e, el, lev, j, k, ii;
	int el, lev, j, k, ii;

	const int dims = E->mesh.nsd;
	const int n = loc_mat_size[E->mesh.nsd];

	if(E->control.verbose && E->parallel.me == 0)
		fprintf(stderr, "storing elt k matrices\n");
	if(E->parallel.me == 0)
		fprintf(stderr, "storing elt k matrices\n");

	for(lev = E->mesh.levmin; lev <= E->mesh.levmax; lev++)
	{

		E->parallel.idb = 1;

		for(el = 1; el <= E->lmesh.NEL[lev]; el++)
		{

			get_elt_k(E, el, E->elt_k[lev][el].k, lev, 0);	/* not for penalty */

			if(E->control.augmented_Lagr)
				get_aug_k(E, el, E->elt_k[lev][el].k, lev);

			build_diagonal_of_K(E, el, E->elt_k[lev][el].k, lev);


		}

		exchange_id_d20(E, E->BI[lev], lev);	/*correct BI   */

		for(j = 0; j < E->lmesh.NEQ[lev]; j++)
		{
			if(E->BI[lev][j] == 0.0)
				fprintf(stderr, "me= %d level %d, equation %d/%d has zero diagonal term\n", E->parallel.me, lev, j, E->mesh.NEQ[lev]);
			assert(E->BI[lev][j] != 0 /* diagonal of matrix = 0, not acceptable */ );
			E->BI[lev][j] = (float)1.0 / E->BI[lev][j];
		}
	}

	if(E->control.verbose)
		for(lev = E->mesh.levmin; lev <= E->mesh.levmax; lev++)
			for(el = 1; el <= E->lmesh.NEL[lev]; el++)
				for(j = 1; j <= enodes[E->mesh.nsd]; j++)
					for(k = 1; k <= enodes[E->mesh.nsd]; k++)
					{
						ii = (j * n + k) * dims - (dims * n + dims);
						/*  fprintf(E->fp,"stiff_for_e %d %d %d %g %g %g %g \n",el,j,k,E->elt_k[lev][el].k[ii],E->elt_k[lev][el].k[ii+1],E->elt_k[lev][el].k[ii+n],E->elt_k[lev][el].k[ii+n+1]);      */
					}

	return;
}



void construct_elt_gs(struct All_variables *E)
{
	//int el, lev, a;
	int el, lev;

	//const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	//const int ends = enodes[dims];

	if(E->control.verbose && E->parallel.me == 0)
		fprintf(stderr, "storing elt g matrices\n");

	for(lev = E->mesh.levmin; lev <= E->mesh.levmax; lev++)
		for(el = 1; el <= E->lmesh.NEL[lev]; el++)
		{
			get_elt_g(E, el, E->elt_del[lev][el].g, lev);
		}

	return;
}



void construct_mat_group(struct All_variables *E)
{
	float rz_botm, rz_top, x3, t_b, slope1, slope2;
	float *Xtmp[4];

	//int i, j, k, el, lev, a, nodea, llayer, nslab, nwz, crit2;
	int i, el, a, nodea, llayer;

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int ends = enodes[dims];

	t_b = 3.0 * E->viscosity.zlith;

	if(E->control.Rsphere)
	{
		rz_botm = E->sphere.ri;
		rz_top = E->sphere.ro;
		for(i = 1; i <= E->mesh.nsd; i++)
			Xtmp[i] = E->SX[i];
	}
	else if(E->control.CART3D)
	{
		rz_botm = 0.0;
		rz_top = 1.0;
		for(i = 1; i <= E->mesh.nsd; i++)
			Xtmp[i] = E->X[i];
	}

	for(el = 1; el <= E->lmesh.nel; el++)
	{
		E->mat[el] = 1;
		x3 = 0;
		for(a = 1; a <= ends; a++)
		{
			nodea = E->ien[el].node[a];
			x3 += Xtmp[3][nodea];
		}
		x3 = x3 / ends;

		llayer = layers(E, x3);
		if(llayer)
		{						/* for layers:1-lith,2-upper and 3-lower mantle */
			E->mat[el] = llayer;
		}


		/*    if (E->mat[el]==1 || E->mat[el]==2)   {
		 * for (a=1;a<=ends;a++) {
		 * nodea = E->ien[el].node[a];
		 * llayer = weak_zones(E,nodea,t_b);
		 * if (llayer)  {
		 * E->mat[el] = 4;
		 * break;
		 * }
		 * }
		 * }
		 */

	}


	E->data.therm_exp_factor = 1.0 / E->data.therm_exp_factor;

	slope1 = (1.0 - E->data.therm_exp_factor) / (rz_top - rz_botm);
	slope2 = (1.0 - E->data.therm_diff_factor) / (rz_top - rz_botm);

	for(i = 1; i <= E->lmesh.noz; i++)
	{
		E->expansivity[i] = (slope1 * (Xtmp[3][i] - rz_botm) + E->data.therm_exp_factor);
		E->diffusivity[i] = (slope2 * (Xtmp[3][i] - rz_botm) + E->data.therm_diff_factor);
	}


	for(i = 1; i <= E->lmesh.noz; i++)
		fprintf(E->fp, "%d  %g %g\n", i, E->expansivity[i], E->diffusivity[i]);


/*
  for (el=1; el<=E->lmesh.nel; el++)
    fprintf(E->fp,"me= %d  mat[%d]= %d \n",E->parallel.me,el,E->mat[el]);
*/

	fflush(E->fp);

	return;
}


/* routine for constructing stiffness and node_maps */

void construct_stiffness_B_matrix(struct All_variables *E)
{
	static int been_here = 0;
	static int been_here0 = 0;

	int i;

	if(been_here0 == 0)
	{
		for(i = E->mesh.levmin; i <= E->mesh.levmax; i++)
			if(!E->control.NMULTIGRID && !E->control.NASSEMBLE)
			{
				E->elt_k[i] = (struct EK *)malloc((E->lmesh.NEL[i] + 1) * sizeof(struct EK));
			}
	}

	if(been_here0 == 0 || (E->viscosity.update_allowed && E->monitor.solution_cycles % E->control.KERNEL == 0))
	{
		/* do the following for the 1st time or update_allowed is true */

		if(E->control.NMULTIGRID)
			project_viscosity(E);

		construct_elt_gs(E);

		if(E->control.NMULTIGRID || E->control.NASSEMBLE)
		{
			if(been_here == 0)
			{					/* node_maps only built once */
				construct_node_maps(E);
				been_here = 1;
			}
			construct_node_ks(E);
		}
		else
		{
			construct_elt_ks(E);
		}

		build_diagonal_of_Ahat(E);




		if(E->control.NMULTIGRID || E->control.NASSEMBLE)
			rebuild_BI_on_boundary(E);


		been_here0 = 1;

	}

	return;
}


void rebuild_BI_on_boundary(struct All_variables *E)
{
	//int m, level, i, j;
	int level, i, j;
	int eqn1, eqn2, eqn3;

	higher_precision *B1, *B2, *B3;
	int *C;

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;

	const int max_eqn = max_eqn_interaction[dims];

	for(level = E->mesh.levmax; level >= E->mesh.levmin; level--)
	{
		for(j = 0; j < E->lmesh.NEQ[level]; j++)
			E->temp[j] = 0.0;

		for(i = 1; i <= E->lmesh.NNO[level]; i++)
			if(E->NODE[level][i] & LIDN)
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
					E->temp[eqn1] += fabs(B1[j]);
					E->temp[eqn2] += fabs(B2[j]);
					E->temp[eqn3] += fabs(B3[j]);
				}

				for(j = 0; j < max_eqn; j++)
					E->temp[C[j]] += fabs(B1[j]) + fabs(B2[j]) + fabs(B3[j]);

			}

		exchange_id_d20(E, E->temp, level);


		for(i = 0; i < E->lmesh.NEQ[level]; i++)
			E->temp[i] = E->temp[i] - 1.0 / E->BI[level][i];
		for(i = 1; i <= E->lmesh.NNO[level]; i++)
			if(E->NODE[level][i] & OFFSIDE)
			{
				eqn1 = E->ID[level][i].doff[1];
				eqn2 = E->ID[level][i].doff[2];
				eqn3 = E->ID[level][i].doff[3];
				E->BI[level][eqn1] = (double)1.0 / E->temp[eqn1];
				E->BI[level][eqn2] = (double)1.0 / E->temp[eqn2];
				E->BI[level][eqn3] = (double)1.0 / E->temp[eqn3];
			}
	}							/* end for level */

	return;
}
