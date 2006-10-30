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

/* Functions to assemble the element k matrices and the element f vector.
   Note that for the regular grid case the calculation of k becomes repetitive 
   to the point of redundancy. */

#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"
#include <sys/time.h>
#include <sys/resource.h>


/* *INDENT-OFF* */
static double Dl[5][5] = {  {0.0, 0.0, 0.0, 0.0, 0.0},
                            {0.0, 1.0, 1.0, 0.0, 1.0},
                            {0.0, 1.0, 1.0, 0.0, 1.0},
                            {0.0, 0.0, 0.0, 0.0, 0.0},
                            {0.0, 1.0, 1.0, 0.0, 1.0}  };

static double Dm[5][5] = {  {0.0, 0.0, 0.0, 0.0, 0.0},
                            {0.0, 2.0, 0.0, 0.0, 0.0},
                            {0.0, 0.0, 2.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0, 1.0, 0.0},
                            {0.0, 0.0, 0.0, 0.0, 2.0}  };
/* *INDENT-ON* */


/* ================================================================
   Function to assemble the global  F vector.
                     +
   Function to get the global H vector (mixed method driving terms)
   ================================================================ */

void assemble_forces(struct All_variables *E, int penalty)
{
	//double elt_f[24], elt_h[1];
	double elt_f[24];
	//int el, p, i, a, a1, a2, a3, e, ii, jj, kk, elx, ely, elz, node, temp_dims;
	int p, a, a1, a2, a3, e;

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int ends = enodes[E->mesh.nsd];
	const int neq = E->lmesh.neq;
	//const int npno = E->lmesh.npno;
	const int nel = E->lmesh.nel;
	const int lev = E->mesh.levmax;

	for(a = 0; a < neq; a++)
		E->F[a] = 0.0;

	/*for(a=0;a<npno;a++)
	 * E->H[a]=0.0; */

	for(e = 1; e <= nel; e++)
	{

		get_elt_f(E, e, elt_f, penalty, 1);
		/*get_elt_h(E,e,elt_h,penalty);
		 * E->H[e] = elt_h[0]; *//* due to single pressure node per element */

		for(a = 1; a <= ends; a++)
		{
			a1 = E->lm[e].node[a].doff[1];
			p = (a - 1) * dims;
			E->F[a1] += elt_f[p];
			a2 = E->lm[e].node[a].doff[2];
			E->F[a2] += elt_f[p + 1];
			a3 = E->lm[e].node[a].doff[3];
			E->F[a3] += elt_f[p + 2];
		}

	}

/*
  for(a=0;a<neq;a++)
   fprintf(E->fp,"bb %d %g\n",a,  E->F[a]); 
*/

	exchange_id_d20(E, E->F, lev);

	strip_bcs_from_residual(E, E->F, lev);

	return;
}



/*==============================================================
  Function to supply the element k matrix for a given element e.
  ==============================================================  */

void get_elt_k(struct All_variables *E, int el, double elt_k[24 * 24], int lev, int iconv)
{
	double bdbmu[4][4];
	//double bdbl[4][4];

	double rtf[4][9], W[9], ra[9], si[9], ct[9];
	//struct Shape_function GN;
	//struct Shape_function_dA dOmega;
	//struct Shape_function_dx GNx;
	static struct CC Cc;
	static struct CCX Ccx;


	//int p1[9], pn, qn, ad, bd;
	int pn, qn, ad, bd;

	//int nodea, nodeb, a, b, i, j, k, p, q, nint;
	int a, b, i, j, k;
	//double RM2[9], RMP[9], r;
	//double Visc, visc[9], temp;
	double temp;
	float gnx0, gnx1, gnx2;
	double shp, cc1, cc2, cc3;


	double ba[7][4][9][9];

	const int n = loc_mat_size[E->mesh.nsd];
	const int vpts = vpoints[E->mesh.nsd];
	//const int ppts = ppoints[E->mesh.nsd];
	const int ends = enodes[E->mesh.nsd];
	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	//const int sphere_key = 1;
	//const double zero = 0.0;
	const double one = 1.0;
	const double two = 2.0;

	for(k = 1; k <= vpts; k++)
	{
		W[k] = g_point[k].weight[dims - 1] * E->GDA[lev][el].vpt[k] * E->EVI[lev][(el - 1) * vpts + k];
	}

	if(E->control.Rsphere)
	{
		get_rtf(E, el, 0, rtf, lev);

		if(iconv == 1 || ((iconv == 0) && (el - 1) % E->lmesh.ELZ[lev] == 0))
			construct_c3x3matrix_el(E, el, &Cc, &Ccx, lev, 0);

		for(k = 1; k <= vpts; k++)
		{
			ra[k] = rtf[3][k];
			si[k] = one / sin(rtf[1][k]);
			ct[k] = cos(rtf[1][k]) * si[k];
		}

		for(a = 1; a <= ends; a++)
			for(i = 1; i <= dims; i++)
				for(k = 1; k <= VPOINTS3D; k++)
				{
					gnx0 = E->GNX[lev][el].vpt[GNVXINDEX(0, a, k)];
					gnx1 = E->GNX[lev][el].vpt[GNVXINDEX(1, a, k)];
					gnx2 = E->GNX[lev][el].vpt[GNVXINDEX(2, a, k)];
					shp = E->N.vpt[GNVINDEX(a, k)];
					cc1 = Cc.vpt[BVINDEX(1, i, a, k)];
					cc2 = Cc.vpt[BVINDEX(2, i, a, k)];
					cc3 = Cc.vpt[BVINDEX(3, i, a, k)];

					ba[1][i][a][k] = (gnx0 * cc1 + shp * Ccx.vpt[BVXINDEX(1, i, 1, a, k)] + shp * cc3) * ra[k];

					ba[2][i][a][k] = (shp * cc1 * ct[k] + shp * cc3 + (gnx1 * cc2 + shp * Ccx.vpt[BVXINDEX(2, i, 2, a, k)]) * si[k]) * ra[k];

					ba[3][i][a][k] = gnx2 * cc3;

					ba[4][i][a][k] = (gnx0 * cc2 + shp * Ccx.vpt[BVXINDEX(2, i, 1, a, k)] - shp * cc2 * ct[k] + (gnx1 * cc1 + shp * Ccx.vpt[BVXINDEX(1, i, 2, a, k)]) * si[k]) * ra[k];

					ba[5][i][a][k] = gnx2 * cc1 + (gnx0 * cc3 + shp * (Ccx.vpt[BVXINDEX(3, i, 1, a, k)] - cc1)) * ra[k];

					ba[6][i][a][k] = gnx2 * cc2 - ra[k] * shp * cc2 + (gnx1 * cc3 + shp * Ccx.vpt[BVXINDEX(3, i, 2, a, k)]) * si[k] * ra[k];
				}


	}							/* end for sphere */


	for(a = 1; a <= ends; a++)
		for(b = a; b <= ends; b++)
		{
			bdbmu[1][1] = bdbmu[1][2] = bdbmu[1][3] = bdbmu[2][1] = bdbmu[2][2] = bdbmu[2][3] = bdbmu[3][1] = bdbmu[3][2] = bdbmu[3][3] = 0.0;

			if(E->control.CART3D)
			{
				for(k = 1; k <= VPOINTS3D; k++)
				{
					bdbmu[1][1] += W[k] * E->GNX[lev][el].vpt[GNVXINDEX(0, a, k)] * E->GNX[lev][el].vpt[GNVXINDEX(0, b, k)];
					bdbmu[1][2] += W[k] * E->GNX[lev][el].vpt[GNVXINDEX(1, a, k)] * E->GNX[lev][el].vpt[GNVXINDEX(0, b, k)];
					bdbmu[1][3] += W[k] * E->GNX[lev][el].vpt[GNVXINDEX(2, a, k)] * E->GNX[lev][el].vpt[GNVXINDEX(0, b, k)];
					bdbmu[2][1] += W[k] * E->GNX[lev][el].vpt[GNVXINDEX(0, a, k)] * E->GNX[lev][el].vpt[GNVXINDEX(1, b, k)];
					bdbmu[2][2] += W[k] * E->GNX[lev][el].vpt[GNVXINDEX(1, a, k)] * E->GNX[lev][el].vpt[GNVXINDEX(1, b, k)];
					bdbmu[2][3] += W[k] * E->GNX[lev][el].vpt[GNVXINDEX(2, a, k)] * E->GNX[lev][el].vpt[GNVXINDEX(1, b, k)];
					bdbmu[3][1] += W[k] * E->GNX[lev][el].vpt[GNVXINDEX(0, a, k)] * E->GNX[lev][el].vpt[GNVXINDEX(2, b, k)];
					bdbmu[3][2] += W[k] * E->GNX[lev][el].vpt[GNVXINDEX(1, a, k)] * E->GNX[lev][el].vpt[GNVXINDEX(2, b, k)];
					bdbmu[3][3] += W[k] * E->GNX[lev][el].vpt[GNVXINDEX(2, a, k)] * E->GNX[lev][el].vpt[GNVXINDEX(2, b, k)];
				}

				temp = 0.0;
				for(k = 1; k <= VPOINTS3D; k++)
					temp += W[k] * (E->GNX[lev][el].vpt[GNVXINDEX(0, a, k)] * E->GNX[lev][el].vpt[GNVXINDEX(0, b, k)] + E->GNX[lev][el].vpt[GNVXINDEX(1, a, k)] * E->GNX[lev][el].vpt[GNVXINDEX(1, b, k)] + E->GNX[lev][el].vpt[GNVXINDEX(2, a, k)] * E->GNX[lev][el].vpt[GNVXINDEX(2, b, k)]);

				bdbmu[1][1] += temp;
				bdbmu[2][2] += temp;
				bdbmu[3][3] += temp;

			}					/* end for Cart3D */

			else if(E->control.Rsphere)
			{
				for(i = 1; i <= dims; i++)
					for(j = 1; j <= dims; j++)
						for(k = 1; k <= VPOINTS3D; k++)
						{
							bdbmu[i][j] += W[k] * (two * (ba[1][i][a][k] * ba[1][j][b][k] + ba[2][i][a][k] * ba[2][j][b][k] + ba[3][i][a][k] * ba[3][j][b][k]) + ba[4][i][a][k] * ba[4][j][b][k] + ba[5][i][a][k] * ba[5][j][b][k] + ba[6][i][a][k] * ba[6][j][b][k]);
						}
			}					/* end for Sphere */

			ad = dims * (a - 1);
			bd = dims * (b - 1);
			pn = ad * n + bd;
			qn = bd * n + ad;

			elt_k[pn] = bdbmu[1][1];	/* above */
			elt_k[pn + 1] = bdbmu[1][2];
			elt_k[pn + 2] = bdbmu[1][3];
			elt_k[pn + n] = bdbmu[2][1];
			elt_k[pn + n + 1] = bdbmu[2][2];
			elt_k[pn + n + 2] = bdbmu[2][3];
			elt_k[pn + 2 * n] = bdbmu[3][1];
			elt_k[pn + 2 * n + 1] = bdbmu[3][2];
			elt_k[pn + 2 * n + 2] = bdbmu[3][3];

			elt_k[qn] = bdbmu[1][1];	/* below diag */
			elt_k[qn + 1] = bdbmu[2][1];
			elt_k[qn + 2] = bdbmu[3][1];
			elt_k[qn + n] = bdbmu[1][2];
			elt_k[qn + n + 1] = bdbmu[2][2];
			elt_k[qn + n + 2] = bdbmu[3][2];
			elt_k[qn + 2 * n] = bdbmu[1][3];
			elt_k[qn + 2 * n + 1] = bdbmu[2][3];
			elt_k[qn + 2 * n + 2] = bdbmu[3][3];
		 /**/}					/*  Sum over all the a,b's to obtain full  elt_k matrix */

	return;
}


	/* =============================================
	 * General calling function for del_squared: 
	 * according to whether it should be element by
	 * element or node by node.
	 * ============================================= */

void assemble_del2_u(struct All_variables *E, double *u, double *Au, int level, int strip_bcs)
{
	if(E->control.NMULTIGRID || E->control.NASSEMBLE)
		n_assemble_del2_u(E, u, Au, level, strip_bcs);
	else
		e_assemble_del2_u(E, u, Au, level, strip_bcs);

	return;
}

	/* ======================================
	 * Assemble del_squared_u vector el by el
	 * ======================================   */

void e_assemble_del2_u(struct All_variables *E, double *u, double *Au, int level, int strip_bcs)
{
	//int el, e, i, a, b, a1, a2, a3, ii;
	int e, i, a, b, a1, a2, a3, ii;
	//double elt_k[24 * 24], U[24], AU[24], alpha = 1.0, beta = 0.0;

	//int indx[24];
	//char uplo = 'U';

	//double U1[24];

	const int n = loc_mat_size[E->mesh.nsd];
	const int ends = enodes[E->mesh.nsd];
	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int nel = E->lmesh.NEL[level];
	const int neq = E->lmesh.NEQ[level];

	for(i = 0; i < neq; i++)
		Au[i] = 0.0;

	for(e = 1; e <= nel; e++)
	{

		for(a = 1; a <= ends; a++)
		{
			a1 = E->LMD[level][e].node[a].doff[1];
			a2 = E->LMD[level][e].node[a].doff[2];
			a3 = E->LMD[level][e].node[a].doff[3];
			for(b = 1; b <= ends; b++)
			{
				ii = (a * n + b) * dims - (dims * n + dims);
				/* i=1, j=1,2,3 */
				Au[a1] += E->elt_k[level][e].k[ii] * u[E->LMD[level][e].node[b].doff[1]] + E->elt_k[level][e].k[ii + 1] * u[E->LMD[level][e].node[b].doff[2]] + E->elt_k[level][e].k[ii + 2] * u[E->LMD[level][e].node[b].doff[3]];
				/* i=2, j=1,2,3 */
				Au[a2] += E->elt_k[level][e].k[ii + n] * u[E->LMD[level][e].node[b].doff[1]] + E->elt_k[level][e].k[ii + n + 1] * u[E->LMD[level][e].node[b].doff[2]] + E->elt_k[level][e].k[ii + n + 2] * u[E->LMD[level][e].node[b].doff[3]];
				/* i=3, j=1,2,3 */
				Au[a3] += E->elt_k[level][e].k[ii + n + n] * u[E->LMD[level][e].node[b].doff[1]] + E->elt_k[level][e].k[ii + n + n + 1] * u[E->LMD[level][e].node[b].doff[2]] + E->elt_k[level][e].k[ii + n + n + 2] * u[E->LMD[level][e].node[b].doff[3]];

			}					/* end for loop b */
		}						/* end for loop a */

	}

	exchange_id_d20(E, Au, level);


	if(strip_bcs)
		strip_bcs_from_residual(E, Au, level);

	return;
}


	/* ======================================================
	 * Assemble Au using stored, nodal coefficients.
	 * ====================================================== */

void n_assemble_del2_u(struct All_variables *E, double *u, double *Au, int level, int strip_bcs)
{
	//int node, e, i, eqn1, eqn2, eqn3, loc0, loc1, loc2, loc3;
	int e, i, eqn1, eqn2, eqn3, loc0; 

	double U1, U2, U3, UU;

	static int been_here = 0;

	int *C;
	higher_precision *B1, *B2, *B3;

	const int neq = E->lmesh.NEQ[level];
	const int nno = E->lmesh.NNO[level];
	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int max_eqn = max_eqn_interaction[dims];

	for(e = 0; e <= neq; e++)
		Au[e] = 0.0;

	loc0 = 1;

	for(e = 1; e <= nno; e++)
	{

		eqn1 = E->ID[level][e].doff[1];
		eqn2 = E->ID[level][e].doff[2];
		eqn3 = E->ID[level][e].doff[3];

		U1 = u[eqn1];
		U2 = u[eqn2];
		U3 = u[eqn3];

		C = E->Node_map[level] + (e - 1) * max_eqn;
		B1 = E->Eqn_k1[level] + (e - 1) * max_eqn;
		B2 = E->Eqn_k2[level] + (e - 1) * max_eqn;
		B3 = E->Eqn_k3[level] + (e - 1) * max_eqn;

		for(i = 3; i < max_eqn; i++)
		{
			UU = u[C[i]];
			Au[eqn1] += B1[i] * UU;
			Au[eqn2] += B2[i] * UU;
			Au[eqn3] += B3[i] * UU;
		}
		/* contributions to the current node or equation from other
		 * adjacent nodes or d.o.f.. Use horizontal entries of the
		 * stiffness matrix. Not completed yet, since we only store
		 * half of the matrix */

		for(i = 0; i < max_eqn; i++)
			Au[C[i]] += B1[i] * U1 + B2[i] * U2 + B3[i] * U3;

		/* contributions from the current node to other nodes or
		 * eqns including current node (i=0 to max_eqn). 
		 * Use vertical entries of the stiffness matrix and the
		 * symmetry. This ultimately will complete the assembly */
	}

	exchange_id_d20(E, Au, level);

	if(strip_bcs)
		strip_bcs_from_residual(E, Au, level);

	been_here++;

	return;
}


void build_diagonal_of_K(struct All_variables *E, int el, double elt_k[24 * 24], int level)

{
	int a, a1, a2, p;

	const int n = loc_mat_size[E->mesh.nsd];
	const int dims = E->mesh.nsd;
	const int ends = enodes[E->mesh.nsd];

	for(a = 1; a <= ends; a++)
	{
		/* dirn 1 */
		a1 = E->LMD[level][el].node[a].doff[1];
		p = (a - 1) * dims;
		E->BI[level][a1] += elt_k[p * n + p];

		/* dirn 2 */
		a2 = E->LMD[level][el].node[a].doff[2];
		p = (a - 1) * dims + 1;
		E->BI[level][a2] += elt_k[p * n + p];

		/* dirn 3 */
		a1 = E->LMD[level][el].node[a].doff[3];
		p = (a - 1) * dims + 2;
		E->BI[level][a1] += elt_k[p * n + p];
	}

	return;
}

void build_diagonal_of_Ahat(struct All_variables *E)
{
	double BU;
	int e, npno, neq;
	int level;
	//float time, time0;

	for(level = E->mesh.levmax; level >= E->mesh.levmin; level--)
	{

		npno = E->lmesh.NPNO[level];
		neq = E->lmesh.NEQ[level];

		for(e = 1; e <= npno; e++)
			E->BPI[level][e] = 1.0;

		if(!E->control.precondition)
			return;

		for(e = 1; e <= npno; e++)
		{
			BU = assemble_dAhatp_entry(E, e, level);
			if(BU != 0.0)
				E->BPI[level][e] = 1.0 / BU;
			else
				E->BPI[level][e] = 1.0;
		}

	}

	return;
}

	/* ==========================================
	 * Assemble a div_u vector element by element
	 * ==========================================  */

void assemble_div_u(struct All_variables *E, double *U, double *divU, int level)
{
	//int e, j1, j2, j3, p, a, b;
	int e, j1, j2, j3, p, a;
	//higher_precision elt_g[24][1];

	const int nel = E->lmesh.NEL[level];
	const int ends = enodes[E->mesh.nsd];
	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int npno = E->lmesh.NPNO[level];

	for(e = 1; e <= npno; e++)
		divU[e] = 0.0;

	for(a = 1; a <= ends; a++)
	{
		p = (a - 1) * dims;
		for(e = 1; e <= nel; e++)
		{
			j1 = E->LMD[level][e].node[a].doff[1];
			j2 = E->LMD[level][e].node[a].doff[2];
			j3 = E->LMD[level][e].node[a].doff[3];
			/* for(b=0;b<ploc_mat_size[E->mesh.nsd];b++) */
			divU[e] += E->elt_del[level][e].g[p][0] * U[j1] + E->elt_del[level][e].g[p + 1][0] * U[j2] + E->elt_del[level][e].g[p + 2][0] * U[j3];
		}
	}

	return;
}


	/* ==========================================
	 * Assemble a grad_P vector element by element
	 * ==========================================  */

void assemble_grad_p(struct All_variables *E, double *P, double *gradP, int lev)
{
	//int el, e, i, j1, j2, p, a, nel, neq;
	int e, i, j1, j2, p, a, nel, neq;
	//higher_precision elt_g[24][1];

	const int ends = enodes[E->mesh.nsd];
	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;

	nel = E->lmesh.NEL[lev];
	neq = E->lmesh.NEQ[lev];

	for(i = 0; i <= neq; i++)
		gradP[i] = 0.0;

	for(e = 1; e <= nel; e++)
	{

		if(0.0 == P[e])
			continue;

		for(a = 1; a <= ends; a++)
		{
			p = (a - 1) * dims;
			j1 = E->LMD[lev][e].node[a].doff[1];
			j2 = E->LMD[lev][e].node[a].doff[2];
			/*for(b=0;b<ploc_mat_size[E->mesh.nsd];b++)  */
			gradP[j1] += E->elt_del[lev][e].g[p][0] * P[e];
			gradP[j2] += E->elt_del[lev][e].g[p + 1][0] * P[e];

			j1 = E->LMD[lev][e].node[a].doff[3];
			gradP[j1] += E->elt_del[lev][e].g[p + 2][0] * P[e];
		}

	}

	exchange_id_d20(E, gradP, lev);	/*  correct gradP   */

	strip_bcs_from_residual(E, gradP, lev);

	return;
}

double assemble_dAhatp_entry(struct All_variables *E, int e, int level)
{
	//int i, j, p, a, b, node, ee, element, lnode, npno;
	int i, j, p, a, b, npno;
	//higher_precision elt_g[24][1];

	double gradP[81], divU;

	const int ends = enodes[E->mesh.nsd];
	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;

	npno = E->lmesh.NPNO[level];

	for(i = 0; i < 81; i++)
		gradP[i] = 0.0;

	divU = 0.0;


	for(a = 1; a <= ends; a++)
	{
		p = (a - 1) * dims;
		j = E->LMD[level][e].node[a].doff[1];
		gradP[p] += E->BI[level][j] * E->elt_del[level][e].g[p][0];

		j = E->LMD[level][e].node[a].doff[2];
		gradP[p + 1] += E->BI[level][j] * E->elt_del[level][e].g[p + 1][0];

		j = E->LMD[level][e].node[a].doff[3];
		gradP[p + 2] += E->BI[level][j] * E->elt_del[level][e].g[p + 2][0];
	}


	/* calculate div U from the same thing .... */

	/* only need to run over nodes with non-zero grad P, i.e. the ones in
	 * the element accessed above, BUT it is only necessary to update the
	 * value in the original element, because the diagonal is all we use at
	 * the end ... */

	for(b = 1; b <= ends; b++)
	{
		p = (b - 1) * dims;
		divU += E->elt_del[level][e].g[p][0] * gradP[p];
		divU += E->elt_del[level][e].g[p + 1][0] * gradP[p + 1];
		divU += E->elt_del[level][e].g[p + 2][0] * gradP[p + 2];
	}




	return (divU);
}


/*==============================================================
  Function to supply the element g matrix for a given element e.
  ==============================================================  */

void get_elt_g(struct All_variables *E, int el, higher_precision elt_del[24][1], int lev)
{
	//double dGNdash[3];
	//double recip_radius, temp;
	double temp;
	//int p, a, nint, es, d, i, j, k;
	int p, a, i;
	double ra, ct, si, x[4], rtf[4][9];
	//int lmsize;

	//struct Shape_function GN;
	//struct Shape_function_dA dOmega;
	//struct Shape_function_dx GNx;
	static struct CC Cc;
	static struct CCX Ccx;

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int ends = enodes[dims];
	//const int vpts = vpoints[dims];
	//const int sphere_key = 1;

	/* Special case, 4/8 node bilinear cartesian square/cube element -> 1 pressure point */

	temp = p_point[1].weight[dims - 1] * E->GDA[lev][el].ppt[1];

	/* unroll etc some more */

	if(E->control.Rsphere)
	{
		if((el - 1) % E->lmesh.ELZ[lev] == 0)
			construct_c3x3matrix_el(E, el, &Cc, &Ccx, lev, 1);

		get_rtf(E, el, 2, rtf, lev);

		ra = rtf[3][1];
		si = 1.0 / sin(rtf[1][1]);
		ct = cos(rtf[1][1]) * si;

		for(a = 1; a <= ends; a++)
		{
			for(i = 1; i <= dims; i++)
				x[i] = E->GNX[lev][el].ppt[GNPXINDEX(2, a, 1)] * Cc.ppt[BPINDEX(3, i, a, 1)]
					+ 2.0 * ra * E->N.ppt[GNPINDEX(a, 1)] * Cc.ppt[BPINDEX(3, i, a, 1)] + ra * (E->GNX[lev][el].ppt[GNPXINDEX(0, a, 1)] * Cc.ppt[BPINDEX(1, i, a, 1)] + E->N.ppt[GNPINDEX(a, 1)] * Ccx.ppt[BPXINDEX(1, i, 1, a, 1)] + ct * E->N.ppt[GNPINDEX(a, 1)] * Cc.ppt[BPINDEX(1, i, a, 1)] + si * (E->GNX[lev][el].ppt[GNPXINDEX(1, a, 1)] * Cc.ppt[BPINDEX(2, i, a, 1)] + E->N.ppt[GNPINDEX(a, 1)] * Ccx.ppt[BPXINDEX(2, i, 2, a, 1)]));

			p = dims * (a - 1);
			elt_del[p][0] = -x[1] * temp;
			elt_del[p + 1][0] = -x[2] * temp;
			elt_del[p + 2][0] = -x[3] * temp;

		}

	}
	else if(E->control.CART3D)
	{

		for(a = 1; a <= ends; a++)
		{
			p = dims * (a - 1);
			elt_del[p][0] = -E->GNX[lev][el].ppt[GNPXINDEX(0, a, 1)] * temp;
			elt_del[p + 1][0] = -E->GNX[lev][el].ppt[GNPXINDEX(1, a, 1)] * temp;
			elt_del[p + 2][0] = -E->GNX[lev][el].ppt[GNPXINDEX(2, a, 1)] * temp;
		}
	}

	return;
}




	/* ===============================================================
	 * Function to create the element pressure-forcing vector (due
	 * to imposed velocity boundary conditions, mixed method).
	 * =============================================================== */

void get_elt_h(struct All_variables *E, int el, double elt_h[1], int penalty)
{
	//int aid, i, p, a, b, d, j, k, q, global, got_g;
	int i, p, a, b, q, got_g;
	unsigned int type;
	higher_precision elt_g[24][1];

	for(p = 0; p < 1; p++)
		elt_h[p] = 0.0;

	if(penalty)
		return;					/* no h term at all */

	got_g = 0;

	type = VBX;
	for(i = 1; i <= E->mesh.nsd; i++)
	{
		for(a = 1; a <= enodes[E->mesh.nsd]; a++)
		{
			if(E->node[E->ien[el].node[a]] & type)
			{
				if(!got_g)
				{
					get_elt_g(E, el, elt_g, E->mesh.levmax);
					got_g++;
				}

				p = E->mesh.nsd * (a - 1) + i - 1;
				for(b = 1; b <= pnodes[E->mesh.nsd]; b++)
				{
					q = b - 1;
					elt_h[q] -= elt_g[p][q] * E->VB[i][E->ien[el].node[a]];
				}
			}
		}
		type *= (unsigned int)2;	/* !!! depends on how x, y, and z are ordered */
	}
	return;
}

	/*=================================================================
	  Function to create the element force vector (allowing for b.c.'s)
	  ================================================================= */

void get_elt_f(struct All_variables *E, int el, double elt_f[24], int penalty, int bcs)
{

	//int aid, i, p, a, b, d, j, k, q, es;
	int i, p, a, b, j, k, q;
	//int node[5], back_front, got_elt_k, nodea, nodeb;
	int got_elt_k, nodea, nodeb;
	unsigned int type[4];
	//static int been_here = 0;

	//double force[9], force_at_gs[9], stress[9];
	double force[9], force_at_gs[9];
	//double vector[4], magnitude;
	//double tmp, rtf[4][9];
	//double rtf[4][9];
	double elt_k[24 * 24];

	//struct Shape_function GN;
	//struct Shape_function_dA dOmega;
	//struct Shape_function_dx GNx;
	//struct Shape_function1 GM;
	//struct Shape_function1_dA dGammax;
	static struct CC Cc;
	static struct CCX Ccx;

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int n = loc_mat_size[dims];
	const int ends = enodes[dims];
	const int vpts = vpoints[dims];
	//const int sphere_key = 1;

	type[1] = VBX;
	type[2] = VBY;
	type[3] = VBZ;


	for(p = 0; p < n; p++)
		elt_f[p] = 0.0;

	for(p = 1; p <= ends; p++)
		force[p] = E->buoyancy[E->ien[el].node[p]];

	for(j = 1; j <= vpts; j++)
	{							/*compute force at each int point */
		force_at_gs[j] = 0.0;
		for(k = 1; k <= ends; k++)
			force_at_gs[j] += force[k] * E->N.vpt[GNVINDEX(k, j)];
	}

	if(E->control.CART3D)
	{
		for(i = 1; i <= dims; i++)
		{
			for(a = 1; a <= ends; a++)
			{
				nodea = E->ien[el].node[a];
				p = dims * (a - 1) + i - 1;

				if(i == 3)
					for(j = 1; j <= vpts; j++)	/*compute sum(Na(j)*F(j)*det(j)) */
						elt_f[p] += force_at_gs[j] * E->N.vpt[GNVINDEX(a, j)] * E->gDA[el].vpt[j] * g_point[j].weight[dims - 1];

				/* imposed velocity terms */

				if(bcs)
				{
					got_elt_k = 0;
					for(j = 1; j <= dims; j++)
					{
						for(b = 1; b <= ends; b++)
						{
							nodeb = E->ien[el].node[b];
							if((E->node[nodeb] & type[j]) && (E->VB[j][nodeb] != 0.0))
							{
								if(!got_elt_k)
								{
									get_elt_k(E, el, elt_k, E->mesh.levmax, 1);
									got_elt_k = 1;
								}
								q = dims * (b - 1) + j - 1;
								if(p != q)
								{
									elt_f[p] -= elt_k[p * n + q] * E->VB[j][nodeb];
/*	 fprintf(E->fp,"el %d: dirn=%d, vbc found at node %d of %g\n",el,j,b,E->VB[j][nodeb]); */
								}
							}
						}		/* end for b */
					}			/* end for j */
				}				/* end if for if bcs */

			}
		}						/*  Complete the loops for a,i    */

	}							/* end for CART3D */

	else if(E->control.Rsphere)
	{

		//get_rtf(E, el, 0, rtf, E->mesh.levmax);
		if((el - 1) % E->lmesh.elz == 0)
			construct_c3x3matrix_el(E, el, &Cc, &Ccx, E->mesh.levmax, 0);

		for(i = 1; i <= dims; i++)
		{
			for(a = 1; a <= ends; a++)
			{
				nodea = E->ien[el].node[a];
				p = dims * (a - 1) + i - 1;

				for(j = 1; j <= vpts; j++)	/*compute sum(Na(j)*F(j)*det(j)) */
					elt_f[p] += force_at_gs[j] * E->N.vpt[GNVINDEX(a, j)] * E->gDA[el].vpt[j] * g_point[j].weight[dims - 1] * Cc.vpt[BVINDEX(3, i, a, j)];;

				/* imposed velocity terms */

				if(bcs)
				{
					got_elt_k = 0;
					for(j = 1; j <= dims; j++)
					{
						for(b = 1; b <= ends; b++)
						{
							nodeb = E->ien[el].node[b];
							if((E->node[nodeb] & type[j]) && (E->VB[j][nodeb] != 0.0))
							{
								if(!got_elt_k)
								{
									get_elt_k(E, el, elt_k, E->mesh.levmax, 1);
									got_elt_k = 1;
								}
								q = dims * (b - 1) + j - 1;
								if(p != q)
								{
									elt_f[p] -= elt_k[p * n + q] * E->VB[j][nodeb];
/*	 fprintf(E->fp,"el %d: dirn=%d, vbc found at node %d of %g\n",el,j,b,E->VB[j][nodeb]); */
								}
							}
						}		/* end for b */
					}			/* end for j */
				}				/* end if for if bcs */

			}
		}						/*  Complete the loops for a,i    */

	}							/* end for sphere */
	return;
}

	/* =================================================================
	 * subroutine to get augmented lagrange part of stiffness matrix
	 * ================================================================== */

void get_aug_k(struct All_variables *E, int el, double elt_k[24 * 24], int level)
{
	//int i, j, k, p[9], a, b, nodea, nodeb;
	int i, p[9], a, b;
	double Visc;

	const int n = loc_mat_size[E->mesh.nsd];
	const int ends = enodes[E->mesh.nsd];
	const int vpts = vpoints[E->mesh.nsd];
	const int dims = E->mesh.nsd;

	Visc = 0.0;
	for(a = 1; a <= vpts; a++)
	{
		p[a] = (a - 1) * dims;
		Visc += E->EVI[level][(el - 1) * vpts + a];
	}
	Visc = Visc / vpts;

	for(a = 1; a <= ends; a++)
		for(b = 1; b <= ends; b++)
		{
			i = (a - 1) * n * dims + (b - 1) * dims;
			elt_k[i] += Visc * E->control.augmented * E->elt_del[level][el].g[p[a]][0] * E->elt_del[level][el].g[p[b]][0];	/*for 11 */
			elt_k[i + 1] += Visc * E->control.augmented * E->elt_del[level][el].g[p[a]][0] * E->elt_del[level][el].g[p[b] + 1][0];	/* for 12 */
			elt_k[i + n] += Visc * E->control.augmented * E->elt_del[level][el].g[p[a] + 1][0] * E->elt_del[level][el].g[p[b]][0];	/* for 21 */
			elt_k[i + n + 1] += Visc * E->control.augmented * E->elt_del[level][el].g[p[a] + 1][0] * E->elt_del[level][el].g[p[b] + 1][0];	/* for 22 */

			elt_k[i + 2] += Visc * E->control.augmented * E->elt_del[level][el].g[p[a]][0] * E->elt_del[level][el].g[p[b] + 2][0];	/* for 13 */
			elt_k[i + n + 2] += Visc * E->control.augmented * E->elt_del[level][el].g[p[a] + 1][0] * E->elt_del[level][el].g[p[b] + 2][0];	/* for 23 */
			elt_k[i + n + n] += Visc * E->control.augmented * E->elt_del[level][el].g[p[a] + 2][0] * E->elt_del[level][el].g[p[b]][0];	/* for 31 */
			elt_k[i + n + n + 1] += Visc * E->control.augmented * E->elt_del[level][el].g[p[a] + 2][0] * E->elt_del[level][el].g[p[b] + 1][0];	/* for 32 */
			elt_k[i + n + n + 2] += Visc * E->control.augmented * E->elt_del[level][el].g[p[a] + 2][0] * E->elt_del[level][el].g[p[b] + 2][0];	/* for 33 */
		}

	return;
}
