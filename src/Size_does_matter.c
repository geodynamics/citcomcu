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

/*   	This is where the scaling functions and grid related things are kept.

		Louis Moresi aka LUIGI   6.xii.1989                */

/* Some relevant functions to spherical geometry were added here in 2004 */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include <mpi.h>


void twiddle_thumbs(struct All_variables *yawn, int scratch_groin)
{								/* Do nothing, just sit back and relax.
								 * Take it easy for a while, maybe size
								 * doesn't matter after all. There, there
								 * that's better. Now ... */

	return;
}

void get_global_shape_fn(struct All_variables *E, int el, int pressure, double rtf[4][9], int sphere, int level)
{
	int i, j, k, d, e;
	//double scale1, scale2, scale3;
	//double area;
	double jacobian;

	double LGNX[4][9];

	double dxda[4][4], cof[4][4], x[4], bc[4][4];


	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int ends = enodes[dims];
	const int vpts = vpoints[dims];
	const int ppts = ppoints[dims];
	//const int spts = spoints[dims];


	if(pressure < 2)
	{
		for(k = 1; k <= vpts; k++)
		{						/* all of the vpoints */
			for(d = 1; d <= dims; d++)
				for(e = 1; e <= dims; e++)
					dxda[d][e] = 0.0;

			for(i = 1; i <= ends; i++)
				for(d = 1; d <= dims; d++)
					for(e = 1; e <= dims; e++)
						dxda[d][e] += E->XX[level][e][E->IEN[level][el].node[i]] * E->Nx.vpt[GNVXINDEX(d - 1, i, k)];	/* This is Shijie's change (d<->e) */

			jacobian = determinant(dxda, E->mesh.nsd);
			E->GDA[level][el].vpt[k] = jacobian;

			for(d = 1; d <= dims; d++)
				for(e = 1; e <= dims; e++)
					cof[d][e] = cofactor(dxda, d, e, dims);

			for(j = 1; j <= ends; j++)
				for(d = 1; d <= dims; d++)
				{
					E->GNX[level][el].vpt[GNVXINDEX(d - 1, j, k)] = 0.0;
					for(e = 1; e <= dims; e++)
						E->GNX[level][el].vpt[GNVXINDEX(d - 1, j, k)] += E->Nx.vpt[GNVXINDEX(e - 1, j, k)] * cof[e][d];	/* switch e and d for cof  -- Shijie's  */
					E->GNX[level][el].vpt[GNVXINDEX(d - 1, j, k)] /= jacobian;
				}

			if(sphere)
			{
				for(d = 1; d <= dims; d++)
				{
					x[d] = 0;
					for(i = 1; i <= ends; i++)
						x[d] += E->XX[level][d][E->IEN[level][el].node[i]] * E->N.vpt[GNVINDEX(i, k)];
				}

				form_rtf_bc(k, x, rtf, bc);
				for(j = 1; j <= ends; j++)
					for(d = 1; d <= dims; d++)
						LGNX[d][j] = bc[d][1] * E->GNX[level][el].vpt[GNVXINDEX(0, j, k)] + bc[d][2] * E->GNX[level][el].vpt[GNVXINDEX(1, j, k)] + bc[d][3] * E->GNX[level][el].vpt[GNVXINDEX(2, j, k)];
				for(j = 1; j <= ends; j++)
					for(d = 1; d <= dims; d++)
						E->GNX[level][el].vpt[GNVXINDEX(d - 1, j, k)] = LGNX[d][j];
			}

		}
	}

	if(pressure > 0 && pressure < 3)
	{
		for(k = 1; k <= ppts; k++)
		{						/* all of the ppoints */
			for(d = 1; d <= dims; d++)
				for(e = 1; e <= dims; e++)
					dxda[d][e] = 0.0;

			for(i = 1; i <= ends; i++)
				for(d = 1; d <= dims; d++)
					for(e = 1; e <= dims; e++)
						dxda[d][e] += E->XX[level][e][E->IEN[level][el].node[i]] * E->Nx.ppt[GNPXINDEX(d - 1, i, k)];

			jacobian = determinant(dxda, E->mesh.nsd);
			E->GDA[level][el].ppt[k] = jacobian;

			for(d = 1; d <= dims; d++)
				for(e = 1; e <= dims; e++)
					cof[d][e] = cofactor(dxda, d, e, E->mesh.nsd);

			for(j = 1; j <= ends; j++)
				for(d = 1; d <= dims; d++)
				{
					E->GNX[level][el].ppt[GNPXINDEX(d - 1, j, k)] = 0.0;
					for(e = 1; e <= dims; e++)
						E->GNX[level][el].ppt[GNPXINDEX(d - 1, j, k)] += E->Nx.ppt[GNPXINDEX(e - 1, j, k)] * cof[e][d];
					E->GNX[level][el].ppt[GNPXINDEX(d - 1, j, k)] /= jacobian;
				}
			if(sphere)
			{
				for(d = 1; d <= dims; d++)
				{
					x[d] = 0;
					for(i = 1; i <= ends; i++)
						x[d] += E->XX[level][d][E->IEN[level][el].node[i]] * E->N.ppt[GNPINDEX(i, k)];
				}
				form_rtf_bc(k, x, rtf, bc);
				for(j = 1; j <= ends; j++)
					for(d = 1; d <= dims; d++)
						LGNX[d][j] = bc[d][1] * E->GNX[level][el].ppt[GNPXINDEX(0, j, k)] + bc[d][2] * E->GNX[level][el].ppt[GNPXINDEX(1, j, k)] + bc[d][3] * E->GNX[level][el].ppt[GNPXINDEX(2, j, k)];
				for(j = 1; j <= ends; j++)
					for(d = 1; d <= dims; d++)
						E->GNX[level][el].ppt[GNPXINDEX(d - 1, j, k)] = LGNX[d][j];
			}
		}
	}


	return;
}


void form_rtf_bc(int k, double x[4], double rtf[4][9], double bc[4][4])
{
	rtf[3][k] = 1.0 / sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
	rtf[1][k] = acos(x[3] * rtf[3][k]);
	rtf[2][k] = myatan(x[2], x[1]);

	bc[1][1] = x[3] * cos(rtf[2][k]);
	bc[1][2] = x[3] * sin(rtf[2][k]);
	bc[1][3] = -sin(rtf[1][k]) / rtf[3][k];
	bc[2][1] = -x[2];
	bc[2][2] = x[1];
	bc[2][3] = 0.0;
	bc[3][1] = x[1] * rtf[3][k];
	bc[3][2] = x[2] * rtf[3][k];
	bc[3][3] = x[3] * rtf[3][k];

	return;
}



void get_rtf(struct All_variables *E, int el, int pressure, double rtf[4][9], int lev)
{
	//int i, j, k, d, e;
	int i, k, d;
	double x[4];

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int ends = enodes[dims];
	const int vpts = vpoints[dims];
	const int ppts = ppoints[dims];

	if(pressure == 0)
	{
		for(k = 1; k <= vpts; k++)
		{						/* all of the vpoints */
			for(d = 1; d <= dims; d++)
				x[d] = 0.0;

			for(d = 1; d <= dims; d++)
				for(i = 1; i <= ends; i++)
					x[d] += E->XX[lev][d][E->IEN[lev][el].node[i]] * E->N.vpt[GNVINDEX(i, k)];

			rtf[3][k] = 1.0 / sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
			rtf[1][k] = acos(x[3] * rtf[3][k]);
			rtf[2][k] = myatan(x[2], x[1]);
		}
	}
	else
	{
		for(k = 1; k <= ppts; k++)
		{						/* all of the ppoints */
			for(d = 1; d <= dims; d++)
				x[d] = 0.0;

			for(d = 1; d <= dims; d++)
				for(i = 1; i <= ends; i++)
					x[d] += E->XX[lev][d][E->IEN[lev][el].node[i]] * E->N.ppt[GNPINDEX(i, k)];

			rtf[3][k] = 1.0 / sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
			rtf[1][k] = acos(x[3] * rtf[3][k]);
			rtf[2][k] = myatan(x[2], x[1]);
		}
	}


	return;
}

void construct_c3x3matrix_el(struct All_variables *E, int el, struct CC *cc, struct CCX *ccx, int lev, int pressure)
{
	//int a, i, j, k, d, e, es, nel_surface;
	int a, i, j, k, d;
	double x[4], u[4][4], ux[3][4][4], ua[4][4];
	double costt, cosff, sintt, sinff, rr, tt, ff;

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int ends = enodes[dims];
	const int vpts = vpoints[dims];
	const int ppts = ppoints[dims];

	if(pressure == 0)
	{
		for(k = 1; k <= vpts; k++)
		{						/* all of the vpoints */
			for(d = 1; d <= dims; d++)
				x[d] = 0.0;

			for(d = 1; d <= dims; d++)
				for(a = 1; a <= ends; a++)
					x[d] += E->XX[lev][d][E->IEN[lev][el].node[a]] * E->N.vpt[GNVINDEX(a, k)];

			rr = sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
			tt = acos(x[3] / rr);
			ff = myatan(x[2], x[1]);

			costt = cos(tt);
			cosff = cos(ff);
			sintt = sin(tt);
			sinff = sin(ff);

			u[1][1] = costt * cosff;
			u[1][2] = costt * sinff;
			u[1][3] = -sintt;
			u[2][1] = -sinff;
			u[2][2] = cosff;
			u[2][3] = 0.0;
			u[3][1] = sintt * cosff;
			u[3][2] = sintt * sinff;
			u[3][3] = costt;

			ux[1][1][1] = -sintt * cosff;
			ux[1][1][2] = -sintt * sinff;
			ux[1][1][3] = -costt;
			ux[2][1][1] = -costt * sinff;
			ux[2][1][2] = costt * cosff;
			ux[2][1][3] = 0.0;
			ux[1][2][1] = 0.0;
			ux[1][2][2] = 0.0;
			ux[1][2][3] = 0.0;
			ux[2][2][1] = -cosff;
			ux[2][2][2] = -sinff;
			ux[2][2][3] = 0.0;
			ux[1][3][1] = costt * cosff;
			ux[1][3][2] = costt * sinff;
			ux[1][3][3] = -sintt;
			ux[2][3][1] = -sintt * sinff;
			ux[2][3][2] = sintt * cosff;
			ux[2][3][3] = 0.0;

			for(a = 1; a <= ends; a++)
			{
				tt = E->SXX[lev][1][E->IEN[lev][el].node[a]];
				ff = E->SXX[lev][2][E->IEN[lev][el].node[a]];
				costt = cos(tt);
				cosff = cos(ff);
				sintt = sin(tt);
				sinff = sin(ff);

				ua[1][1] = costt * cosff;
				ua[1][2] = costt * sinff;
				ua[1][3] = -sintt;
				ua[2][1] = -sinff;
				ua[2][2] = cosff;
				ua[2][3] = 0.0;
				ua[3][1] = sintt * cosff;
				ua[3][2] = sintt * sinff;
				ua[3][3] = costt;

				for(i = 1; i <= dims; i++)
					for(j = 1; j <= dims; j++)
					{
						cc->vpt[BVINDEX(i, j, a, k)] = ua[j][1] * u[i][1] + ua[j][2] * u[i][2] + ua[j][3] * u[i][3];
						ccx->vpt[BVXINDEX(i, j, 1, a, k)] = ua[j][1] * ux[1][i][1] + ua[j][2] * ux[1][i][2] + ua[j][3] * ux[1][i][3];
						ccx->vpt[BVXINDEX(i, j, 2, a, k)] = ua[j][1] * ux[2][i][1] + ua[j][2] * ux[2][i][2] + ua[j][3] * ux[2][i][3];
/*            cc->vpt[BVINDEX(i,j,a,k)] = 0.0;
              if (i==j)cc->vpt[BVINDEX(i,j,a,k)] = 1.0;
              ccx->vpt[BVXINDEX(i,j,1,a,k)] = 0.0;
              ccx->vpt[BVXINDEX(i,j,2,a,k)] = 0.0;
*/
					}
			}					/* end for local node */

		}						/* end for int points */
	}							/* end if */

	else if(pressure)
	{

		for(k = 1; k <= ppts; k++)
		{						/* all of the ppoints */
			for(d = 1; d <= dims; d++)
				x[d] = 0.0;

			for(d = 1; d <= dims; d++)
				for(a = 1; a <= ends; a++)
					x[d] += E->XX[lev][d][E->IEN[lev][el].node[a]] * E->N.ppt[GNPINDEX(a, k)];

			rr = sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
			tt = acos(x[3] / rr);
			ff = myatan(x[2], x[1]);

			costt = cos(tt);
			cosff = cos(ff);
			sintt = sin(tt);
			sinff = sin(ff);

			u[1][1] = costt * cosff;
			u[1][2] = costt * sinff;
			u[1][3] = -sintt;
			u[2][1] = -sinff;
			u[2][2] = cosff;
			u[2][3] = 0.0;
			u[3][1] = sintt * cosff;
			u[3][2] = sintt * sinff;
			u[3][3] = costt;

			ux[1][1][1] = -sintt * cosff;
			ux[1][1][2] = -sintt * sinff;
			ux[1][1][3] = -costt;
			ux[2][1][1] = -costt * sinff;
			ux[2][1][2] = costt * cosff;
			ux[2][1][3] = 0.0;
			ux[1][2][1] = 0.0;
			ux[1][2][2] = 0.0;
			ux[1][2][3] = 0.0;
			ux[2][2][1] = -cosff;
			ux[2][2][2] = -sinff;
			ux[2][2][3] = 0.0;
			ux[1][3][1] = costt * cosff;
			ux[1][3][2] = costt * sinff;
			ux[1][3][3] = -sintt;
			ux[2][3][1] = -sintt * sinff;
			ux[2][3][2] = sintt * cosff;
			ux[2][3][3] = 0.0;

			for(a = 1; a <= ends; a++)
			{
				tt = E->SXX[lev][1][E->IEN[lev][el].node[a]];
				ff = E->SXX[lev][2][E->IEN[lev][el].node[a]];
				costt = cos(tt);
				cosff = cos(ff);
				sintt = sin(tt);
				sinff = sin(ff);

				ua[1][1] = costt * cosff;
				ua[1][2] = costt * sinff;
				ua[1][3] = -sintt;
				ua[2][1] = -sinff;
				ua[2][2] = cosff;
				ua[2][3] = 0.0;
				ua[3][1] = sintt * cosff;
				ua[3][2] = sintt * sinff;
				ua[3][3] = costt;

				for(i = 1; i <= dims; i++)
					for(j = 1; j <= dims; j++)
					{
						cc->ppt[BPINDEX(i, j, a, k)] = ua[j][1] * u[i][1] + ua[j][2] * u[i][2] + ua[j][3] * u[i][3];
						ccx->ppt[BPXINDEX(i, j, 1, a, k)] = ua[j][1] * ux[1][i][1] + ua[j][2] * ux[1][i][2] + ua[j][3] * ux[1][i][3];
						ccx->ppt[BPXINDEX(i, j, 2, a, k)] = ua[j][1] * ux[2][i][1] + ua[j][2] * ux[2][i][2] + ua[j][3] * ux[2][i][3];
/*              cc->ppt[BPINDEX(i,j,a,k)] = 0.0;
              if (i==j)cc->ppt[BPINDEX(i,j,a,k)] = 1.0;
              ccx->ppt[BPXINDEX(i,j,1,a,k)] =0.0;
              ccx->ppt[BPXINDEX(i,j,2,a,k)] =0.0;
 */ }

			}					/* end for local node */

		}						/* end for int points */


	}							/* end if pressure  */

	return;
}




/*   ======================================================================
     ======================================================================  */
void get_global_1d_shape_fn(struct All_variables *E, int el, struct Shape_function1 *GM, struct Shape_function1_dA *dGammax, int top)
{
	int i, k, d, e, ii;
	//int dirn, locn, node;
	int node;
	//int collapsed_dirn[2];
	//double scale[4];

	double jacobian;
	const int oned = onedvpoints[E->mesh.nsd];


	//double to, fo, dxdy[4][4], avet, aver, dxda[4][4], cof[4][4], xx[4][5];
	double to, fo, dxdy[4][4], dxda[4][4], xx[4][5];

	for(ii = 0; ii <= top; ii++)
	{							/* ii=0 for bottom ii=1 for top */

		if(E->control.CART3D)
		{
			for(i = 1; i <= onedvpoints[E->mesh.nsd]; i++)
			{					/* nodes */
				e = i + ii * onedvpoints[E->mesh.nsd];
				xx[1][i] = E->X[1][E->ien[el].node[e]];
				xx[2][i] = E->X[2][E->ien[el].node[e]];
			}
			for(k = 1; k <= onedvpoints[E->mesh.nsd]; k++)
			{					/* all of the vpoints */

				for(d = 1; d <= E->mesh.nsd - 1; d++)
					for(e = 1; e <= E->mesh.nsd - 1; e++)
						dxda[d][e] = 0.0;

				for(i = 1; i <= onedvpoints[E->mesh.nsd]; i++)	/* nodes */
					for(d = 1; d <= E->mesh.nsd - 1; d++)
						for(e = 1; e <= E->mesh.nsd - 1; e++)
							dxda[d][e] += xx[e][i] * E->Mx.vpt[GMVXINDEX(d - 1, i, k)];

				jacobian = determinant(dxda, E->mesh.nsd - 1);
				dGammax->vpt[GMVGAMMA(ii, k)] = jacobian;
			}
		}
		else if(E->control.Rsphere)
		{
			to = E->eco[el].centre[1];
			fo = E->eco[el].centre[2];

			dxdy[1][1] = cos(to) * cos(fo);
			dxdy[1][2] = cos(to) * sin(fo);
			dxdy[1][3] = -sin(to);
			dxdy[2][1] = -sin(fo);
			dxdy[2][2] = cos(fo);
			dxdy[2][3] = 0.0;
			dxdy[3][1] = sin(to) * cos(fo);
			dxdy[3][2] = sin(to) * sin(fo);
			dxdy[3][3] = cos(to);
			for(i = 1; i <= oned; i++)
			{					/* nodes */
				e = i + ii * oned;
				node = E->ien[el].node[e];
				xx[1][i] = E->X[1][node] * dxdy[1][1] + E->X[2][node] * dxdy[1][2] + E->X[3][node] * dxdy[1][3];
				xx[2][i] = E->X[1][node] * dxdy[2][1] + E->X[2][node] * dxdy[2][2] + E->X[3][node] * dxdy[2][3];
				xx[3][i] = E->X[1][node] * dxdy[3][1] + E->X[2][node] * dxdy[3][2] + E->X[3][node] * dxdy[3][3];
			}
			for(k = 1; k <= oned; k++)
			{					/* all of the vpoints */
				for(d = 1; d <= E->mesh.nsd - 1; d++)
					for(e = 1; e <= E->mesh.nsd - 1; e++)
						dxda[d][e] = 0.0;

				for(i = 1; i <= oned; i++)	/* nodes */
					for(d = 1; d <= E->mesh.nsd - 1; d++)
						for(e = 1; e <= E->mesh.nsd - 1; e++)
							dxda[d][e] += xx[e][i] * E->Mx.vpt[GMVXINDEX(d - 1, i, k)];

				jacobian = determinant(dxda, E->mesh.nsd - 1);
				dGammax->vpt[GMVGAMMA(ii, k)] = jacobian;


			}
		}
	}

	return;
}

void get_global_1d_shape_fn1(struct All_variables *E, int el, struct Shape_function1 *GM, struct Shape_function1_dA *dGammax, int top)
{
	int i, k, d, e, ii;
	//int dirn, locn, node[5];
	//int collapsed_dirn[2];
	//double scale[4];

	double jacobian;

	//double avet, aver, dxda[4][4], cof[4][4], xx[4][5];
	double avet, aver, dxda[4][4], xx[4][5];

	for(ii = 0; ii <= top; ii++)
	{							/* ii=0 for bottom ii=1 for top */

		if(E->control.CART3D)
		{
			for(i = 1; i <= onedvpoints[E->mesh.nsd]; i++)
			{					/* nodes */
				e = i + ii * onedvpoints[E->mesh.nsd];
				xx[1][i] = E->X[1][E->ien[el].node[e]];
				xx[2][i] = E->X[2][E->ien[el].node[e]];
			}
			for(k = 1; k <= onedvpoints[E->mesh.nsd]; k++)
			{					/* all of the vpoints */

				for(d = 1; d <= E->mesh.nsd - 1; d++)
					for(e = 1; e <= E->mesh.nsd - 1; e++)
						dxda[d][e] = 0.0;

				for(i = 1; i <= onedvpoints[E->mesh.nsd]; i++)	/* nodes */
					for(d = 1; d <= E->mesh.nsd - 1; d++)
						for(e = 1; e <= E->mesh.nsd - 1; e++)
							dxda[d][e] += xx[e][i] * E->Mx.vpt[GMVXINDEX(d - 1, i, k)];

				jacobian = determinant(dxda, E->mesh.nsd - 1);
				dGammax->vpt[GMVGAMMA(ii, k)] = jacobian;
			}
		}
		else if(E->control.Rsphere)
		{
			avet = 0;
			aver = 0;
			for(i = 1; i <= onedvpoints[E->mesh.nsd]; i++)
			{					/* nodes */
				e = i + ii * onedvpoints[E->mesh.nsd];
				xx[1][i] = E->SX[1][E->ien[el].node[e]];
				xx[2][i] = E->SX[2][E->ien[el].node[e]];
				avet += xx[1][i];
				aver += E->SX[3][E->ien[el].node[e]];
			}
			avet = avet / onedvpoints[E->mesh.nsd];
			aver = aver / onedvpoints[E->mesh.nsd];

			for(k = 1; k <= onedvpoints[E->mesh.nsd]; k++)
			{					/* all of the vpoints */

				for(d = 1; d <= E->mesh.nsd - 1; d++)
					for(e = 1; e <= E->mesh.nsd - 1; e++)
						dxda[d][e] = 0.0;

				for(i = 1; i <= onedvpoints[E->mesh.nsd]; i++)	/* nodes */
					for(d = 1; d <= E->mesh.nsd - 1; d++)
						for(e = 1; e <= E->mesh.nsd - 1; e++)
							dxda[d][e] += xx[e][i] * E->Mx.vpt[GMVXINDEX(d - 1, i, k)];

				jacobian = determinant(dxda, E->mesh.nsd - 1);
				dGammax->vpt[GMVGAMMA(ii, k)] = jacobian * aver * aver * sin(avet);
				dGammax->vpt[GMVGAMMA(ii, k)] = jacobian;
			}
		}
	}

	return;
}

/*  ==========================================
    construct the lumped mass matrix. The full
    matrix is the FE integration of the density 
    field. The lumped version is the diagonal
    matrix obtained by letting the shape function
    Na be delta(a,b)
    ========================================== */


void mass_matrix(struct All_variables *E)
{
	//int k, n[9], node, el, i, nint, e, lv;
	int n[9], node, el, i, nint, e, lv;
	double temp[9], area, centre[4], rtf[4][9];
	//float dx1, dx2, dx3, xlowmean, normlow, normhigh, xhighmean;
	float dx1, dx2, dx3;
	//struct Shape_function GN;
	//struct Shape_function_dA dOmega;
	//struct Shape_function_dx GNx;

	//const int ppts = ppoints[E->mesh.nsd];
	const int vpts = vpoints[E->mesh.nsd];

	/* ECO .size can also be defined here */

	for(lv = E->mesh.levmax; lv >= E->mesh.levmin; lv--)
	{

		if(E->control.CART3D)
			for(el = 1; el <= E->lmesh.NEL[lv]; el++)
			{
				get_global_shape_fn(E, el, 0, rtf, 0, lv);
				get_global_shape_fn(E, el, 2, rtf, 0, lv);
			}
		else if(E->control.Rsphere)
			for(el = 1; el <= E->lmesh.NEL[lv]; el++)
			{
				get_global_shape_fn(E, el, 2, rtf, 1, lv);
				get_global_shape_fn(E, el, 0, rtf, 1, lv);
			}

		for(node = 1; node <= E->lmesh.NNO[lv]; node++)
			E->MASS[lv][node] = 0.0;

		for(e = 1; e <= E->lmesh.NEL[lv]; e++)
		{

			for(node = 1; node <= enodes[E->mesh.nsd]; node++)
				n[node] = E->IEN[lv][e].node[node];

			area = centre[1] = centre[2] = centre[3] = 0.0;

			if(E->control.Rsphere)
			{
				for(i = 1; i <= E->mesh.nsd; i++)
				{
					for(node = 1; node <= enodes[E->mesh.nsd]; node++)
						centre[i] += E->XX[lv][i][E->IEN[lv][e].node[node]];

					centre[i] = centre[i] / enodes[E->mesh.nsd];
				}

				dx3 = sqrt(centre[1] * centre[1] + centre[2] * centre[2] + centre[3] * centre[3]);
				dx1 = acos(centre[3] / dx3);
				dx2 = myatan(centre[2], centre[1]);

				E->ECO[lv][e].centre[1] = dx1;
				E->ECO[lv][e].centre[2] = dx2;
				E->ECO[lv][e].centre[3] = dx3;

				dx1 = max(fabs(E->SXX[lv][1][n[3]] - E->SXX[lv][1][n[1]]), fabs(E->SXX[lv][1][n[2]] - E->SXX[lv][1][n[4]]));
				E->ECO[lv][e].size[1] = dx1 * E->ECO[lv][e].centre[3];

				dx2 = max(fabs(E->SXX[lv][2][n[3]] - E->SXX[lv][2][n[1]]), fabs(E->SXX[lv][2][n[2]] - E->SXX[lv][2][n[4]]));
				E->ECO[lv][e].size[2] = dx2 * E->ECO[lv][e].centre[3] * sin(E->ECO[lv][e].centre[1]);

				dx3 = 0.25 * (E->SXX[lv][3][n[5]] + E->SXX[lv][3][n[6]] + E->SXX[lv][3][n[7]] + E->SXX[lv][3][n[8]] - E->SXX[lv][3][n[1]] - E->SXX[lv][3][n[2]] - E->SXX[lv][3][n[3]] - E->SXX[lv][3][n[4]]);
				E->ECO[lv][e].size[3] = sqrt(dx3 * dx3);
			}
			else if(E->control.CART3D)
			{
				for(i = 1; i <= E->mesh.nsd; i++)
				{
					for(node = 1; node <= enodes[E->mesh.nsd]; node++)
						centre[i] += E->XX[lv][i][E->IEN[lv][e].node[node]];

					E->ECO[lv][e].centre[i] = centre[i] / enodes[E->mesh.nsd];
				}				/* end loop for dof */

				dx2 = 0.25 * (E->XX[lv][2][n[3]] + E->XX[lv][2][n[4]] + E->XX[lv][2][n[7]] + E->XX[lv][2][n[8]] - E->XX[lv][2][n[1]] - E->XX[lv][2][n[2]] - E->XX[lv][2][n[5]] - E->XX[lv][2][n[6]]);
				E->ECO[lv][e].size[2] = sqrt(dx2 * dx2);

				dx1 = 0.25 * (E->XX[lv][1][n[2]] + E->XX[lv][1][n[3]] + E->XX[lv][1][n[6]] + E->XX[lv][1][n[7]] - E->XX[lv][1][n[1]] - E->XX[lv][1][n[4]] - E->XX[lv][1][n[5]] - E->XX[lv][1][n[8]]);
				E->ECO[lv][e].size[1] = sqrt(dx1 * dx1);

				dx3 = 0.25 * (E->XX[lv][3][n[5]] + E->XX[lv][3][n[6]] + E->XX[lv][3][n[7]] + E->XX[lv][3][n[8]] - E->XX[lv][3][n[1]] - E->XX[lv][3][n[2]] - E->XX[lv][3][n[3]] - E->XX[lv][3][n[4]]);
				E->ECO[lv][e].size[3] = sqrt(dx3 * dx3);

			}



			for(nint = 1; nint <= vpts; nint++)
				area += g_point[nint].weight[E->mesh.nsd - 1] * E->GDA[lv][e].vpt[nint];
			E->ECO[lv][e].area = area;

			for(node = 1; node <= enodes[E->mesh.nsd]; node++)
			{
				temp[node] = 0.0;
				for(nint = 1; nint <= vpts; nint++)
					temp[node] += E->GDA[lv][e].vpt[nint] * g_point[nint].weight[E->mesh.nsd - 1] * E->N.vpt[GNVINDEX(node, nint)];	/* int Na dV */
			}

			for(node = 1; node <= enodes[E->mesh.nsd]; node++)
				E->MASS[lv][E->IEN[lv][e].node[node]] += temp[node];

			for(node = 1; node <= enodes[E->mesh.nsd]; node++)
			{
				E->TWW[lv][e].node[node] = temp[node];
			}

		}						/* over element */


		exchange_node_f20(E, E->MASS[lv], lv);


		for(node = 1; node <= E->lmesh.NNO[lv]; node++)
		{
			E->MASS[lv][node] = 1.0 / E->MASS[lv][node];
		}


	}


	E->lmesh.volume = 0;
	E->mesh.volume = 0;

	for(e = 1; e <= E->lmesh.nel; e++)
	{
		E->lmesh.volume += E->eco[e].area;
	}

	MPI_Allreduce(&E->lmesh.volume, &E->mesh.volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


	return;
}
