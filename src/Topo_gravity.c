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

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

#define c_re(a) a.real
#define c_im(a) a.imag
typedef struct compl
{
	double real;
	double imag;
} COMPLEX;



/* ===================================================================
   Consistent boundary flux method for stress ... Zhong,Gurnis,Hulbert 

   Solve for the stress as the code defined it internally, rather than
   what was intended to be solved. This is more appropriate.

   Note also that the routine is dependent on the method 
   used to solve for the velocity in the first place.
   ===================================================================  */



/* call this only for top and bottom processors */
void get_CBF_topo(struct All_variables *E, float *H, float *HB)
{
	//int el, elb, els, node, nodeb, nodes, i, j, k, l, m, n, count;
	int el, elb, els, node, nodeb, i, l, m, n;
	//int nodel, nodem, nodesl, nodesm, lnsf, nel2;
	int lnsf;

	struct Shape_function1 GM, GMb;
	struct Shape_function1_dA dGammax, dGammabx;

	float *temp, *eltTU, *eltTL, *SU, *SL, *RU, *RL;

	double eltk[24 * 24], eltf[24];
	double eltkb[24 * 24], eltfb[24];
	double res[24], resb[24], eu[24], eub[24];
	higher_precision eltg[24][1], eltgb[24][1];

	const int dims = E->mesh.nsd;
	const int Tsize = 5;		/* maximum values, applicable to 3d, harmless for 2d */
	//const int Ssize = 4;
	const int ends = enodes[dims];
	//const int noz = E->lmesh.noz;
	//const int noy = E->lmesh.noy;
	const int nno = E->lmesh.nno;
	const int onedv = onedvpoints[dims];
	//const int snode1 = 1, snode2 = 4, snode3 = 5, snode4 = 8;
	const int elz = E->lmesh.elz;
	//const int ely = E->lmesh.ely;
	const int lev = E->mesh.levmax;

	lnsf = E->lmesh.nsf;

	eltTU = (float *)malloc((1 + Tsize) * sizeof(float));
	eltTL = (float *)malloc((1 + Tsize) * sizeof(float));
	SU = (float *)malloc((1 + lnsf) * sizeof(float));
	SL = (float *)malloc((1 + lnsf) * sizeof(float));
	RU = (float *)malloc((1 + lnsf) * sizeof(float));
	RL = (float *)malloc((1 + lnsf) * sizeof(float));
	temp = (float *)malloc((1 + nno) * sizeof(float));

	for(i = 0; i <= nno; i++)
		temp[i] = 0;

	for(i = 0; i <= lnsf; i++)
		RU[i] = RL[i] = SU[i] = SL[i] = 0.0;

	/* calculate the element residuals */

	for(els = 1; els <= E->lmesh.snel; els++)
	{
		el = E->surf_element[els];
		elb = el + elz - 1;

		for(m = 0; m < ends; m++)
		{						/* for bottom elements */
			nodeb = E->ien[elb].node[m + 1];
			eub[m * dims] = E->V[1][nodeb];
			eub[m * dims + 1] = E->V[2][nodeb];
			if(3 == dims)
				eub[m * dims + 2] = E->V[3][nodeb];
		}

		for(m = 0; m < ends; m++)
		{
			node = E->ien[el].node[m + 1];
			eu[m * dims] = E->V[1][node];
			eu[m * dims + 1] = E->V[2][node];
			if(3 == dims)
				eu[m * dims + 2] = E->V[3][node];
		}

		get_elt_k(E, el, eltk, lev, 1);
		get_elt_k(E, elb, eltkb, lev, 1);
		get_elt_f(E, el, eltf, 0, 0);
		get_elt_f(E, elb, eltfb, 0, 0);
		get_elt_g(E, el, eltg, lev);
		get_elt_g(E, elb, eltgb, lev);

		for(m = 0; m < dims * ends; m++)
		{
			res[m] = eltf[m] - E->elt_del[lev][el].g[m][0] * E->P[el];
			resb[m] = eltfb[m] - E->elt_del[lev][elb].g[m][0] * E->P[elb];
		}

		for(m = 0; m < dims * ends; m++)
			for(l = 0; l < dims * ends; l++)
			{
				res[m] -= eltk[ends * dims * m + l] * eu[l];
				resb[m] -= eltkb[ends * dims * m + l] * eub[l];
			}

		/* Put relevant (vertical & surface) parts of element residual into surface residual */

		for(m = 1; m <= ends; m++)
		{						/* for bottom elements */
			switch (m)
			{
			case 2:
				RL[E->sien[els].node[1]] += resb[(m - 1) * dims + 1];
				break;
			case 3:
				RL[E->sien[els].node[2]] += resb[(m - 1) * dims + 1];
				break;
			case 7:
				RL[E->sien[els].node[3]] += resb[(m - 1) * dims + 1];
				break;
			case 6:
				RL[E->sien[els].node[4]] += resb[(m - 1) * dims + 1];
				break;
			}
		}


		for(m = 1; m <= ends; m++)
		{
			switch (m)
			{
			case 1:
				RU[E->sien[els].node[1]] += res[(m - 1) * dims + 1];
				break;
			case 4:
				RU[E->sien[els].node[2]] += res[(m - 1) * dims + 1];
				break;
			case 8:
				RU[E->sien[els].node[3]] += res[(m - 1) * dims + 1];
				break;
			case 5:
				RU[E->sien[els].node[4]] += res[(m - 1) * dims + 1];
				break;
			}
		}
	}

	/* calculate the LHS */

	for(els = 1; els <= E->lmesh.snel; els++)
	{
		el = E->surf_element[els];
		elb = el + elz - 1;

		get_global_1d_shape_fn(E, el, &GM, &dGammax, 1);
		get_global_1d_shape_fn(E, elb, &GMb, &dGammabx, 1);

		for(m = 1; m <= onedv; m++)
		{
			eltTU[m - 1] = 0.0;
			eltTL[m - 1] = 0.0;
			for(n = 1; n <= onedv; n++)
			{
				eltTU[m - 1] += dGammax.vpt[GMVGAMMA(1, n)] * l_1d[n].weight[dims - 1] * E->L.vpt[GMVINDEX(m, n)] * E->L.vpt[GMVINDEX(m, n)];
				eltTL[m - 1] += dGammabx.vpt[GMVGAMMA(1 + dims, n)] * l_1d[n].weight[dims - 1] * E->L.vpt[GMVINDEX(m, n)] * E->L.vpt[GMVINDEX(m, n)];
			}
		}

		for(m = 1; m <= onedv; m++)	/* for bottom */
			SL[E->sien[els].node[m]] += eltTL[m - 1];

		for(m = 1; m <= onedv; m++)
			SU[E->sien[els].node[m]] += eltTU[m - 1];
	}

	for(i = 1; i <= E->lmesh.nsf; i++)
	{
		node = E->surf_node[i];
		temp[node] = RU[i];
	}
	exchange_node_f20(E, temp, E->mesh.levmax);
	for(i = 1; i <= E->lmesh.nsf; i++)
	{
		node = E->surf_node[i];
		RU[i] = temp[node];
	}
	for(i = 1; i <= E->lmesh.nsf; i++)
	{
		node = E->surf_node[i];
		temp[node] = SU[i];
	}
	exchange_node_f20(E, temp, E->mesh.levmax);
	for(i = 1; i <= E->lmesh.nsf; i++)
	{
		node = E->surf_node[i];
		SU[i] = temp[node];
	}

	for(i = 1; i <= E->lmesh.nsf; i++)
	{
		node = E->surf_node[i];
		temp[node] = RL[i];
	}
	exchange_node_f20(E, temp, E->mesh.levmax);
	for(i = 1; i <= E->lmesh.nsf; i++)
	{
		node = E->surf_node[i];
		RL[i] = temp[node];
	}
	for(i = 1; i <= E->lmesh.nsf; i++)
	{
		node = E->surf_node[i];
		temp[node] = SL[i];
	}
	exchange_node_f20(E, temp, E->mesh.levmax);
	for(i = 1; i <= E->lmesh.nsf; i++)
	{
		node = E->surf_node[i];
		SL[i] = temp[node];
	}

	if(E->parallel.me_loc[3] == 0)
	{
		for(i = 1; i <= E->lmesh.nsf; i++)
			H[i] = -RU[i] / SU[i];
	}

	if(E->parallel.me_loc[3] == E->parallel.nprocz - 1)
	{
		for(i = 1; i <= E->lmesh.nsf; i++)
			HB[i] = -RL[i] / SL[i];
	}

	free((void *)eltTU);
	free((void *)eltTL);
	free((void *)SU);
	free((void *)SL);
	free((void *)RU);
	free((void *)RL);
	free((void *)temp);
	return;

}

void get_STD_topo(struct All_variables *E, float *tpg, float *tpgb, int ii)
{
	//int i, j, k, e, nel2, snode, node;
	int i, j, e, snode, node;

	float *SZZ, *SXX, *SYY, *SXY, *SXZ, *SZY, VZ[9], VY[9], VX[9], Szz, Sxx, Syy, Sxy, Sxz, Szy;
	float Vzz[9], Vxx[9], Vyy[9], Vxy[9], Vxz[9], Vzy[9];
	//float pre[9], el_volume, tww[9], Visc, a, b;
	float pre[9];
	double rtf[4][9];

	const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int vpts = vpoints[dims];
	//const int ppts = ppoints[dims];
	const int ends = enodes[dims];
	const int nno = E->lmesh.nno;
	const int nel = E->lmesh.nel;
	const int lev = E->mesh.levmax;

	SXX = (float *)malloc((nno + 1) * sizeof(float));
	SYY = (float *)malloc((nno + 1) * sizeof(float));
	SXY = (float *)malloc((nno + 1) * sizeof(float));
	SXZ = (float *)malloc((nno + 1) * sizeof(float));
	SZY = (float *)malloc((nno + 1) * sizeof(float));
	SZZ = (float *)malloc((nno + 1) * sizeof(float));

	for(i = 1; i <= nno; i++)
	{
		SZZ[i] = 0.0;
		SXX[i] = 0.0;
		SYY[i] = 0.0;
		SXY[i] = 0.0;
		SXZ[i] = 0.0;
		SZY[i] = 0.0;
	}

	for(e = 1; e <= nel; e++)
	{
		Szz = 0.0;
		Sxx = 0.0;
		Syy = 0.0;
		Sxy = 0.0;
		Sxz = 0.0;
		Szy = 0.0;
		get_rtf(E, e, 0, rtf, lev);

		for(j = 1; j <= vpts; j++)
		{
			pre[j] = E->EVI[lev][(e - 1) * vpts + j] * E->gDA[e].vpt[j];
			Vzz[j] = 0.0;
			Vxx[j] = 0.0;
			Vyy[j] = 0.0;
			Vxy[j] = 0.0;
			Vxz[j] = 0.0;
			Vzy[j] = 0.0;
		}

		for(j = 1; j <= ends; j++)
		{
			VX[j] = E->V[1][E->ien[e].node[j]];
			VZ[j] = E->V[3][E->ien[e].node[j]];
			VY[j] = E->V[2][E->ien[e].node[j]];
		}

		for(i = 1; i <= vpts; i++)
		{
			if(E->control.CART3D)
				for(j = 1; j <= ends; j++)
				{
					Vzz[i] += VZ[j] * E->gNX[e].vpt[GNVXINDEX(2, j, i)];
					Vxx[i] += VX[j] * E->gNX[e].vpt[GNVXINDEX(0, j, i)];
					Vxz[i] += (VX[j] * E->gNX[e].vpt[GNVXINDEX(2, j, i)] + VZ[j] * E->gNX[e].vpt[GNVXINDEX(0, j, i)]);
					Vyy[i] += VY[j] * E->gNX[e].vpt[GNVXINDEX(1, j, i)];
					Vxy[i] += (VX[j] * E->gNX[e].vpt[GNVXINDEX(1, j, i)] + VY[j] * E->gNX[e].vpt[GNVXINDEX(0, j, i)]);
					Vzy[i] += (VY[j] * E->gNX[e].vpt[GNVXINDEX(2, j, i)] + VZ[j] * E->gNX[e].vpt[GNVXINDEX(1, j, i)]);
				}
			else if(E->control.Rsphere)
				for(j = 1; j <= ends; j++)
				{
					Vzz[i] += VZ[j] * E->gNX[e].vpt[GNVXINDEX(2, j, i)];
					Vxx[i] += (VX[j] * E->gNX[e].vpt[GNVXINDEX(0, j, i)] + VZ[j] * E->N.vpt[GNVINDEX(j, i)]) * rtf[3][i];
					Vxz[i] += VX[j] * E->gNX[e].vpt[GNVXINDEX(2, j, i)] + rtf[3][i] * (VZ[j] * E->gNX[e].vpt[GNVXINDEX(0, j, i)] - VX[j] * E->N.vpt[GNVINDEX(j, i)]);
					Vyy[i] += ((VY[j] * E->gNX[e].vpt[GNVXINDEX(1, j, i)] + VX[j] * E->N.vpt[GNVINDEX(j, i)] * cos(rtf[1][i])) / sin(rtf[1][i]) + VZ[j] * E->N.vpt[GNVINDEX(j, i)]) * rtf[3][i];
					Vxy[i] += ((VX[j] * E->gNX[e].vpt[GNVXINDEX(1, j, i)] - VY[j] * E->N.vpt[GNVINDEX(j, i)] * cos(rtf[1][i])) / sin(rtf[1][i]) + VY[j] * E->gNX[e].vpt[GNVXINDEX(0, j, i)]) * rtf[3][i];
					Vzy[i] += VY[j] * E->gNX[e].vpt[GNVXINDEX(2, j, i)] + rtf[3][i] * (VZ[j] * E->gNX[e].vpt[GNVXINDEX(1, j, i)] / sin(rtf[1][i]) - VY[j] * E->N.vpt[GNVINDEX(j, i)]);

				}

			Sxx += 2.0 * pre[i] * Vxx[i];
			Syy += 2.0 * pre[i] * Vyy[i];
			Szz += 2.0 * pre[i] * Vzz[i];
			Sxy += pre[i] * Vxy[i];
			Sxz += pre[i] * Vxz[i];
			Szy += pre[i] * Vzy[i];
		}


		Sxx /= E->eco[e].area;
		Syy /= E->eco[e].area;
		Szz /= E->eco[e].area;
		Sxz /= E->eco[e].area;
		Sxy /= E->eco[e].area;
		Szy /= E->eco[e].area;

		Szz -= E->P[e];			/* add the pressure term */
		Sxx -= E->P[e];			/* add the pressure term */
		Syy -= E->P[e];			/* add the pressure term */

		for(j = 1; j <= ends; j++)
		{
			node = E->ien[e].node[j];
			SZZ[node] += E->TWW[E->mesh.levmax][e].node[j] * Szz;
			SXX[node] += E->TWW[E->mesh.levmax][e].node[j] * Sxx;
			SXZ[node] += E->TWW[E->mesh.levmax][e].node[j] * Sxz;
			SYY[node] += E->TWW[E->mesh.levmax][e].node[j] * Syy;
			SXY[node] += E->TWW[E->mesh.levmax][e].node[j] * Sxy;
			SZY[node] += E->TWW[E->mesh.levmax][e].node[j] * Szy;
		}

	}

	exchange_node_f20(E, SXX, lev);
	exchange_node_f20(E, SZZ, lev);
	exchange_node_f20(E, SXZ, lev);
	exchange_node_f20(E, SXY, lev);
	exchange_node_f20(E, SYY, lev);
	exchange_node_f20(E, SZY, lev);

	for(i = 1; i <= nno; i++)
	{
		SZZ[i] = SZZ[i] * E->Mass[i];
		SXX[i] = SXX[i] * E->Mass[i];
		SXZ[i] = SXZ[i] * E->Mass[i];
		SYY[i] = SYY[i] * E->Mass[i];
		SXY[i] = SXY[i] * E->Mass[i];
		SZY[i] = SZY[i] * E->Mass[i];
	}

	for(snode = 1; snode <= E->lmesh.nsf; snode++)
	{
		node = E->surf_node[snode];
		tpg[snode] = -2 * SZZ[node] + SZZ[node - 1];
		tpgb[snode] = 2 * SZZ[node - E->lmesh.noz + 1] - SZZ[node - E->lmesh.noz + 2];
	}

	free((void *)SXX);
	free((void *)SYY);
	free((void *)SXY);
	free((void *)SXZ);
	free((void *)SZY);
	free((void *)SZZ);

	return;
}
