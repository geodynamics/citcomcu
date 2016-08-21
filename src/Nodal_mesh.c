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

/* Functions relating to the building and use of mesh locations ... */


#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

extern int Emergency_stop;

/* =================================================
   Standard node positions including mesh refinement 

   =================================================  */

void node_locations(struct All_variables *E)
{
	//int lev, nodel, i, j, k, ii, d, node;
	int lev, nodel, i, j, k, d;
	float *XX[4], *XG[4], dx[4], dxx[40];
	//int n00, nox, noz, noy, fn;
	int nox, noz, noy;
	double rad_conv;
	int m;
	//const int dims = E->mesh.nsd;

	m = E->parallel.me;

	rad_conv = M_PI / 180;

	nox = E->mesh.nox;
	noz = E->mesh.noz;
	noy = E->mesh.noy;

	for(d = 1; d <= E->mesh.nsd; d++)
	{
		XX[d] = (float *)safe_malloc((2 + E->mesh.nnx[d]) * sizeof(float));
		XG[d] = (float *)safe_malloc((2 + 1) * sizeof(float));
	}

	if(E->control.CART3D)
		for(d = 1; d <= E->mesh.nsd; d++)
		{						/* for even space */
			dx[d] = E->mesh.layer[d] / (E->mesh.nnx[d] - 1);
			XX[d][1] = 0.0;
			XX[d][E->mesh.nnx[d]] = E->mesh.layer[d];
			for(i = 2; i < E->mesh.nnx[d]; i++)
				XX[d][i] = XX[d][i - 1] + dx[d];
		}
	else if(E->control.Rsphere)
		for(d = 1; d <= E->mesh.nsd; d++)
		{						/* for even space */
			dx[d] = (E->sphere.corner[2][d] - E->sphere.corner[1][d]) / (E->mesh.nnx[d] - 1);
			XX[d][1] = E->sphere.corner[1][d];
			XX[d][E->mesh.nnx[d]] = E->sphere.corner[2][d];
			for(i = 2; i < E->mesh.nnx[d]; i++)
				XX[d][i] = XX[d][i - 1] + dx[d];
		}


	if(E->control.CART3D)
	{
		input_int("z_grid_layers", &(E->segment.zlayers), "1", m);
		input_float_vector("zz", E->segment.zlayers, (E->segment.zzlayer), m);
		input_int_vector("nz", E->segment.zlayers, (E->segment.nzlayer), m);

		input_int("x_grid_layers", &(E->segment.xlayers), "1", m);
		input_float_vector("xx", E->segment.xlayers, (E->segment.xxlayer), m);
		input_int_vector("nx", E->segment.xlayers, (E->segment.nxlayer), m);

		input_int("y_grid_layers", &(E->segment.ylayers), "1", m);
		input_float_vector("yy", E->segment.ylayers, (E->segment.yylayer), m);
		input_int_vector("ny", E->segment.ylayers, (E->segment.nylayer), m);
	}
	else if(E->control.Rsphere)
	{
		input_int("r_grid_layers", &(E->segment.zlayers), "1", m);
		input_float_vector("rr", E->segment.zlayers, (E->segment.zzlayer), m);
		input_int_vector("nr", E->segment.zlayers, (E->segment.nzlayer), m);

		input_int("t_grid_layers", &(E->segment.xlayers), "1", m);
		input_float_vector("tt", E->segment.xlayers, (E->segment.xxlayer), m);
		input_int_vector("nt", E->segment.xlayers, (E->segment.nxlayer), m);

		input_int("f_grid_layers", &(E->segment.ylayers), "1", m);
		input_float_vector("ff", E->segment.ylayers, (E->segment.yylayer), m);
		input_int_vector("nf", E->segment.ylayers, (E->segment.nylayer), m);
	}


/* 1st direction  */

	for(j = 1; j < E->segment.xlayers; j++)
		dxx[j] = (E->segment.xxlayer[j] - E->segment.xxlayer[j - 1]) / (E->segment.nxlayer[j] - E->segment.nxlayer[j - 1]);
	j = 1;
	for(i = 2; i < E->mesh.nnx[1]; i++)
	{
		if(i <= E->segment.nxlayer[j])
			XX[1][i] = XX[1][i - 1] + dxx[j];
		if(i == E->segment.nxlayer[j])
			j++;
	}


/* 3rd dimension vertical */

	for(j = 1; j < E->segment.zlayers; j++)
		dxx[j] = (E->segment.zzlayer[j] - E->segment.zzlayer[j - 1]) / (E->segment.nzlayer[j] - E->segment.nzlayer[j - 1]);
	j = 1;
	for(i = 2; i < E->mesh.nnx[3]; i++)
	{
		if(i <= E->segment.nzlayer[j])
			XX[3][i] = XX[3][i - 1] + dxx[j];
		if(i == E->segment.nzlayer[j])
			j++;
	}

/* 2nd dimension */

	for(j = 1; j < E->segment.ylayers; j++)
		dxx[j] = (E->segment.yylayer[j] - E->segment.yylayer[j - 1]) / (E->segment.nylayer[j] - E->segment.nylayer[j - 1]);
	j = 1;
	for(i = 2; i < E->mesh.nnx[2]; i++)
	{
		if(i <= E->segment.nylayer[j])
			XX[2][i] = XX[2][i - 1] + dxx[j];
		if(i == E->segment.nylayer[j])
			j++;
	}



	lev = E->mesh.levmax;
	nox = E->lmesh.NOX[lev];
	noy = E->lmesh.NOY[lev];
	noz = E->lmesh.NOZ[lev];


	if(E->control.CART3D)
	{
		for(k = 1; k <= noy; k++)
			for(i = 1; i <= nox; i++)
				for(j = 1; j <= noz; j++)
				{
					nodel = j + (i - 1) * noz + (k - 1) * noz * nox;
					E->XX[lev][1][nodel] = XX[1][i - 1 + E->lmesh.NXS[lev]];
					E->XX[lev][2][nodel] = XX[2][k - 1 + E->lmesh.NYS[lev]];
					E->XX[lev][3][nodel] = XX[3][j - 1 + E->lmesh.NZS[lev]];
				}
	}
	else if(E->control.Rsphere)
	{
		for(k = 1; k <= noy; k++)
			for(i = 1; i <= nox; i++)
				for(j = 1; j <= noz; j++)
				{
					nodel = j + (i - 1) * noz + (k - 1) * noz * nox;
					E->SXX[lev][1][nodel] = XX[1][i - 1 + E->lmesh.NXS[lev]] * rad_conv;
					E->SXX[lev][2][nodel] = XX[2][k - 1 + E->lmesh.NYS[lev]] * rad_conv;
					E->SXX[lev][3][nodel] = XX[3][j - 1 + E->lmesh.NZS[lev]];
					E->XX[lev][1][nodel] = E->SXX[lev][3][nodel] * sin(E->SXX[lev][1][nodel]) * cos(E->SXX[lev][2][nodel]);
					E->XX[lev][2][nodel] = E->SXX[lev][3][nodel] * sin(E->SXX[lev][1][nodel]) * sin(E->SXX[lev][2][nodel]);
					E->XX[lev][3][nodel] = E->SXX[lev][3][nodel] * cos(E->SXX[lev][1][nodel]);
				}
	}
	else
	{
		fprintf(stderr, "Your requested geometry is not supported! Sorry!\n");
		parallel_process_termination();
	}


	if(E->control.CART3D)
	{
		for(j = 1; j <= E->lmesh.noz; j++)
			E->XP[3][j] = E->XX[E->mesh.levmax][3][j];
		for(k = 1; k <= E->lmesh.noy; k++)
			E->XP[2][k] = E->XX[E->mesh.levmax][2][1 + (k - 1) * E->lmesh.noz * E->lmesh.nox];
		for(i = 1; i <= E->lmesh.nox; i++)
			E->XP[1][i] = E->XX[E->mesh.levmax][1][1 + (i - 1) * E->lmesh.noz];
		E->XG1[1] = CITCOM_TRACER_EPS_MARGIN;
		E->XG1[2] = CITCOM_TRACER_EPS_MARGIN;
		E->XG1[3] = CITCOM_TRACER_EPS_MARGIN;
		E->XG2[1] = E->mesh.layer[1] - CITCOM_TRACER_EPS_MARGIN;
		E->XG2[2] = E->mesh.layer[2] - CITCOM_TRACER_EPS_MARGIN;
		E->XG2[3] = E->mesh.layer[3] - CITCOM_TRACER_EPS_MARGIN;
	}
	else if(E->control.Rsphere)
	{
		for(j = 1; j <= E->lmesh.noz; j++)
			E->XP[3][j] = E->SXX[E->mesh.levmax][3][j];
		for(k = 1; k <= E->lmesh.noy; k++)
			E->XP[2][k] = E->SXX[E->mesh.levmax][2][1 + (k - 1) * E->lmesh.noz * E->lmesh.nox];
		for(i = 1; i <= E->lmesh.nox; i++)
			E->XP[1][i] = E->SXX[E->mesh.levmax][1][1 + (i - 1) * E->lmesh.noz];


		E->XG1[1] = E->sphere.ti * rad_conv + CITCOM_TRACER_EPS_MARGIN;
		E->XG1[2] = E->sphere.fi * rad_conv + CITCOM_TRACER_EPS_MARGIN;
		E->XG1[3] = E->sphere.ri + CITCOM_TRACER_EPS_MARGIN;
		E->XG2[1] = E->sphere.to * rad_conv - CITCOM_TRACER_EPS_MARGIN;
		E->XG2[2] = E->sphere.fo * rad_conv - CITCOM_TRACER_EPS_MARGIN;
		E->XG2[3] = E->sphere.ro - CITCOM_TRACER_EPS_MARGIN;
	}

	pre_interpolation(E);


	if(E->control.CART3D)
	{
		for(lev = E->mesh.levmax; lev > E->mesh.levmin; lev--)
			inject_node_fvector(E, lev, E->XX[lev], E->XX[lev - 1]);
	}
	else if(E->control.Rsphere)
	{
		for(lev = E->mesh.levmax; lev > E->mesh.levmin; lev--)
		{
			inject_node_fvector(E, lev, E->XX[lev], E->XX[lev - 1]);
			inject_node_fvector(E, lev, E->SXX[lev], E->SXX[lev - 1]);
		}
	}


	if(E->control.verbose)
	{
		for(lev = E->mesh.levmax; lev >= E->mesh.levmin; lev--)
		{
			fprintf(E->fp, "output_coordinates \n");
			if(E->control.Rsphere)
				for(i = 1; i <= E->lmesh.NNO[lev]; i++)
					fprintf(E->fp, "%d %g %g %g\n", i, E->SXX[lev][1][i], E->SXX[lev][2][i], E->SXX[lev][3][i]);
			else
				for(i = 1; i <= E->lmesh.NNO[lev]; i++)
					fprintf(E->fp, "%d %g %g %g\n", i, E->XX[lev][1][i], E->XX[lev][2][i], E->XX[lev][3][i]);
		}
	}

	for(d = 1; d <= E->mesh.nsd; d++)
	{
		free((void *)XX[d]);
		free((void *)XG[d]);
	}

	return;
}

void pre_interpolation(struct All_variables *E)
{
	//int i, j, k, e;
	int j, e;

	for(j = 1; j <= E->lmesh.rnoz; j++)
		E->XRG[3][j] = E->XP[3][1] + (j - 1) * (E->XP[3][E->lmesh.noz] - E->XP[3][1]) / (E->lmesh.rnoz - 1);

	E->XRG[3][E->lmesh.rnoz] = E->XP[3][E->lmesh.noz];
	E->XRG[3][1] = E->XP[3][1];

	for(j = 1; j < E->lmesh.rnoz; j++)
	{
		E->RG[3][j] = 0;
		for(e = 1; e < E->lmesh.noz; e++)
		{
			if(E->XRG[3][j + 1] <= E->XP[3][e + 1] && E->XRG[3][j] >= E->XP[3][e])
			{
				E->RG[3][j] = e;
				break;
			}
		}
//      fprintf(E->fp,"aaa %d %d %g\n",j,E->RG[3][j],E->XRG[3][j]);
	}

	return;
}



void dlogical_mesh_to_real(struct All_variables *E, double *data, int level)

{
	int i, j, n1, n2;

	if(E->mesh.periodic_x)
		for(i = 1; i <= E->mesh.NOZ[level]; i++)
			for(j = 1; j <= E->mesh.NOY[level]; j++)
			{
				n1 = i + (j - 1) * E->mesh.NOX[level] * E->mesh.NOZ[level];
				n2 = n1 + (E->mesh.NOX[level] - 1) * E->mesh.NOZ[level];

				data[n2] = data[n1];

			}

	if(E->mesh.periodic_y)
		for(i = 1; i <= E->mesh.NOZ[level]; i++)
			for(j = 1; j <= E->mesh.NOX[level]; j++)
			{
				n1 = i + (j - 1) * E->mesh.NOZ[level];
				n2 = n1 + (E->mesh.NOY[level] - 1) * E->mesh.NOZ[level] * E->mesh.NOX[level];

				data[n2] = data[n1];

			}

	if(E->mesh.periodic_y && E->mesh.periodic_x)	/* then need to do the 1st one again */
		for(i = 1; i <= E->mesh.NOZ[level]; i++)
			for(j = 1; j <= E->mesh.NOY[level]; j++)
			{
				n1 = i + (j - 1) * E->mesh.NOX[level] * E->mesh.NOZ[level];
				n2 = n1 + (E->mesh.NOX[level] - 1) * E->mesh.NOZ[level];

				data[n2] = data[n1];

			}

	return;
}


void flogical_mesh_to_real(struct All_variables *E, float *data, int level)
{
	int i, j, n1, n2;

	if(E->mesh.periodic_x)
		for(i = 1; i <= E->mesh.NOZ[level]; i++)
			for(j = 1; j <= E->mesh.NOY[level]; j++)
			{
				n1 = i + (j - 1) * E->mesh.NOX[level] * E->mesh.NOZ[level];
				n2 = n1 + (E->mesh.NOX[level] - 1) * E->mesh.NOZ[level];

				data[n2] = data[n1];

			}

	if(E->mesh.periodic_y)
		for(i = 1; i <= E->mesh.NOZ[level]; i++)
			for(j = 1; j <= E->mesh.NOX[level]; j++)
			{
				n1 = i + (j - 1) * E->mesh.NOZ[level];
				n2 = n1 + (E->mesh.NOY[level] - 1) * E->mesh.NOZ[level] * E->mesh.NOX[level];

				data[n2] = data[n1];

			}

	if(E->mesh.periodic_y && E->mesh.periodic_x)	/* then need to do the 1st one again */
		for(i = 1; i <= E->mesh.NOZ[level]; i++)
			for(j = 1; j <= E->mesh.NOY[level]; j++)
			{
				n1 = i + (j - 1) * E->mesh.NOX[level] * E->mesh.NOZ[level];
				n2 = n1 + (E->mesh.NOX[level] - 1) * E->mesh.NOZ[level];

				data[n2] = data[n1];

			}

	return;
}

void p_to_nodes(struct All_variables *E, double *P, float *PN, int lev)
{
	//int e, element, node, j;
	int element, node, j;

	for(node = 1; node <= E->lmesh.NNO[lev]; node++)
		PN[node] = 0.0;

	for(element = 1; element <= E->lmesh.NEL[lev]; element++)
	{

		for(j = 1; j <= enodes[E->mesh.nsd]; j++)
		{
			node = E->IEN[lev][element].node[j];
			PN[node] += P[element] * E->TW[lev][node];
		}

	}

	exchange_node_f20(E, PN, lev);

	return;
}
void e2_to_nodes(struct All_variables *E, float *e2, float *e2n, int lev)
{
  int element, node, j;
  
  for(node = 1; node <= E->lmesh.NNO[lev]; node++)
    e2n[node] = 0.0;
  for(element = 1; element <= E->lmesh.NEL[lev]; element++){
    for(j = 1; j <= enodes[E->mesh.nsd]; j++){
      node = E->IEN[lev][element].node[j];
      e2n[node] += e2[element] * E->TW[lev][node];
    }
  }
  
  exchange_node_f20(E, e2n, lev);
  
  return;
}



void p_to_centres(struct All_variables *E, float *PN, double *P, int lev)
{
	//int p, element, node, j;
	int p, j;
	double weight;

	for(p = 1; p <= E->lmesh.NEL[lev]; p++)
		P[p] = 0.0;

	weight = 1.0 / ((double)enodes[E->mesh.nsd]);

	for(p = 1; p <= E->lmesh.NEL[lev]; p++)
		for(j = 1; j <= enodes[E->mesh.nsd]; j++)
			P[p] += PN[E->IEN[lev][p].node[j]] * weight;

	return;
}


void v_to_intpts(struct All_variables *E, float *VN, float *VE, int lev)
{
	//int e, i, j, k;
	int e, i, j;
	const int nsd = E->mesh.nsd;
	const int vpts = vpoints[nsd];
	const int ends = enodes[nsd];

	for(e = 1; e <= E->lmesh.NEL[lev]; e++)
		for(i = 1; i <= vpts; i++)
		{
			VE[(e - 1) * vpts + i] = 0.0;
			for(j = 1; j <= ends; j++)
				VE[(e - 1) * vpts + i] += VN[E->IEN[lev][e].node[j]] * E->N.vpt[GNVINDEX(j, i)];
		}
}

void v_to_nodes(struct All_variables *E, float *VE, float *VN, int lev)
{
	//int e, i, j, k, n;
	int e, i, j, n;
	const int nsd = E->mesh.nsd;
	const int vpts = vpoints[nsd];
	const int ends = enodes[nsd];
	for(i = 1; i <= E->lmesh.NNO[lev]; i++)
		VN[i] = 0.0;

	for(e = 1; e <= E->lmesh.NEL[lev]; e++)
		for(j = 1; j <= ends; j++)
		{
			n = E->IEN[lev][e].node[j];
			for(i = 1; i <= vpts; i++)
				VN[n] += E->N.vpt[GNVINDEX(j, i)] * E->TW[lev][n] * VE[(e - 1) * vpts + i];
		}
	flogical_mesh_to_real(E, VN, E->mesh.levmax);
	return;
}

void visc_to_intpts(struct All_variables *E, float *VN, float *VE, int lev)
{
	//int e, i, j, k;
	int e, i, j;
	const int nsd = E->mesh.nsd;
	const int vpts = vpoints[nsd];
	const int ends = enodes[nsd];

	for(e = 1; e <= E->lmesh.NEL[lev]; e++)
		for(i = 1; i <= vpts; i++)
		{
			VE[(e - 1) * vpts + i] = 0.0;
			for(j = 1; j <= ends; j++)
				VE[(e - 1) * vpts + i] += log(VN[E->IEN[lev][e].node[j]]) * E->N.vpt[GNVINDEX(j, i)];
			VE[(e - 1) * vpts + i] = exp(VE[(e - 1) * vpts + i]);
		}
	return;
}


void visc_to_nodes(struct All_variables *E, float *VE, float *VN, int lev)
{
	//int e, i, j, k, n;
	int e, i, j, n;
	const int nsd = E->mesh.nsd;
	const int vpts = vpoints[nsd];
	const int ends = enodes[nsd];
	double temp_visc;

	for(i = 1; i <= E->lmesh.NNO[lev]; i++)
		VN[i] = 0.0;

	for(e = 1; e <= E->lmesh.NEL[lev]; e++)
		for(j = 1; j <= ends; j++)
		{
			n = E->IEN[lev][e].node[j];
			temp_visc = 0.0;
			for(i = 1; i <= vpts; i++)
				temp_visc += E->TW[lev][n] * log(E->N.vpt[GNVINDEX(j, i)] * VE[(e - 1) * vpts + i]);
			VN[n] += exp(temp_visc);
		}
	
	return;
}

void visc_from_ele_to_gint(struct All_variables *E, float *VN, float *VE, int lev)
{
	//int m, e, i, j, k, n;
	int e, i;
	const int nsd = E->mesh.nsd;
	const int vpts = vpoints[nsd];
	//const int ends = enodes[nsd];

	for(e = 1; e <= E->lmesh.NEL[lev]; e++)
		for(i = 1; i <= vpts; i++)
		{
			VE[(e - 1) * vpts + i] = VN[e];
		}

	return;
}


void visc_from_gint_to_ele(struct All_variables *E, float *VE, float *VN, int lev)
{
	//int m, e, i, j, k, n;
	int e, i;
	const int nsd = E->mesh.nsd;
	const int vpts = vpoints[nsd];
	//const int ends = enodes[nsd];
	double temp_visc;

	for(e = 1; e <= E->lmesh.NEL[lev]; e++)
	{
		temp_visc = 0.0;
		for(i = 1; i <= vpts; i++)
			temp_visc += VE[(e - 1) * vpts + i];
		temp_visc = temp_visc / vpts;

		VN[e] = temp_visc;
	}

	return;
}



void visc_from_gint_to_nodes(struct All_variables *E, float *VE, float *VN, int lev)
{
	//int m, e, i, j, k, n;
	int e, i, j, n;
	const int nsd = E->mesh.nsd;
	const int vpts = vpoints[nsd];
	const int ends = enodes[nsd];
	double temp_visc;

	for(i = 1; i <= E->lmesh.NNO[lev]; i++)
		VN[i] = 0.0;

	for(e = 1; e <= E->lmesh.NEL[lev]; e++)
	{
		temp_visc = 0.0;
		for(i = 1; i <= vpts; i++)
			temp_visc += VE[(e - 1) * vpts + i];
		temp_visc = temp_visc / vpts;

		for(j = 1; j <= ends; j++)
		{
			n = E->IEN[lev][e].node[j];
			VN[n] += E->TWW[lev][e].node[j] * temp_visc;
		}
	}

	exchange_node_f20(E, VN, lev);

	for(n = 1; n <= E->lmesh.NNO[lev]; n++)
		VN[n] *= E->MASS[lev][n];

	return;
}

void visc_from_nodes_to_gint(struct All_variables *E, float *VN, float *VE, int lev)
{
	//int m, e, i, j, k, n;
	int e, i, j;
	const int nsd = E->mesh.nsd;
	const int vpts = vpoints[nsd];
	const int ends = enodes[nsd];
	double temp_visc;
	for(e = 1; e <= E->lmesh.NEL[lev]; e++)
		for(i = 1; i <= vpts; i++)
			VE[(e - 1) * vpts + i] = 0.0;

	for(e = 1; e <= E->lmesh.NEL[lev]; e++)
		for(i = 1; i <= vpts; i++)
		{
			temp_visc = 0.0;
			for(j = 1; j <= ends; j++)
				temp_visc += E->N.vpt[GNVINDEX(j, i)] * VN[E->IEN[lev][e].node[j]];

			VE[(e - 1) * vpts + i] = temp_visc;
		}

	return;
}
