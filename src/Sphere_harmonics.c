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
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include <stdlib.h>

/*   ======================================================================
    ======================================================================  */

void set_sphere_harmonics(struct All_variables *E)
{
	int node, ll, mm, i, j;
	double dth, dfi;

	for(i = 0; i <= E->sphere.llmax; i++)
		E->sphere.hindex[i] = (int *)malloc((E->sphere.llmax + 3) * sizeof(int));

	E->sphere.elx = E->sphere.nox - 1;
	E->sphere.ely = E->sphere.noy - 1;
	E->sphere.snel = E->sphere.ely * E->sphere.elx;
	E->sphere.nsf = E->sphere.noy * E->sphere.nox;

	E->sphere.lelx = E->sphere.elx / E->parallel.nprocx;
	E->sphere.lely = E->sphere.elx / E->parallel.nprocy;
	E->sphere.lsnel = E->sphere.lely * E->sphere.lelx;
	E->sphere.lnox = E->sphere.lelx + 1;
	E->sphere.lnoy = E->sphere.lely + 1;
	E->sphere.lnsf = E->sphere.lnox * E->sphere.lnoy;
	E->sphere.lexs = E->sphere.lelx * E->parallel.nprocx;
	E->sphere.leys = E->sphere.lely * E->parallel.nprocy;


	E->sphere.sx[1] = (double *)malloc((E->sphere.nsf + 1) * sizeof(double));
	E->sphere.sx[2] = (double *)malloc((E->sphere.nsf + 1) * sizeof(double));

	i = 0;
	for(ll = 0; ll <= E->sphere.llmax; ll++)
		for(mm = 0; mm <= ll; mm++)
		{
			E->sphere.hindex[ll][mm] = i;
			i++;
		}
	E->sphere.hindice = i;

	E->sphere.con = (double *)malloc((E->sphere.hindice + 3) * sizeof(double));

	for(i = 1; i <= E->sphere.lelx; i++)
		E->sphere.tableplm[i] = (double *)malloc((E->sphere.hindice + 3) * sizeof(double));
	for(i = 1; i <= E->sphere.lely; i++)
	{
		E->sphere.tablecosf[i] = (double *)malloc((E->sphere.llmax + 3) * sizeof(double));
		E->sphere.tablesinf[i] = (double *)malloc((E->sphere.llmax + 3) * sizeof(double));
	}

	for(i = 1; i <= E->sphere.lnox; i++)
		E->sphere.tableplm_n[i] = (double *)malloc((E->sphere.hindice + 3) * sizeof(double));
	for(i = 1; i <= E->sphere.lnoy; i++)
	{
		E->sphere.tablecosf_n[i] = (double *)malloc((E->sphere.llmax + 3) * sizeof(double));
		E->sphere.tablesinf_n[i] = (double *)malloc((E->sphere.llmax + 3) * sizeof(double));
	}
	E->sphere.sien = (struct SIEN *)malloc((E->sphere.lsnel + 1) * sizeof(struct SIEN));



	for(i = 0; i <= 1; i++)
	{
		E->sphere.harm_tpgt[i] = (float *)malloc((E->sphere.hindice + 2) * sizeof(float));
		E->sphere.harm_tpgb[i] = (float *)malloc((E->sphere.hindice + 2) * sizeof(float));
		E->sphere.harm_velp[i] = (float *)malloc((E->sphere.hindice + 2) * sizeof(float));
		E->sphere.harm_velt[i] = (float *)malloc((E->sphere.hindice + 2) * sizeof(float));
		E->sphere.harm_divg[i] = (float *)malloc((E->sphere.hindice + 2) * sizeof(float));
		E->sphere.harm_vort[i] = (float *)malloc((E->sphere.hindice + 2) * sizeof(float));
		E->sphere.harm_visc[i] = (float *)malloc((E->sphere.hindice + 2) * sizeof(float));
		E->sphere.harm_geoid[i] = (float *)malloc((E->sphere.hindice + 2) * sizeof(float));
	}



	for(ll = 0; ll <= E->sphere.llmax; ll++)
		for(mm = 0; mm <= ll; mm++)
		{
			E->sphere.con[E->sphere.hindex[ll][mm]] = sqrt((2.0 - ((mm == 0) ? 1.0 : 0.0)) * (2 * ll + 1) / (4.0 * M_PI)) * sqrt_multis(ll + mm, ll - mm);	/* which is sqrt((ll-mm)!/(ll+mm)!) */
		}

	dth = M_PI / E->sphere.elx;
	dfi = 2.0 * M_PI / E->sphere.ely;

	for(j = 1; j <= E->sphere.noy; j++)
		for(i = 1; i <= E->sphere.nox; i++)
		{
			node = i + (j - 1) * E->sphere.nox;
			E->sphere.sx[1][node] = dth * (i - 1);
			E->sphere.sx[2][node] = dfi * (j - 1);
		}



	compute_sphereh_table(E);

	return;
}

/*   ======================================================================
    ======================================================================  */

void sphere_harmonics_layer(struct All_variables *E, float **T, float *sphc, float *sphs, int iprint, char *filen)
{
	//FILE *fp;
	//char output_file[255];
	//int i, node, j, ll, mm, printt, proc_loc;
	int printt, proc_loc;
	//float minx, maxx, t, f, rad;
	float minx, maxx, rad;
	//static int been_here = 0;
	float *TG;

	rad = 180.0 / M_PI;

	maxx = -1.e6;
	minx = 1.e6;

	printt = 0;

	if(E->parallel.me_loc[3] == E->parallel.nprocz - 1 && iprint == 1)
		printt = 1;
	if(E->parallel.me_loc[3] == 0 && iprint == 0)
		printt = 1;

	TG = (float *)malloc((E->sphere.nsf + 1) * sizeof(float));

	proc_loc = E->parallel.me_loc[3];

	sphere_interpolate(E, T, TG);

	sphere_expansion(E, TG, sphc, sphs);

	if(printt)
		print_field_spectral_regular(E, TG, sphc, sphs, proc_loc, filen);


	parallel_process_sync();

	free((void *)TG);

	return;
}

/* ===================================================================
  do the interpolation on sphere for data T, which is needed for both
  spherical harmonic expansion and graphics 
 =================================================================== */

void sphere_interpolate(struct All_variables *E, float **T, float *TG)
{
	//int ii, jj, es, i, j, m, el, node;
	int i, j, node;
	//double x[4], t, f;

	//const int ends = 4;

	TG[0] = 0.0;
	for(i = 1; i <= E->sphere.nox; i++)
		for(j = 1; j <= E->sphere.noy; j++)
		{
			node = i + (j - 1) * E->sphere.nox;
			TG[node] = 0.0;
			/* first find which cap this node (i,j) is in  */

		}						/* end for i and j */

	gather_TG_to_me0(E, TG);

	return;
}

/* =========================================================
  ========================================================= */
void sphere_expansion(struct All_variables *E, float *TG, float *sphc, float *sphs)
{
	int p, i, j, es, mm, ll;
	//double temp, area, t, f;
	double temp, area;
	const double pt25 = 0.25;
	//static int been_here = 0;

	for(i = 0; i < E->sphere.hindice; i++)
	{
		sphc[i] = 0.0;
		sphs[i] = 0.0;
	}

	area = 2.0 * M_PI * M_PI / (E->sphere.elx * E->sphere.ely);

	for(ll = 0; ll <= E->sphere.llmax; ll++)
		for(mm = 0; mm <= ll; mm++)
		{

			p = E->sphere.hindex[ll][mm];

			for(j = 1; j <= E->sphere.lely; j++)
				for(i = 1; i <= E->sphere.lelx; i++)
				{

					es = i + (j - 1) * E->sphere.lelx;

					temp = pt25 * (TG[E->sphere.sien[es].node[1]] + TG[E->sphere.sien[es].node[2]] + TG[E->sphere.sien[es].node[3]] + TG[E->sphere.sien[es].node[4]]);

					sphc[p] += temp * E->sphere.tableplm[i][p] * E->sphere.tablecosf[j][mm] * E->sphere.tablesint[i];
					sphs[p] += temp * E->sphere.tableplm[i][p] * E->sphere.tablesinf[j][mm] * E->sphere.tablesint[i];

				}

			sphc[p] *= area;
			sphs[p] *= area;

		}						/* end for ll and mm  */


	sum_across_surf_sph1(E, sphc, sphs);


	return;
}

 /* =========================================================== */
void inv_sphere_harmonics(struct All_variables *E, float *sphc, float *sphs, float *TG, int proc_loc)
{
	//int k, ll, mm, node, i, j, p, noz, snode;
	int ll, mm, node, i, j, p;
	//float t1, f1, rad;

	if(E->parallel.me_loc[3] == proc_loc)
	{

		for(j = 1; j <= E->sphere.noy; j++)
			for(i = 1; i <= E->sphere.nox; i++)
			{
				node = i + (j - 1) * E->sphere.nox;
				TG[node] = 0.0;
			}

		for(ll = 0; ll <= E->sphere.llmax; ll++)
			for(mm = 0; mm <= ll; mm++)
			{

				p = E->sphere.hindex[ll][mm];

				for(i = 1; i <= E->sphere.lnox; i++)
					for(j = 1; j <= E->sphere.lnoy; j++)
					{
						node = i + E->sphere.lexs + (j + E->sphere.leys - 1) * E->sphere.nox;
						TG[node] += (sphc[p] * E->sphere.tableplm_n[i][p] * E->sphere.tablecosf_n[j][mm] + sphs[p] * E->sphere.tableplm_n[i][p] * E->sphere.tablesinf_n[j][mm]);
					}
			}

		gather_TG_to_me0(E, TG);

	}

	parallel_process_sync();

	return;
}

/* ==================================================*/
void compute_sphereh_table(struct All_variables *E)
{
	int rr, node, ends, ll, mm, es, i, j, p;
	double t, f;
	const double pt25 = 0.25;

	ends = 4;

	for(j = 1; j <= E->sphere.lely; j++)
		for(i = 1; i <= E->sphere.lelx; i++)
		{
			es = i + (j - 1) * E->sphere.lelx;
			node = E->sphere.lexs + i + (E->sphere.leys + j - 1) * E->sphere.nox;
			for(rr = 1; rr <= ends; rr++)
				E->sphere.sien[es].node[rr] = node + offset[rr].vector[1] + offset[rr].vector[2] * E->sphere.nox;
		}

	for(j = 1; j <= E->sphere.lely; j++)
	{
		es = 1 + (j - 1) * E->sphere.lelx;
		f = pt25 * (E->sphere.sx[2][E->sphere.sien[es].node[1]] + E->sphere.sx[2][E->sphere.sien[es].node[2]] + E->sphere.sx[2][E->sphere.sien[es].node[3]] + E->sphere.sx[2][E->sphere.sien[es].node[4]]);
		for(mm = 0; mm <= E->sphere.llmax; mm++)
		{
			E->sphere.tablecosf[j][mm] = cos((double)(mm) * f);
			E->sphere.tablesinf[j][mm] = sin((double)(mm) * f);
		}
	}

	for(i = 1; i <= E->sphere.lelx; i++)
	{
		es = i + (1 - 1) * E->sphere.lelx;
		t = pt25 * (E->sphere.sx[1][E->sphere.sien[es].node[1]] + E->sphere.sx[1][E->sphere.sien[es].node[2]] + E->sphere.sx[1][E->sphere.sien[es].node[3]] + E->sphere.sx[1][E->sphere.sien[es].node[4]]);
		E->sphere.tablesint[i] = sin(t);
		for(ll = 0; ll <= E->sphere.llmax; ll++)
			for(mm = 0; mm <= ll; mm++)
			{
				p = E->sphere.hindex[ll][mm];
				E->sphere.tableplm[i][p] = modified_plgndr_a(ll, mm, t);
			}
	}


	for(j = 1; j <= E->sphere.lnoy; j++)
	{
		node = E->sphere.lexs + 1 + (E->sphere.leys + j - 1) * E->sphere.nox;
		f = E->sphere.sx[2][node];
		for(mm = 0; mm <= E->sphere.llmax; mm++)
		{
			E->sphere.tablecosf_n[j][mm] = cos((double)(mm) * f);
			E->sphere.tablesinf_n[j][mm] = sin((double)(mm) * f);
		}
	}

	for(i = 1; i <= E->sphere.lnox; i++)
	{
		node = E->sphere.lexs + i + (E->sphere.leys + 1 - 1) * E->sphere.nox;
		t = E->sphere.sx[1][node];
		/* fprintf(E->fp_out,"kkk %d %d %d %d %.5e\n",i,node,E->sphere.lexs,E->sphere.leys,t); */
		for(ll = 0; ll <= E->sphere.llmax; ll++)
			for(mm = 0; mm <= ll; mm++)
			{
				p = E->sphere.hindex[ll][mm];
				E->sphere.tableplm_n[i][p] = modified_plgndr_a(ll, mm, t);
				/* fprintf(E->fp_out,"%d %d %.5e\n",ll,mm,E->sphere.tableplm_n[i][p]); */
			}
	}

	return;
}
