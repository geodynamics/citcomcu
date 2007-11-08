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

/* Functions relating to the determination of viscosity field either
   as a function of the run, as an initial condition or as specified from
   a previous file */

#include <mpi.h>
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"

void viscosity_parameters(struct All_variables *E)
{
	int i, l;
	float temp;

	/* default values .... */
	E->viscosity.update_allowed = 0;
	E->viscosity.SDEPV = E->viscosity.TDEPV = 0;
	E->viscosity.EXPX = 0;


	for(i = 0; i < 40; i++)
	{
		E->viscosity.N0[i] = 1.0;
		E->viscosity.T[i] = 0.0;
		E->viscosity.Z[i] = 0.0;
		E->viscosity.E[i] = 0.0;
		E->viscosity.T0[i] = 0.0;
	}

	/* read in information */
	input_int("rheol", &(E->viscosity.RHEOL), "essential");
	input_int("num_mat", &(E->viscosity.num_mat), "1");

	input_float_vector("viscT", E->viscosity.num_mat, (E->viscosity.T));	/* redundant */
	input_float_vector("viscT1", E->viscosity.num_mat, (E->viscosity.T));
	input_float_vector("viscZ", E->viscosity.num_mat, (E->viscosity.Z));
	input_float_vector("viscE", E->viscosity.num_mat, (E->viscosity.E));
	input_float_vector("viscT0", E->viscosity.num_mat, (E->viscosity.T0));
	input_float_vector("visc0", E->viscosity.num_mat, (E->viscosity.N0));	/* redundant */
	input_float_vector("viscN0", E->viscosity.num_mat, (E->viscosity.N0));

	input_boolean("TDEPV", &(E->viscosity.TDEPV), "on");
	input_boolean("SDEPV", &(E->viscosity.SDEPV), "off");

	input_float("sdepv_misfit", &(E->viscosity.sdepv_misfit), "0.001");
	input_float_vector("sdepv_expt", E->viscosity.num_mat, (E->viscosity.sdepv_expt));
	input_float_vector("sdepv_trns", E->viscosity.num_mat, (E->viscosity.sdepv_trns));

	input_boolean("TDEPV_AVE", &(E->viscosity.TDEPV_AVE), "off");
	input_boolean("VFREEZE", &(E->viscosity.FREEZE), "off");
	input_boolean("VMAX", &(E->viscosity.MAX), "off");
	input_boolean("VMIN", &(E->viscosity.MIN), "off");
	input_boolean("VISC_UPDATE", &(E->viscosity.update_allowed), "on");

	input_float("freeze_thresh", &(E->viscosity.freeze_thresh), "0.0");
	input_float("freeze_value", &(E->viscosity.freeze_value), "1.0");
	input_float("visc_max", &(E->viscosity.max_value), "nodefault");
	input_float("visc_min", &(E->viscosity.min_value), "nodefault");

	input_boolean("VISC_GUESS", &(E->viscosity.guess), "off");
	input_string("visc_old_file", E->viscosity.old_file, " ");

	return;
}

void get_viscosity_option(struct All_variables *E)
{
	/* general, essential default */

	input_string("Viscosity", E->viscosity.STRUCTURE, NULL);	/* Which form of viscosity */

	input_boolean("VISC_EQUIVDD", &(E->viscosity.EQUIVDD), "off");	/* Whether to average it */
	input_int("equivdd_opt", &(E->viscosity.equivddopt), "1");
	input_int("equivdd_x", &(E->viscosity.proflocx), "1");
	input_int("equivdd_y", &(E->viscosity.proflocy), "1");

	input_int("update_every_steps", &(E->control.KERNEL), "1");

	input_boolean("VISC_SMOOTH", &(E->viscosity.SMOOTH), "off");
	input_int("visc_smooth_cycles", &(E->viscosity.smooth_cycles), "0");

	if(strcmp(E->viscosity.STRUCTURE, "system") == 0)	/* Interpret */
	{
		if(E->parallel.me == 0)
			fprintf(E->fp, "Viscosity derived from system state\n");
		E->viscosity.FROM_SYSTEM = 1;
		viscosity_for_system(E);
	}

	return;

}

/* ============================================ */


void viscosity_for_system(struct All_variables *E)
{

	get_system_viscosity(E, 1, E->EVI[E->mesh.levmax], E->VI[E->mesh.levmax]);

	return;
}


void get_system_viscosity(struct All_variables *E, int propogate, float *evisc, float *visc)
{
	int i, j;
	//float *visc_old, *evisc_old;

	const int vpts = vpoints[E->mesh.nsd];


	if(E->viscosity.TDEPV)
		visc_from_T(E, visc, evisc, propogate);
	else
		visc_from_mat(E, visc, evisc);

	if(E->viscosity.SDEPV)
		visc_from_S(E, visc, evisc, propogate);

	if(E->viscosity.SMOOTH)
		apply_viscosity_smoother(E, visc, evisc);

	if(E->viscosity.MAX)
	{
		for(i = 1; i <= E->lmesh.nel; i++)
			for(j = 1; j <= vpts; j++)
			{
				if(evisc[(i - 1) * vpts + j] > E->viscosity.max_value)
					evisc[(i - 1) * vpts + j] = E->viscosity.max_value;
			}
	}

	if(E->viscosity.MIN)
	{
		for(i = 1; i <= E->lmesh.nel; i++)
			for(j = 1; j <= vpts; j++)
				if(evisc[(i - 1) * vpts + j] < E->viscosity.min_value)
					evisc[(i - 1) * vpts + j] = E->viscosity.min_value;
	}

	/*   v_to_nodes(E,evisc,visc,E->mesh.levmax);  */

	return;
}


void apply_viscosity_smoother(struct All_variables *E, float *visc, float *evisc)

{
	double *ViscCentre;
	int i;

	ViscCentre = (double *)malloc((E->lmesh.nno + 10) * sizeof(double));

	for(i = 1; i <= E->viscosity.smooth_cycles; i++)
	{
		p_to_centres(E, visc, ViscCentre, E->mesh.levmax);
		p_to_nodes(E, ViscCentre, visc, E->mesh.levmax);
	}

	free((void *)ViscCentre);

	return;
}

void visc_from_mat(struct All_variables *E, float *Eta, float *EEta)
{
	//int i, j, k, l, z, jj, kk;
	int i, jj;

	for(i = 1; i <= E->lmesh.nel; i++)
		for(jj = 1; jj <= vpoints[E->mesh.nsd]; jj++)
		{
			EEta[(i - 1) * vpoints[E->mesh.nsd] + jj] = E->viscosity.N0[E->mat[i] - 1];
		}

	return;
}

void visc_from_T(struct All_variables *E, float *Eta, float *EEta, int propogate)
{
	//int i, j, k, l, z, e, jj, kk, imark;
	int i, l, e, jj, kk;
	//float c1, c2, c3, zero, e_6, one, eta0, Tave, depth, temp, tempa, TT[9];
	float zero, one, temp, tempa, TT[9];
	//double ztop, zbotm, zz, visc1, area1, temp1, temp2, temp3, temp4;
	double ztop, zbotm, zz, temp1, temp2, temp3, temp4;
	float *Xtmp[4];
	static int visits = 0;
	static float *Tadi;
	static double slope = 0;
	const int vpts = vpoints[E->mesh.nsd];
	const int ends = enodes[E->mesh.nsd];
	const int nel = E->lmesh.nel;
	const int noz = E->lmesh.noz;

	one = 1.0;
	zero = 0.0;
	temp3 = temp4 = 0.0;

	if(visits == 0)
	{
		fprintf(E->fp, "\tRheological option : %d\n", E->viscosity.RHEOL);

		for(l = 1; l <= E->viscosity.num_mat; l++)
		{
			fprintf(E->fp, "\tlayer %d/%d: E=%g T1=%g N0=%g Z0=%g\n", l, E->viscosity.num_mat, E->viscosity.E[l - 1], E->viscosity.T[l - 1], E->viscosity.N0[l - 1], E->viscosity.Z[l - 1]);
		}
		fflush(E->fp);

	}

	if(E->control.Rsphere)
	{
		for(i = 1; i <= E->mesh.nsd; i++)
			Xtmp[i] = E->SX[i];
		ztop = E->sphere.ro;
		zbotm = E->sphere.ri;
	}
	else if(E->control.CART3D)
	{
		for(i = 1; i <= E->mesh.nsd; i++)
			Xtmp[i] = E->X[i];
		ztop = 1.0;
		zbotm = 0.0;
	}


	if(visits == 0 || E->monitor.solution_cycles % E->control.KERNEL == 0)
	{

		if(E->viscosity.RHEOL == 0)
		{
			for(i = 1; i <= nel; i++)
			{
				l = E->mat[i];
				e = (i - 1) % E->lmesh.elz + 1;

				tempa = E->viscosity.N0[l - 1];

				for(kk = 1; kk <= ends; kk++)
					TT[kk] = E->T[E->ien[i].node[kk]];

				for(jj = 1; jj <= vpts; jj++)
				{
					temp = 1.0e-32;
					for(kk = 1; kk <= ends; kk++)
					{
						temp += max(zero, TT[kk]) * E->N.vpt[GNVINDEX(kk, jj)];;
					}
					EEta[(i - 1) * vpts + jj] = tempa * exp(E->viscosity.E[l - 1] * (one - temp ));
				}
			}
		}
		else if(E->viscosity.RHEOL == 1)
		{
			for(i = 1; i <= nel; i++)
			{
				l = E->mat[i];
				tempa = E->viscosity.N0[l - 1];

				for(kk = 1; kk <= ends; kk++)
					TT[kk] = E->T[E->ien[i].node[kk]];

				for(jj = 1; jj <= vpts; jj++)
				{
					temp = 1.0e-32;
					for(kk = 1; kk <= ends; kk++)
					{
						temp += max(zero, TT[kk]) * E->N.vpt[GNVINDEX(kk, jj)];;
					}
					EEta[(i - 1) * vpts + jj] = tempa * exp((E->viscosity.E[l - 1]) / (temp + E->viscosity.T[l - 1]));
				}
			}
		}
		else if(E->viscosity.RHEOL == 2)
		{
			for(i = 1; i <= nel; i++)
			{
				l = E->mat[i];
				tempa = E->viscosity.N0[l - 1];

				for(kk = 1; kk <= ends; kk++)
					TT[kk] = E->T[E->ien[i].node[kk]];

				for(jj = 1; jj <= vpts; jj++)
				{
					temp = 1.0e-32;
					zz = 0;
					for(kk = 1; kk <= ends; kk++)
					{
						temp += max(zero, TT[kk]) * E->N.vpt[GNVINDEX(kk, jj)];;
						zz += Xtmp[3][E->ien[i].node[kk]] * E->N.vpt[GNVINDEX(kk, jj)];
					}
					EEta[(i - 1) * vpts + jj] = tempa * exp((E->viscosity.E[l - 1] + (1 - zz) * E->viscosity.Z[l - 1]) /(temp + E->viscosity.T[l - 1]) );
				}
			}
		}
		else if(E->viscosity.RHEOL == 3)
		{
		}

		visits++;

	}

for(i = 1; i <= E->lmesh.nel; i++) {

    /* get temperature field on vpts */
        for(jj = 1; jj <= vpts; jj++) {
                TT[jj] = 0.0;
                        for(kk = 1; kk <= ends; kk++)
                                     TT[jj] += E->T[E->ien[i].node[kk]] * E->N.vpt[GNVINDEX(kk,jj)];
                                         }

                                             /* write to the log files */
                                                 for(jj = 1; jj <= vpts; jj++) {
                                                         fprintf(E->fp, "%d %d %e %e\n", i, jj,
                                                                         E->EVI[E->mesh.levmax][(i - 1) * vpts + jj],
                                                                                         TT[jj]);
                                                                                             }
                                                                                              
                                                                                                
                                                                                                }




/*
	fprintf(E->fp,"aaa\n");
	for(i=1;i<=nel;i++)
		fprintf(E->fp,"%d %d %g\n",i,E->mat[i],EEta[(i-1)*vpts+1]);
*/

	return;
}

void visc_from_S(struct All_variables *E, float *Eta, float *EEta, int propogate)
{
	static int visits = 0;
	//float one, two, scale, stress_magnitude, depth, exponent1;
	float one, two, scale, exponent1;
	float *eedot;

	//int e, l, z, jj, kk;
	int e, jj;

	const int vpts = vpoints[E->mesh.nsd];
	const int nel = E->lmesh.nel;

	eedot = (float *)malloc((2 + nel) * sizeof(float));
	one = 1.0;
	two = 2.0;

	if(visits == 0)
	{
		for(e = 1; e <= nel; e++)
			eedot[e] = one;
	}
	else
		strain_rate_2_inv(E, eedot, 1);

	for(e = 1; e <= nel; e++)
	{
		exponent1 = one - one / E->viscosity.sdepv_expt[E->mat[e] - 1];
		scale = pow(two * eedot[e] / E->viscosity.sdepv_trns[E->mat[e] - 1], exponent1);
		for(jj = 1; jj <= vpts; jj++)
			EEta[(e - 1) * vpts + jj] = two * EEta[(e - 1) * vpts + jj] / (one + scale * pow(EEta[(e - 1) * vpts + jj], exponent1));
	}


	visits++;

	free((void *)eedot);
	return;
}


void strain_rate_2_inv(struct All_variables *E, float *EEDOT, int SQRT)
{
	double edot[4][4], dudx[4][4], rtf[4][9];
	float VV[4][9], Vxyz[9][9];

	//int e, i, j, p, q, n, nel, k;
	int e, i, j, p, q, n, nel;

	const int dims = E->mesh.nsd;
	const int ends = enodes[dims];
	const int lev = E->mesh.levmax;
	//const int nno = E->lmesh.nno;
	//const int vpts = vpoints[dims];
	const int ppts = ppoints[dims];
	//const int sphere_key = 1;

	nel = E->lmesh.nel;

	if(E->control.Rsphere)
	{

		for(e = 1; e <= nel; e++)
		{

			get_rtf(E, e, 2, rtf, lev);

			for(i = 1; i <= ends; i++)
			{
				n = E->ien[e].node[i];
				VV[1][i] = E->V[1][n];
				VV[2][i] = E->V[2][n];
				VV[3][i] = E->V[3][n];
			}

			for(j = 1; j <= ppts; j++)
			{
				Vxyz[1][j] = 0.0;
				Vxyz[2][j] = 0.0;
				Vxyz[3][j] = 0.0;
				Vxyz[4][j] = 0.0;
				Vxyz[5][j] = 0.0;
				Vxyz[6][j] = 0.0;
			}

			for(j = 1; j <= ppts; j++)
			{
				for(i = 1; i <= ends; i++)
				{
					Vxyz[1][j] += (VV[1][i] * E->gNX[e].ppt[GNPXINDEX(0, i, j)] + VV[3][i] * E->N.ppt[GNPINDEX(i, j)]) * rtf[3][j];
					Vxyz[2][j] += ((VV[2][i] * E->gNX[e].ppt[GNPXINDEX(1, i, j)] + VV[1][i] * E->N.ppt[GNPINDEX(i, j)] * cos(rtf[1][j])) / sin(rtf[1][j]) + VV[3][i] * E->N.ppt[GNPINDEX(i, j)]) * rtf[3][j];
					Vxyz[3][j] += VV[3][i] * E->gNX[e].ppt[GNPXINDEX(2, i, j)];

					Vxyz[4][j] += ((VV[1][i] * E->gNX[e].ppt[GNPXINDEX(1, i, j)] - VV[2][i] * E->N.ppt[GNPINDEX(i, j)] * cos(rtf[1][j])) / sin(rtf[1][j]) + VV[2][i] * E->gNX[e].ppt[GNPXINDEX(0, i, j)]) * rtf[3][j];
					Vxyz[5][j] += VV[1][i] * E->gNX[e].ppt[GNPXINDEX(2, i, j)] + rtf[3][j] * (VV[3][i] * E->gNX[e].ppt[GNPXINDEX(0, i, j)] - VV[1][i] * E->N.ppt[GNPINDEX(i, j)]);
					Vxyz[6][j] += VV[2][i] * E->gNX[e].ppt[GNPXINDEX(2, i, j)] + rtf[3][j] * (VV[3][i] * E->gNX[e].ppt[GNPXINDEX(1, i, j)] / sin(rtf[1][j]) - VV[2][i] * E->N.ppt[GNPINDEX(i, j)]);
				}
				edot[1][1] = 2.0 * Vxyz[1][j];
				edot[2][2] = 2.0 * Vxyz[2][j];
				edot[3][3] = 2.0 * Vxyz[3][j];
				edot[1][2] = Vxyz[4][j];
				edot[1][3] = Vxyz[5][j];
				edot[2][3] = Vxyz[6][j];
			}

			EEDOT[e] = edot[1][1] * edot[1][1] + edot[1][2] * edot[1][2] * 2.0 + edot[2][2] * edot[2][2] + edot[2][3] * edot[2][3] * 2.0 + edot[3][3] * edot[3][3] + edot[1][3] * edot[1][3] * 2.0;

		}

	}

	else if(E->control.CART3D)
	{

		for(e = 1; e <= nel; e++)
		{

			for(i = 1; i <= ends; i++)
			{
				n = E->ien[e].node[i];
				VV[1][i] = E->V[1][n];
				VV[2][i] = E->V[2][n];
				if(dims == 3)
					VV[3][i] = E->V[3][n];
			}

			for(p = 1; p <= dims; p++)
				for(q = 1; q <= dims; q++)
					dudx[p][q] = 0.0;

			for(i = 1; i <= ends; i++)
				for(p = 1; p <= dims; p++)
					for(q = 1; q <= dims; q++)
						dudx[p][q] += VV[p][i] * E->gNX[e].ppt[GNPXINDEX(q - 1, i, 1)];

			for(p = 1; p <= dims; p++)
				for(q = 1; q <= dims; q++)
					edot[p][q] = dudx[p][q] + dudx[q][p];

			if(dims == 2)
				EEDOT[e] = edot[1][1] * edot[1][1] + edot[2][2] * edot[2][2] + edot[1][2] * edot[1][2] * 2.0;

			else if(dims == 3)
				EEDOT[e] = edot[1][1] * edot[1][1] + edot[1][2] * edot[1][2] * 2.0 + edot[2][2] * edot[2][2] + edot[2][3] * edot[2][3] * 2.0 + edot[3][3] * edot[3][3] + edot[1][3] * edot[1][3] * 2.0;

		}
	}


	if(SQRT)
		for(e = 1; e <= nel; e++)
			EEDOT[e] = sqrt(0.5 * EEDOT[e]);
	else
		for(e = 1; e <= nel; e++)
			EEDOT[e] *= 0.5;

	return;
}




int layers(struct All_variables *E, float x3)
{
	int llayers = 0;

	if(x3 >= E->viscosity.zlith)
		llayers = 1;
	else if(x3 < E->viscosity.zlith && x3 >= E->viscosity.zlm)
		llayers = 2;
	else if(x3 < E->viscosity.zlm)
		llayers = 3;

	return (llayers);
}



int weak_zones(struct All_variables *E, int node, float t_b) /* very tentative */
{
	int wweak_zones;

	wweak_zones = 0;
	return wweak_zones;
}



float boundary_thickness(struct All_variables *E, float *H)
{
	float thickness;
	int i, j;

	thickness = 0.0;

	for(i = 1; i <= E->lmesh.noz; i++)
		H[i] = -H[i] / E->control.Atemp;

	if(E->parallel.me_loc[1] == 0)
	{
		for(i = 1; i <= E->lmesh.noz; i++)
		{
			if(H[i] > 0.45)
			{
				j = i;
				break;
			}
		}
		thickness = E->X[3][j];
	}


	MPI_Bcast(&thickness, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

	return (thickness);
}
