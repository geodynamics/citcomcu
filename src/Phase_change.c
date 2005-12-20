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
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

void phase_change(struct All_variables *E, float *B6, float *B_b6, float *B4, float *B_b4)
{
	//int i, j, k, n, ns, m;
	int i, j, k, n, ns;
	double e_pressure, pt5, one, temp1;
	float *Xtmp[4];
	static int been_here = 0;
	static float *H;

	if(been_here == 0)
	{
		H = (float *)malloc((E->lmesh.noz + 1) * sizeof(float));

		E->control.width670 = E->monitor.length_scale / E->control.width670;
		E->control.width410 = E->monitor.length_scale / E->control.width410;

		E->control.transT670 = E->control.transT670 / E->data.ref_temperature;
		E->control.transT410 = E->control.transT410 / E->data.ref_temperature;

		E->control.Ra_670 = E->control.Ra_670 * E->control.Atemp / (E->data.density * E->data.therm_exp * E->data.ref_temperature);
		E->control.Ra_410 = E->control.Ra_410 * E->control.Atemp / (E->data.density * E->data.therm_exp * E->data.ref_temperature);

		E->control.clapeyron670 = E->control.clapeyron670 * E->data.ref_temperature / (E->data.density * E->data.grav_acc * E->monitor.length_scale);
		E->control.clapeyron410 = E->control.clapeyron410 * E->data.ref_temperature / (E->data.density * E->data.grav_acc * E->monitor.length_scale);

		fprintf(E->fp, "%g %g %g %g %g %g %g %g\n", E->control.clapeyron410, E->control.clapeyron670, E->control.Ra_410, E->control.Ra_670, E->control.transT410, E->control.transT670, E->control.width410, E->control.width670);
		fflush(E->fp);
	}


	pt5 = 0.5;
	one = 1.0;

	return_horiz_ave(E, E->T, H);

	if(E->control.Rsphere)
	{
		for(i = 1; i <= E->mesh.nsd; i++)
			Xtmp[i] = E->SX[i];
	}
	else if(E->control.CART3D)
	{
		for(i = 1; i <= E->mesh.nsd; i++)
			Xtmp[i] = E->X[i];
	}



	if(been_here == 0 || E->monitor.solution_cycles % 10 == 0)
	{
		temp1 = 0.0;
		for(i = 1; i < E->lmesh.noz; i++)
		{
			if(E->viscosity.zlm <= Xtmp[3][i + 1] && E->viscosity.zlm >= Xtmp[3][i])
			{
				temp1 = H[i] + (H[i + 1] - H[i]) * (E->viscosity.zlm - Xtmp[3][i]) / (Xtmp[3][i + 1] - Xtmp[3][i]);
				break;
			}
		}
		E->control.transT670 = temp1;
		if(E->parallel.nprocz > 1)
			E->control.transT670 = sum_across_depth(E, temp1);

		temp1 = 0.0;
		for(i = 1; i < E->lmesh.noz; i++)
		{
			if(E->viscosity.z410 <= Xtmp[3][i + 1] && E->viscosity.z410 >= Xtmp[3][i])
			{
				temp1 = H[i] + (H[i + 1] - H[i]) * (E->viscosity.z410 - Xtmp[3][i]) / (Xtmp[3][i + 1] - Xtmp[3][i]);
				break;
			}
		}
		E->control.transT410 = temp1;
		if(E->parallel.nprocz > 1)
			E->control.transT410 = sum_across_depth(E, temp1);

	}

	for(i = 1; i <= E->lmesh.nno; i++)
	{
		e_pressure = E->viscosity.zlm - Xtmp[3][i] - E->control.clapeyron670 * (E->T[i] - E->control.transT670);
		B6[i] = pt5 * (one + tanh(E->control.width670 * e_pressure));
	}

	for(i = 1; i <= E->lmesh.nno; i++)
	{
		e_pressure = E->viscosity.z410 - Xtmp[3][i] - E->control.clapeyron410 * (E->T[i] - E->control.transT410);
		B4[i] = pt5 * (one + tanh(E->control.width410 * e_pressure));
	}

	if(E->monitor.solution_cycles % (10 * E->control.record_every) == 0)
	{
		ns = 0;
		for(k = 1; k <= E->lmesh.noy; k++)
			for(j = 1; j <= E->lmesh.nox; j++)
			{
				ns = ns + 1;
				B_b6[ns] = 0.0;
				B_b4[ns] = 0.0;
				for(i = 1; i < E->lmesh.noz; i++)
				{
					n = (k - 1) * E->lmesh.noz * E->lmesh.nox + (j - 1) * E->lmesh.noz + i;
					if(B6[n] >= pt5 && B6[n + 1] <= pt5)
						B_b6[ns] = (Xtmp[3][n + 1] - Xtmp[3][n]) * (pt5 - B6[n]) / (B6[n + 1] - B6[n]) + Xtmp[3][n];
					if(B4[n] >= pt5 && B4[n + 1] <= pt5)
						B_b4[ns] = (Xtmp[3][n + 1] - Xtmp[3][n]) * (pt5 - B4[n]) / (B4[n + 1] - B4[n]) + Xtmp[3][n];
				}
			}
	}


	if(E->monitor.solution_cycles % E->control.record_every == 0)
	{
		fprintf(E->fp, "fas=%g %g %g %g %g %g %g %g %g %g\n", E->control.clapeyron410, E->control.clapeyron670, E->control.Ra_410, E->control.Ra_670, E->control.transT410, E->control.transT670, E->control.width410, E->control.width670, E->viscosity.zlm, E->viscosity.z410);
		fflush(E->fp);
	}

	been_here = 1;


	return;
}
