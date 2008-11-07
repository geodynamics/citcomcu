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

/*  Here are the routines which process the results of each velocity solution, and call
    the relevant output routines. At this point, the velocity and pressure fields have
    been calculated and stored at the nodes. The only properties of the velocity field
    which are already known are those required to check convergence of the iterative
    scheme and so on. */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include <stdlib.h>				/* for "system" command */

#include "element_definitions.h"
#include "global_defs.h"

void process_new_velocity(struct All_variables *E, int ii)
{

  if(E->control.stokes || ((ii % E->control.record_every) == 0))
    {
      /* get_CBF_topo(E,E->slice.tpg,E->slice.tpgb); */
      
      get_STD_topo(E, E->slice.tpg, E->slice.tpgb, ii);

      averages(E);

#ifdef USE_GZDIR
      if(E->control.gzdir)
	output_velo_related_gzdir(E, ii);	/* also topo */
      else
	output_velo_related(E, ii);	/* also topo */
#else
      output_velo_related(E, ii);	/* also topo */
#endif
    }

  return;
}

/* ===============================================   */

void get_surface_velo(struct All_variables *E, float *SV)
{
	//FILE *fp;
	//char output_file[255];

	//int el, els, i, m, node, lev;
	int i, m, node, lev;

	//const int dims = E->mesh.nsd;
	//const int ends = enodes[dims];
	const int nno = E->lmesh.nno;

	lev = E->mesh.levmax;

	m = 0;

	for(node = 1; node <= nno; node++)
		if((node - 1) % E->lmesh.noz == 0)
		{
			i = (node - 1) / E->lmesh.noz + 1;
			SV[(i - 1) * 2 + 1] = E->V[1][node];
			SV[(i - 1) * 2 + 2] = E->V[2][node];
		}

	return;
}

/* ===============================================   */

void get_ele_visc(struct All_variables *E, float *EV)
{
	int el, j, lev;

	const int nel = E->lmesh.nel;
	const int vpts = vpoints[E->mesh.nsd];

	lev = E->mesh.levmax;

	for(el = 1; el <= nel; el++)
	{
		EV[el] = 0.0;
		for(j = 1; j <= vpts; j++)
			EV[el] += E->EVI[lev][(el - 1) * vpts + j];

		EV[el] /= vpts;
	}

	return;
}


void get_surf_stress(struct All_variables *E, float *SXX, float *SYY, float *SZZ, float *SXY, float *SXZ, float *SZY)
{
	int i, node, stride;

	stride = E->lmesh.nsf * 6;

	for(node = 1; node <= E->lmesh.nno; node++)
		if(((node - 1) % E->lmesh.noz) == 0)
		{
			i = (node - 1) / E->lmesh.noz + 1;
			E->stress[(i - 1) * 6 + 1] = SXX[node];
			E->stress[(i - 1) * 6 + 2] = SZZ[node];
			E->stress[(i - 1) * 6 + 3] = SYY[node];
			E->stress[(i - 1) * 6 + 4] = SXY[node];
			E->stress[(i - 1) * 6 + 5] = SXZ[node];
			E->stress[(i - 1) * 6 + 6] = SZY[node];
		}
		else if(((node - 2) % E->lmesh.noz) == 0)
		{
			i = (node - 2) / E->lmesh.noz + 1;
			E->stress[stride + (i - 1) * 6 + 1] = SXX[node];
			E->stress[stride + (i - 1) * 6 + 2] = SZZ[node];
			E->stress[stride + (i - 1) * 6 + 3] = SYY[node];
			E->stress[stride + (i - 1) * 6 + 4] = SXY[node];
			E->stress[stride + (i - 1) * 6 + 5] = SXZ[node];
			E->stress[stride + (i - 1) * 6 + 6] = SZY[node];
		}

	return;
}

void averages(struct All_variables *E)
{
	//int lev, i, j, el;
	int lev, i;
	float *temp, z_thld;

	lev = E->mesh.levmax;

	temp = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));

	visc_from_gint_to_nodes(E, E->EVI[lev], temp, lev);
	return_horiz_ave(E, temp, E->Have.Vi);
	return_horiz_ave(E, E->C, E->Have.C);

	z_thld = -0.1;
//  fprintf(E->fp,"oooo\n");fflush(E->fp);
	E->monitor.Sigma_interior = return_bulk_value(E, E->C, z_thld, 0);

	z_thld = E->viscosity.zcomp + 0.1;
	E->monitor.Sigma_max = return_bulk_value(E, E->C, z_thld, 0);

	if(E->mesh.nsd == 2)
		for(i = 1; i <= E->lmesh.nno; i++)
		{
			temp[i] = E->V[1][i] * E->V[1][i] + E->V[2][i] * E->V[2][i];
		}
	else
		for(i = 1; i <= E->lmesh.nno; i++)
		{
			temp[i] = E->V[1][i] * E->V[1][i] + E->V[2][i] * E->V[2][i] + E->V[3][i] * E->V[3][i];
		}

	return_horiz_ave(E, temp, E->Have.vrms);

	for(i = 1; i <= E->lmesh.noz; i++)
		E->Have.vrms[i] = sqrt(E->Have.vrms[i]);

/*
  plume_buoyancy_flux(E);
*/


	free((void *)temp);

	return;
}
