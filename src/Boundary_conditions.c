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

#include "element_definitions.h"
#include "global_defs.h"
#include <math.h>


/* ========================================== */

void velocity_boundary_conditions(struct All_variables *E)
{
	int lv;
	//int node, d;
	int node;

	for(lv = E->mesh.levmax; lv >= E->mesh.levmin; lv--)
	{
		if(E->mesh.botvbc != 1)
		{
			horizontal_bc(E, E->VB, 1, 1, 0.0, VBX, 0, lv);
			horizontal_bc(E, E->VB, 1, 3, 0.0, VBZ, 1, lv);
			horizontal_bc(E, E->VB, 1, 2, 0.0, VBY, 0, lv);
			horizontal_bc(E, E->VB, 1, 1, E->control.VBXbotval, SBX, 1, lv);
			horizontal_bc(E, E->VB, 1, 3, 0.0, SBZ, 0, lv);
			horizontal_bc(E, E->VB, 1, 2, E->control.VBYbotval, SBY, 1, lv);
		}
		if(E->mesh.topvbc != 1)
		{
			horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 1, 0.0, VBX, 0, lv);
			horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 3, 0.0, VBZ, 1, lv);
			horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 2, 0.0, VBY, 0, lv);
			horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 1, E->control.VBXtopval, SBX, 1, lv);
			horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 3, 0.0, SBZ, 0, lv);
			horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 2, E->control.VBYtopval, SBY, 1, lv);
		}
	}

	if(E->mesh.periodic_x || E->mesh.periodic_y)
		velocity_apply_periodic_bcs(E);
	else
		velocity_refl_vert_bc(E);	/* default */

	for(lv = E->mesh.levmax; lv >= E->mesh.levmin; lv--)
	{
		if(E->mesh.botvbc == 1)
		{
			horizontal_bc(E, E->VB, 1, 1, E->control.VBXbotval, VBX, 1, lv);
			horizontal_bc(E, E->VB, 1, 3, 0.0, VBZ, 1, lv);
			horizontal_bc(E, E->VB, 1, 2, E->control.VBYbotval, VBY, 1, lv);
			horizontal_bc(E, E->VB, 1, 1, 0.0, SBX, 0, lv);
			horizontal_bc(E, E->VB, 1, 3, 0.0, SBZ, 0, lv);
			horizontal_bc(E, E->VB, 1, 2, 0.0, SBY, 0, lv);
		}
		if(E->mesh.topvbc == 1)
		{
			E->control.VBXtopval = E->control.plate_vel;
			horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 1, E->control.VBXtopval, VBX, 1, lv);
			horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 3, 0.0, VBZ, 1, lv);
			horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 2, E->control.VBYtopval, VBY, 1, lv);
			horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 1, 0.0, SBX, 0, lv);
			horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 3, 0.0, SBZ, 0, lv);
			horizontal_bc(E, E->VB, E->mesh.NOZ[lv], 2, 0.0, SBY, 0, lv);
		}
	}


	if(E->control.verbose)
	{
		for(node = 1; node <= E->lmesh.nno; node++)
			fprintf(E->fp, "VB== %d %g %g %g\n", node, E->VB[1][node], E->VB[2][node], E->VB[3][node]);
		for(lv = E->mesh.levmax; lv >= E->mesh.levmin; lv--)
		{
			fprintf(E->fp, "VBB level=%d %d\n", lv, E->lmesh.NNO[lv]);
			for(node = 1; node <= E->lmesh.NNO[lv]; node++)
			{
				fprintf(E->fp, "VB== %d %u %u %u\n", node, E->NODE[lv][node] & VBX, E->NODE[lv][node] & VBY, E->NODE[lv][node] & VBZ);
				fprintf(E->fp, "SB== %d %u %u %u\n", node, E->NODE[lv][node] & SBX, E->NODE[lv][node] & SBY, E->NODE[lv][node] & SBZ);
			}
		}
	}


	/* If any imposed internal velocity structure it goes here */


	return;
}

/* ========================================== */

void temperature_boundary_conditions(struct All_variables *E)
{
	int node;

	if(E->mesh.bottbc == 1)
	{
		horizontal_bc(E, E->TB, 1, 3, E->control.TBCbotval, TBZ, 1, E->mesh.levmax);
		horizontal_bc(E, E->TB, 1, 3, E->control.TBCbotval, FBZ, 0, E->mesh.levmax);
	}
	else
	{
		horizontal_bc(E, E->TB, 1, 3, E->control.TBCbotval, TBZ, 0, E->mesh.levmax);
		horizontal_bc(E, E->TB, 1, 3, E->control.TBCbotval, FBZ, 1, E->mesh.levmax);
	}

	if(E->mesh.toptbc == 1)
	{
		horizontal_bc(E, E->TB, E->mesh.noz, 3, E->control.TBCtopval, TBZ, 1, E->mesh.levmax);
		horizontal_bc(E, E->TB, E->mesh.noz, 3, E->control.TBCtopval, FBZ, 0, E->mesh.levmax);
	}
	else
	{
		horizontal_bc(E, E->TB, E->mesh.noz, 3, E->control.TBCtopval, TBZ, 0, E->mesh.levmax);
		horizontal_bc(E, E->TB, E->mesh.noz, 3, E->control.TBCtopval, FBZ, 1, E->mesh.levmax);
	}



	if(E->mesh.periodic_x || E->mesh.periodic_y)
		temperature_apply_periodic_bcs(E);
	else
		temperature_refl_vert_bc(E);	/* default */


	temperatures_conform_bcs(E);

	if(E->control.verbose)
	{
		for(node = 1; node <= E->lmesh.nno; node++)
			fprintf(E->fp, "TB== %d %g %g %g\n", node, E->TB[1][node], E->TB[2][node], E->TB[3][node]);
		for(node = 1; node <= E->lmesh.nno; node++)
			fprintf(E->fp, "TB== %d %u %u %u\n", node, E->node[node] & TBX, E->node[node] & TBY, E->node[node] & TBZ);
	}

	return;
}

/* ========================================== */

void velocity_refl_vert_bc(struct All_variables *E)
{
	int i, j, ii, jj;
	int node1, node2;
	int level, nox, noy, noz;
	//const int dims = E->mesh.nsd;

	/* except one side with XOZ and y=0, all others are not reflecting BC */
	/* for two YOZ planes if 3-D, or two OZ side walls for 2-D */

	if(E->parallel.me_loc[1] == 0 || E->parallel.me_loc[1] == E->parallel.nprocx - 1)
		for(j = 1; j <= E->lmesh.noy; j++)
			for(i = 1; i <= E->lmesh.noz; i++)
			{
				node1 = i + (j - 1) * E->lmesh.noz * E->lmesh.nox;
				node2 = node1 + (E->lmesh.nox - 1) * E->lmesh.noz;

				ii = i + E->lmesh.nzs - 1;
				if(E->parallel.me_loc[1] == 0)
				{
					E->VB[1][node1] = 0.0;
					if((ii != 1) && (ii != E->mesh.noz))
						E->VB[3][node1] = 0.0;
				}
				if(E->parallel.me_loc[1] == E->parallel.nprocx - 1)
				{
					E->VB[1][node2] = 0.0;
					if((ii != 1) && (ii != E->mesh.noz))
						E->VB[3][node2] = 0.0;
				}
			}					/* end loop for i and j */

	/* for two XOZ planes if 3-D */

	if(E->parallel.me_loc[2] == 0 || E->parallel.me_loc[2] == E->parallel.nprocy - 1)
		for(j = 1; j <= E->lmesh.nox; j++)
			for(i = 1; i <= E->lmesh.noz; i++)
			{
				node1 = i + (j - 1) * E->lmesh.noz;
				node2 = node1 + (E->lmesh.noy - 1) * E->lmesh.noz * E->lmesh.nox;
				ii = i + E->lmesh.nzs - 1;

				if(E->parallel.me_loc[2] == 0)
				{
					E->VB[2][node1] = 0.0;
					if((ii != 1) && (ii != E->mesh.noz))
						E->VB[3][node1] = 0.0;
				}
				if(E->parallel.me_loc[2] == E->parallel.nprocy - 1)
				{
					E->VB[2][node2] = 0.0;
					if((ii != 1) && (ii != E->mesh.noz))
						E->VB[3][node2] = 0.0;
				}
			}					/* end of loop i & j */

	/* all vbc's apply at all levels  */
	for(level = E->mesh.levmax; level >= E->mesh.levmin; level--)
	{
		nox = E->lmesh.NOX[level];
		noz = E->lmesh.NOZ[level];
		noy = E->lmesh.NOY[level];

		if(E->parallel.me_loc[1] == 0 || E->parallel.me_loc[1] == E->parallel.nprocx - 1)
			for(j = 1; j <= noy; j++)
				for(i = 1; i <= noz; i++)
				{
					node1 = i + (j - 1) * noz * nox;
					node2 = node1 + (nox - 1) * noz;
					ii = i + E->lmesh.NZS[level] - 1;
					if(E->parallel.me_loc[1] == 0)
					{
						E->NODE[level][node1] = E->NODE[level][node1] & (~SBX);
						E->NODE[level][node1] = E->NODE[level][node1] | (VBX);
						if((ii != 1) && (ii != E->mesh.NOZ[level]))
						{
							E->NODE[level][node1] = E->NODE[level][node1] & (~VBY);
							E->NODE[level][node1] = E->NODE[level][node1] | SBY;
							E->NODE[level][node1] = E->NODE[level][node1] & (~VBZ);
							E->NODE[level][node1] = E->NODE[level][node1] | SBZ;
						}
					}
					if(E->parallel.me_loc[1] == E->parallel.nprocx - 1)
					{
						E->NODE[level][node2] = E->NODE[level][node2] & (~SBX);
						E->NODE[level][node2] = E->NODE[level][node2] | (VBX);
						if((ii != 1) && (ii != E->mesh.NOZ[level]))
						{
							E->NODE[level][node2] = E->NODE[level][node2] & (~VBY);
							E->NODE[level][node2] = E->NODE[level][node2] | SBY;
							E->NODE[level][node2] = E->NODE[level][node2] & (~VBZ);
							E->NODE[level][node2] = E->NODE[level][node2] | SBZ;
						}
					}
				}				/* end for loop i & j */

		if(E->parallel.me_loc[2] == 0 || E->parallel.me_loc[2] == E->parallel.nprocy - 1)
			for(j = 1; j <= nox; j++)
				for(i = 1; i <= noz; i++)
				{
					node1 = i + (j - 1) * noz;
					node2 = node1 + (noy - 1) * noz * nox;
					ii = i + E->lmesh.NZS[level] - 1;
					jj = j + E->lmesh.NXS[level] - 1;
					if(E->parallel.me_loc[2] == 0)
					{
						E->NODE[level][node1] = E->NODE[level][node1] | VBY;
						E->NODE[level][node1] = E->NODE[level][node1] & (~SBY);
						if((ii != 1) && (ii != E->mesh.NOZ[level]))
						{
							E->NODE[level][node1] = E->NODE[level][node1] & (~VBZ);
							E->NODE[level][node1] = E->NODE[level][node1] | SBZ;
						}
						if((jj != 1) && (jj != E->mesh.NOX[level]) && (ii != 1) && (ii != E->mesh.NOZ[level]))
						{
							E->NODE[level][node1] = E->NODE[level][node1] & (~VBX);
							E->NODE[level][node1] = E->NODE[level][node1] | SBX;
						}
					}
					if(E->parallel.me_loc[2] == E->parallel.nprocy - 1)
					{
						E->NODE[level][node2] = E->NODE[level][node2] | VBY;
						E->NODE[level][node2] = E->NODE[level][node2] & (~SBY);
						if((ii != 1) && (ii != E->mesh.NOZ[level]))
						{
							E->NODE[level][node2] = E->NODE[level][node2] & (~VBZ);
							E->NODE[level][node2] = E->NODE[level][node2] | SBZ;
						}
						if((jj != 1) && (jj != E->mesh.NOX[level]) && (ii != 1) && (ii != E->mesh.NOZ[level]))
						{
							E->NODE[level][node2] = E->NODE[level][node2] & (~VBX);
							E->NODE[level][node2] = E->NODE[level][node2] | SBX;
						}
					}

				}				/* end for loop i & j  */
	}							/* end for loop level */


	return;
}

void temperature_refl_vert_bc(struct All_variables *E)
{
	int i, j;
	int node1, node2;
	//const int dims = E->mesh.nsd;

	/* Temps and bc-values  at top level only */
	/* fixed temperature at x=0 */

	if(E->parallel.me_loc[1] == 0 || E->parallel.me_loc[1] == E->parallel.nprocx - 1)
		for(j = 1; j <= E->lmesh.noy; j++)
			for(i = 1; i <= E->lmesh.noz; i++)
			{
				node1 = i + (j - 1) * E->lmesh.noz * E->lmesh.nox;
				node2 = node1 + (E->lmesh.nox - 1) * E->lmesh.noz;
				if(E->parallel.me_loc[1] == 0)
				{
					E->node[node1] = E->node[node1] & (~TBX);
					E->node[node1] = E->node[node1] | FBX;
					E->TB[1][node1] = 0.0;
				}
				if(E->parallel.me_loc[1] == E->parallel.nprocx - 1)
				{
					E->node[node2] = E->node[node2] & (~TBX);
					E->node[node2] = E->node[node2] | FBX;
					E->TB[1][node2] = 0.0;
				}
			}					/* end for loop i & j */

	if(E->parallel.me_loc[2] == 0 || E->parallel.me_loc[2] == E->parallel.nprocy - 1)
		for(j = 1; j <= E->lmesh.nox; j++)
			for(i = 1; i <= E->lmesh.noz; i++)
			{
				node1 = i + (j - 1) * E->lmesh.noz;
				node2 = node1 + (E->lmesh.noy - 1) * E->lmesh.noz * E->lmesh.nox;
				if(E->parallel.me_loc[2] == 0)
				{
					E->node[node1] = E->node[node1] & (~TBY);
					E->node[node1] = E->node[node1] | FBY;
					E->TB[2][node1] = 0.0;
				}
				if(E->parallel.me_loc[2] == E->parallel.nprocy - 1)
				{
					E->node[node2] = E->node[node2] & (~TBY);
					E->node[node2] = E->node[node2] | FBY;
					E->TB[2][node2] = 0.0;
				}
			}					/* end loop for i and j */

	return;
}


void temperature_imposed_botm_bcs(struct All_variables *E, float *BC[], int dirn)
{
	//int i, j, node, rowl;
	int i, j, node;
	const int dims = E->mesh.nsd;
	const int level = E->mesh.levmax;
	//const int noz = E->lmesh.NOZ[E->mesh.levmax];
	float dT, aa2, rr2;

	aa2 = E->segment.plume_radius * E->segment.plume_radius;


	if(E->parallel.me_loc[3] == E->parallel.nprocz - 1)
	{
		for(j = 1; j <= E->lmesh.NOY[level]; j++)
			for(i = 1; i <= E->lmesh.NOX[level]; i++)
			{
				node = 1 + (i - 1) * E->lmesh.NOZ[level] + (j - 1) * E->lmesh.NOZ[level] * E->lmesh.NOX[level];
				if(dims == 2)
					rr2 = (E->X[1][node] - E->segment.plume_coord[1]) * (E->X[1][node] - E->segment.plume_coord[1]);
				else
					rr2 = (E->X[1][node] - E->segment.plume_coord[1]) * (E->X[1][node] - E->segment.plume_coord[1]) + (E->X[2][node] - E->segment.plume_coord[2]) * (E->X[2][node] - E->segment.plume_coord[2]);
				dT = E->segment.plume_DT * exp(-rr2 / aa2);
				BC[dirn][node] += dT;
			}					/* end for loop i & j */
	}
	return;
}

/*  =========================================================  */


void horizontal_bc(struct All_variables *E, float *BC[], int ROW, int dirn, float value, unsigned int mask, char onoff, int level)

{
	int i, j, node, rowl;
	//const int dims = E->mesh.nsd;

	/* safety feature */
	if(dirn > E->mesh.nsd)
		return;

	if(ROW == 1)
		rowl = 1;
	else
		rowl = E->lmesh.NOZ[level];

	if(ROW == 1 && E->parallel.me_loc[3] == 0 || ROW == E->mesh.NOZ[level] && E->parallel.me_loc[3] == E->parallel.nprocz - 1)
	{

		/* turn bc marker to zero */
		if(onoff == 0)
			for(j = 1; j <= E->lmesh.NOY[level]; j++)
				for(i = 1; i <= E->lmesh.NOX[level]; i++)
				{
					node = rowl + (i - 1) * E->lmesh.NOZ[level] + (j - 1) * E->lmesh.NOZ[level] * E->lmesh.NOX[level];
					E->NODE[level][node] = E->NODE[level][node] & (~mask);
				}				/* end for loop i & j */

		/* turn bc marker to one */
		else
			for(j = 1; j <= E->lmesh.NOY[level]; j++)
				for(i = 1; i <= E->lmesh.NOX[level]; i++)
				{
					node = rowl + (i - 1) * E->lmesh.NOZ[level] + (j - 1) * E->lmesh.NOZ[level] * E->lmesh.NOX[level];
					E->NODE[level][node] = E->NODE[level][node] | (mask);
					if(level == E->mesh.levmax)	/* NB */
						BC[dirn][node] = value;
				}				/* end for loop i & j */

	}							/* end for if */

	return;
}


void velocity_apply_periodic_bcs(struct All_variables *E)
{
	//int n1, n2, level;
	//int i, j, ii, jj;
	//const int dims = E->mesh.nsd;

	fprintf(E->fp, "Periodic boundary conditions\n");

	return;
}

void temperature_apply_periodic_bcs(struct All_variables *E)
{
	//int n1, n2, e1, level;
	//int i, j, ii, jj;
	//const int dims = E->mesh.nsd;

	fprintf(E->fp, "Periodic temperature boundary conditions\n");

	return;
}



void strip_bcs_from_residual(struct All_variables *E, double *Res, int level)
{
	int i;

	//const int dims = E->mesh.nsd;
	//const int dofs = E->mesh.dof;
	const int nno = E->lmesh.NNO[level];
	//const int addi_dof = additional_dof[dims];

	for(i = 1; i <= nno; i++)
	{
		if((E->NODE[level][i] & VBX) != 0)
			Res[E->ID[level][i].doff[1]] = 0.0;
		if((E->NODE[level][i] & VBY) != 0)
			Res[E->ID[level][i].doff[2]] = 0.0;
		if((E->NODE[level][i] & VBZ) != 0)
			Res[E->ID[level][i].doff[3]] = 0.0;

	}

	return;
}

void temperatures_conform_bcs(struct All_variables *E)
{
	int node;
	unsigned int type;

	for(node = 1; node <= E->lmesh.nno; node++)
	{
		type = (E->node[node] & (TBX | TBZ | TBY));

		switch (type)
		{
		case 0:				/* no match, next node */
			break;
		case TBX:
			E->T[node] = E->TB[1][node];
			break;
		case TBZ:
			E->T[node] = E->TB[3][node];
			break;
		case TBY:
			E->T[node] = E->TB[2][node];
			break;
		case (TBX | TBZ):		/* clashes ! */
			E->T[node] = 0.5 * (E->TB[1][node] + E->TB[3][node]);
			break;
		case (TBX | TBY):		/* clashes ! */
			E->T[node] = 0.5 * (E->TB[1][node] + E->TB[2][node]);
			break;
		case (TBZ | TBY):		/* clashes ! */
			E->T[node] = 0.5 * (E->TB[2][node] + E->TB[3][node]);
			break;
		case (TBZ | TBY | TBX):	/* clashes ! */
			E->T[node] = 0.3333333 * (E->TB[1][node] + E->TB[2][node] + E->TB[3][node]);
			break;
		}

		/* next node */
	}

	return;

}


void velocities_conform_bcs(struct All_variables *E, double *U)
{
	//int node, d;
	int node;

	const unsigned int typex = VBX;
	const unsigned int typez = VBZ;
	const unsigned int typey = VBY;

	//const int dofs = E->mesh.dof;
	const int nno = E->lmesh.nno;

	for(node = 1; node <= nno; node++)
	{
		if(E->node[node] & typex)
			U[E->id[node].doff[1]] = E->VB[1][node];
		if(E->node[node] & typez)
			U[E->id[node].doff[3]] = E->VB[3][node];
		if(E->node[node] & typey)
			U[E->id[node].doff[2]] = E->VB[2][node];

	}

	return;
}


void equalize_id_ien_lm(struct All_variables *E)
{

	if(E->mesh.periodic_x && E->parallel.nprocx == 1)
	{

	}
	if(E->mesh.periodic_y && E->parallel.nprocy == 1)
	{

	}

	return;
}
