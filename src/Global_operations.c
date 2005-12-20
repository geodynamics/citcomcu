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
/*  Most of these functions compute average quantities or summations
 *  for parallel operations. They were implemented by SZ from 1997 and on
 *  */


#include <mpi.h>

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

/* ===============================================
   strips horizontal average from nodal field X. 
   Assumes orthogonal mesh, otherwise, horizontals
   aren't & another method is required.
   =============================================== */

void remove_horiz_ave(struct All_variables *E, float *X, float *H, int store_or_not)
{
	//int i, j, k, n, ln, nox, noz, noy;
	int i, j, k, n, nox, noz, noy;

	nox = E->lmesh.nox;
	noy = E->lmesh.noy;
	noz = E->lmesh.noz;

	return_horiz_ave(E, X, H);

/*
	for(i = 1; i <= noz; i++)
		fprintf(E->fp, "%d %.5e ave\n", i, H[i]);
*/

	for(i = 1; i <= noz; i++)
		for(k = 1; k <= noy; k++)
			for(j = 1; j <= nox; j++)
			{
				n = i + (j - 1) * noz + (k - 1) * nox * noz;
				X[n] -= H[i];
			}

	return;
}

void return_horiz_sum(struct All_variables *E, float *X, float *H, int nn)
{
	//const int dims = E->mesh.nsd;
	//int i, j, k, d, nint, noz, nox, noy, el, elz, elx, ely, j1, j2, i1, i2, k1, k2, nproc;
	int i, j, d, nproc;
	//int lnode[5], sizeofH, noz2, iroot;
	//float *Have, *temp;

	int *processors;

	MPI_Comm world, horizon_p;
	MPI_Group world_g, horizon_g;

	processors = (int *)malloc((E->parallel.nprocxy + 2) * sizeof(int));

	/* determine which processors should get the message from me for 
	 * computing the layer averages */

	nproc = 0;
	for(j = 0; j < E->parallel.nprocy; j++)
		for(i = 0; i < E->parallel.nprocx; i++)
		{
			d = E->parallel.me_loc[3] + i * E->parallel.nprocz + j * E->parallel.nprocxz;
			processors[nproc] = d;
			nproc++;
		}

	if(nproc > 1)
	{
		world = MPI_COMM_WORLD;
		MPI_Comm_group(world, &world_g);
		MPI_Group_incl(world_g, nproc, processors, &horizon_g);
		MPI_Comm_create(world, horizon_g, &horizon_p);

		MPI_Allreduce(X, H, nn, MPI_FLOAT, MPI_SUM, horizon_p);

		MPI_Comm_free(&horizon_p);
		MPI_Group_free(&horizon_g);
		MPI_Group_free(&world_g);

	}
	else
		for(i = 0; i < nn; i++)
			H[i] = X[i];

	free((void *)processors);

	return;
}


void return_horiz_ave(struct All_variables *E, float *X, float *H)
{
	//const int dims = E->mesh.nsd;
	//int i, j, k, d, nint, noz, nox, noy, el, elz, elx, ely, j1, j2, i1, i2, k1, k2, nproc;
	int i, j, k, d, nint, noz, nox, noy, el, elz, elx, ely, nproc;
	//int lnode[5], sizeofH, noz2, iroot;
	int lnode[5], noz2;
//	float *Have,*temp;
	double *Have, *temp;
	struct Shape_function1 M;
	struct Shape_function1_dA dGamma;

	int *processors;

	MPI_Comm world, horizon_p;
	MPI_Group world_g, horizon_g;


	processors = (int *)malloc((E->parallel.nprocxy + 2) * sizeof(int));
//  Have = (float *)malloc((2*E->lmesh.noz+2)*sizeof(float));
//  temp = (float *)malloc((2*E->lmesh.noz+2)*sizeof(float));
	Have = (double *)malloc((2 * E->lmesh.noz + 2) * sizeof(double));
	temp = (double *)malloc((2 * E->lmesh.noz + 2) * sizeof(double));

	noz = E->lmesh.noz;
	noy = E->lmesh.noy;
	nox = E->lmesh.nox;
	elz = E->lmesh.elz;
	elx = E->lmesh.elx;
	ely = E->lmesh.ely;
	noz2 = 2 * noz;

	for(i = 1; i <= elz; i++)
	{
		temp[i] = temp[i + noz] = 0.0;
		temp[i + 1] = temp[i + 1 + noz] = 0.0;
		for(j = 1; j <= elx; j++)
			for(k = 1; k <= ely; k++)
			{
				el = i + (j - 1) * elz + (k - 1) * elx * elz;
				get_global_1d_shape_fn(E, el, &M, &dGamma, 0);

				lnode[1] = i + (j - 1) * noz + (k - 1) * nox * noz;
				lnode[2] = i + j * noz + (k - 1) * nox * noz;
				lnode[3] = i + j * noz + k * nox * noz;
				lnode[4] = i + (j - 1) * noz + k * nox * noz;

				for(d = 1; d <= onedvpoints[E->mesh.nsd]; d++)
					for(nint = 1; nint <= onedvpoints[E->mesh.nsd]; nint++)
					{
						temp[i] += X[lnode[d]] * E->M.vpt[GMVINDEX(d, nint)] * dGamma.vpt[GMVGAMMA(0, nint)];
						temp[i + noz] += E->M.vpt[GMVINDEX(d, nint)] * dGamma.vpt[GMVGAMMA(0, nint)];
					}

				if(i == elz)
				{
					lnode[1] = 1 + i + (j - 1) * noz + (k - 1) * nox * noz;
					lnode[2] = 1 + i + j * noz + (k - 1) * nox * noz;
					lnode[3] = 1 + i + j * noz + k * nox * noz;
					lnode[4] = 1 + i + (j - 1) * noz + k * nox * noz;

					for(d = 1; d <= onedvpoints[E->mesh.nsd]; d++)
						for(nint = 1; nint <= onedvpoints[E->mesh.nsd]; nint++)
						{
							temp[i + 1] += X[lnode[d]] * E->M.vpt[GMVINDEX(d, nint)] * dGamma.vpt[GMVGAMMA(0, nint)];
							temp[i + 1 + noz] += E->M.vpt[GMVINDEX(d, nint)] * dGamma.vpt[GMVGAMMA(0, nint)];
						}
				}

			}					/* Done one traverse */

	}							/* Done for i */


	/* determine which processors should get the message from me for 
	 * computing the layer averages */

	nproc = 0;
	for(j = 0; j < E->parallel.nprocy; j++)
		for(i = 0; i < E->parallel.nprocx; i++)
		{
			d = E->parallel.me_loc[3] + i * E->parallel.nprocz + j * E->parallel.nprocxz;
			processors[nproc] = d;
			nproc++;
		}

	if(nproc > 1)
	{
		world = MPI_COMM_WORLD;
		MPI_Comm_group(world, &world_g);
		MPI_Group_incl(world_g, nproc, processors, &horizon_g);
		MPI_Comm_create(world, horizon_g, &horizon_p);

//    MPI_Allreduce(temp,Have,noz2+1,MPI_FLOAT,MPI_SUM,horizon_p);
		MPI_Allreduce(temp, Have, noz2 + 1, MPI_DOUBLE, MPI_SUM, horizon_p);

		MPI_Comm_free(&horizon_p);
		MPI_Group_free(&horizon_g);
		MPI_Group_free(&world_g);

	}
	else
		for(i = 1; i <= noz2; i++)
		{
			Have[i] = temp[i];
		}

	for(i = 1; i <= noz; i++)
	{
		if(Have[i + noz] != 0.0)
			H[i] = Have[i] / Have[i + noz];
	}

	free((void *)Have);
	free((void *)temp);
	free((void *)processors);

	return;
}


float return_bulk_value(struct All_variables *E, float *Z, float z_thld, int average)
{
	//int i, j, k, n, el, elx, ely, elz, i1, i2, j1, j2, k1, k2;
	int i, j, n, el, elx, ely, elz, i1, j1, k1;
	float integral;
//	float volume,volume1,integral1,integral0;
	double volume, volume1, integral1, integral0;

	//struct Shape_function GN;
	//struct Shape_function_dx GNx;
	//struct Shape_function_dA dOmega;

	const int vpts = vpoints[E->mesh.nsd];
	const int ends = enodes[E->mesh.nsd];

	volume1 = 0.0;
	integral1 = 0.0;


	elz = E->lmesh.elz;
	elx = E->lmesh.elx;
	ely = E->lmesh.ely;


	for(i1 = 1; i1 <= elz; i1++)
		if(E->XP[3][i1] >= z_thld)
		{
			for(j1 = 1; j1 <= elx; j1++)
				for(k1 = 1; k1 <= ely; k1++)
				{
					el = i1 + (j1 - 1) * elz + (k1 - 1) * elz * elx;

					for(j = 1; j <= vpts; j++)
						for(i = 1; i <= ends; i++)
						{
							n = E->ien[el].node[i];
							volume1 += E->N.vpt[GNVINDEX(i, j)] * E->gDA[el].vpt[j];
							integral1 += Z[n] * E->N.vpt[GNVINDEX(i, j)] * E->gDA[el].vpt[j];
						}
				}

		}


	/*
	 * MPI_Allreduce(&volume1  ,&volume  ,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
	 * MPI_Allreduce(&integral1,&integral0,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
	 */

	MPI_Allreduce(&volume1, &volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&integral1, &integral0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	integral = integral0;

	if(average && volume != 0.0)
		integral = integral0 / volume;


	return ((float)integral);
}



double global_vdot(struct All_variables *E, double *A, double *B, int lev)
{
	int i, neq;
	double prod, temp;

	neq = E->lmesh.NEQ[lev];

	temp = 0.0;
	for(i = 0; i < neq; i++)
	{
		if(E->parallel.IDD[lev][i])	/* only get the sum from relevant ID */
			temp += A[i] * B[i];
	}

	MPI_Allreduce(&temp, &prod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return (prod);
}


double global_pdot(struct All_variables *E, double *A, double *B, int lev)

{
	int i, npno;
	double prod, temp;

	npno = E->lmesh.NPNO[lev];

	temp = 0.0;
	for(i = 1; i <= npno; i++)
		temp += A[i] * B[i];

	MPI_Allreduce(&temp, &prod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return (prod);
}


float global_tdot(struct All_variables *E, float *A, float *B, int lev)

{
	int i, nno;
	float prod, temp;

	nno = E->lmesh.NNO[lev];

	temp = 0.0;
	for(i = 1; i <= nno; i++)
		temp += A[i] * B[i];

	MPI_Allreduce(&temp, &prod, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	return (prod);
}

float Tmax(struct All_variables *E, float *T)
{
	float temp, temp1;
	//int i, m;
	int i;

	temp = -10.0;
	for(i = 1; i <= E->lmesh.nno; i++)
		temp = max(T[i], temp);

	temp1 = global_fmax(E, temp);
	return (temp1);
}


float global_fmin(struct All_variables *E, float a)
{
	float temp;
	MPI_Allreduce(&a, &temp, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
	return (temp);
}

float global_fmax(struct All_variables *E, float a)
{
	float temp;
	MPI_Allreduce(&a, &temp, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
	return (temp);
}

/* ================================================== */
void sum_across_depth_sph1(struct All_variables *E, float *sphc, float *sphs)
{
	int jumpp, total, j, d;

	static float *sphcs, *temp;
	static int been_here = 0;
	static int *processors, nproc;

	static MPI_Comm world, horizon_p;
	static MPI_Group world_g, horizon_g;

	if(been_here == 0)
	{
		processors = (int *)malloc((E->parallel.nprocz + 2) * sizeof(int));
		temp = (float *)malloc((E->sphere.hindice * 2 + 3) * sizeof(float));
		sphcs = (float *)malloc((E->sphere.hindice * 2 + 3) * sizeof(float));

		nproc = 0;
		for(j = 0; j < E->parallel.nprocz; j++)
		{
			d = E->parallel.me_sph * E->parallel.nprocz + E->parallel.nprocz - 1 - j;
			processors[nproc] = d;
			nproc++;
		}

		if(nproc > 0)
		{
			world = MPI_COMM_WORLD;
			MPI_Comm_group(world, &world_g);
			MPI_Group_incl(world_g, nproc, processors, &horizon_g);
			MPI_Comm_create(world, horizon_g, &horizon_p);
		}

		been_here++;
	}

	total = E->sphere.hindice * 2 + 3;
	jumpp = E->sphere.hindice;
	for(j = 0; j < E->sphere.hindice; j++)
	{
		sphcs[j] = sphc[j];
		sphcs[j + jumpp] = sphs[j];
	}

	if(nproc > 0)
	{

		MPI_Allreduce(sphcs, temp, total, MPI_FLOAT, MPI_SUM, horizon_p);

		for(j = 0; j < E->sphere.hindice; j++)
		{
			sphc[j] = temp[j];
			sphs[j] = temp[j + jumpp];
		}

	}

	return;
}

/* ================================================== */
void sum_across_surface(struct All_variables *E, float *data, int total)
{
	int j, d;
	float *temp;
	static int been_here = 0;
	static int *processors, nproc;

	static MPI_Comm world, horizon_p;
	static MPI_Group world_g, horizon_g;

	if(been_here == 0)
	{
		processors = (int *)malloc((E->parallel.nprocxy + 2) * sizeof(int));

		nproc = 0;
		for(j = 0; j < E->parallel.nprocxy; j++)
		{
			d = E->parallel.me_loc[3] + j * E->parallel.nprocz;
			processors[nproc] = d;
			nproc++;
		}

		if(nproc > 0)
		{
			world = MPI_COMM_WORLD;
			MPI_Comm_group(world, &world_g);
			MPI_Group_incl(world_g, nproc, processors, &horizon_g);
			MPI_Comm_create(world, horizon_g, &horizon_p);
		}

		been_here++;
	}

	if(nproc > 0)
	{

		temp = (float *)malloc((total + 1) * sizeof(float));
		MPI_Allreduce(data, temp, total, MPI_FLOAT, MPI_SUM, horizon_p);

		for(j = 0; j < total; j++)
		{
			data[j] = temp[j];
		}

		free((void *)temp);

	}

	return;
}

/* ================================================== */
void sum_across_surf_sph1(struct All_variables *E, float *sphc, float *sphs)
{
	int jumpp, total, j, d;
	static float *sphcs, *temp;
	static int been_here = 0;
	static int *processors, nproc;

	static MPI_Comm world, horizon_p;
	static MPI_Group world_g, horizon_g;

	if(been_here == 0)
	{
		processors = (int *)malloc((E->parallel.nprocxy + 2) * sizeof(int));
		temp = (float *)malloc((E->sphere.hindice * 2 + 3) * sizeof(float));
		sphcs = (float *)malloc((E->sphere.hindice * 2 + 3) * sizeof(float));

		nproc = 0;
		for(j = 0; j < E->parallel.nprocxy; j++)
		{
			d = E->parallel.me_loc[3] + j * E->parallel.nprocz;
			processors[nproc] = d;
			nproc++;
		}

		if(nproc > 0)
		{
			world = MPI_COMM_WORLD;
			MPI_Comm_group(world, &world_g);
			MPI_Group_incl(world_g, nproc, processors, &horizon_g);
			MPI_Comm_create(world, horizon_g, &horizon_p);
		}

		been_here++;
	}

	jumpp = E->sphere.hindice;
	total = E->sphere.hindice * 2 + 3;
	for(j = 0; j < E->sphere.hindice; j++)
	{
		sphcs[j] = sphc[j];
		sphcs[j + jumpp] = sphs[j];
	}

	if(nproc > 0)
	{

		MPI_Allreduce(sphcs, temp, total, MPI_FLOAT, MPI_SUM, horizon_p);

		for(j = 0; j < E->sphere.hindice; j++)
		{
			sphc[j] = temp[j];
			sphs[j] = temp[j + jumpp];
		}

	}

	return;
}

/* ==========================================================  */
void gather_TG_to_me0(struct All_variables *E, float *TG)
{
	int i, j, nsl, idb, ii, to_everyone, from_proc, mst, me;

	static float *RG[20];
	static int been_here = 0;
	const float e_16 = 1.e-16;

	MPI_Status status[100];
	//MPI_Status status1;
	MPI_Request request[100];

	if(E->parallel.nprocxy == 1)
		return;

	nsl = E->sphere.nsf + 1;
	me = E->parallel.me;
	if(been_here == 0)
	{
		been_here++;
		for(i = 1; i < E->parallel.nprocxy; i++)
		{
			RG[i] = (float *)malloc((E->sphere.nsf + 1) * sizeof(float));
			RG[i][0] = 0.0;
		}
	}


	idb = 0;
	for(i = 1; i <= E->parallel.nprocxy; i++)
	{
		to_everyone = E->parallel.nprocz * (i - 1) + E->parallel.me_loc[3];

		if(me != to_everyone)
		{						/* send TG */
			idb++;
			mst = me;
			MPI_Isend(TG, nsl, MPI_FLOAT, to_everyone, mst, MPI_COMM_WORLD, &request[idb - 1]);
		}
	}


/* parallel_process_sync(); */

	ii = 0;
	for(i = 1; i <= E->parallel.nprocxy; i++)
	{
		from_proc = E->parallel.nprocz * (i - 1) + E->parallel.me_loc[3];
		if(me != from_proc)
		{						/* me==0 receive all TG and add them up */
			mst = from_proc;
			idb++;
			ii++;
			MPI_Irecv(RG[ii], nsl, MPI_FLOAT, from_proc, mst, MPI_COMM_WORLD, &request[idb - 1]);
		}
	}

	MPI_Waitall(idb, request, status);

	for(i = 1; i < E->parallel.nprocxy; i++)
		for(j = 1; j <= E->sphere.nsf; j++)
		{
			if(fabs(TG[j]) < e_16)
				TG[j] += RG[i][j];
		}

/* parallel_process_sync(); */

	return;
}


/* ==========================================================  */
void propogator_down_process(struct All_variables *E, float *Tadi)
{
	//int i, j, noz, idb, ii, to_proc, from_proc, mst, me;
	int i, j, noz, idb, to_proc, from_proc, mst, me;
	float temp;

	static float *RG[20], *SD;
	static int been_here = 0;
	//const float e_16 = 1.e-16;

	MPI_Status status[100];
	//MPI_Status status1;
	MPI_Request request[100];

	if(E->parallel.nprocz == 1)
		return;

	noz = E->lmesh.noz;
	me = E->parallel.me;
	if(been_here == 0)
	{
		been_here++;
		for(i = 0; i < E->parallel.nprocz; i++)
		{
			RG[i] = (float *)malloc((4) * sizeof(float));
			RG[i][0] = 0.0;
		}
		SD = (float *)malloc((4) * sizeof(float));
		SD[0] = 0.0;
	}

	SD[1] = Tadi[1];
	SD[2] = Tadi[noz];

	idb = 0;
	for(i = E->parallel.nprocz - 1; i >= 0; i--)
	{
		to_proc = i + E->parallel.me_loc[1] * E->parallel.nprocz + E->parallel.me_loc[2] * E->parallel.nprocz * E->parallel.nprocx;
		if(to_proc != me)
		{
			idb++;
			mst = (1 + me) * (to_proc + 1);
			MPI_Isend(SD, 3, MPI_FLOAT, to_proc, mst, MPI_COMM_WORLD, &request[idb - 1]);
		}
	}

	for(i = E->parallel.nprocz - 1; i >= 0; i--)
	{
		from_proc = i + E->parallel.me_loc[1] * E->parallel.nprocz + E->parallel.me_loc[2] * E->parallel.nprocz * E->parallel.nprocx;
		if(from_proc != me)
		{
			idb++;
			mst = (1 + me) * (from_proc + 1);
			j = from_proc % E->parallel.nprocz;
			MPI_Irecv(RG[j], 3, MPI_FLOAT, from_proc, mst, MPI_COMM_WORLD, &request[idb - 1]);
		}
	}

	MPI_Waitall(idb, request, status);

	temp = 0;
	for(i = E->parallel.nprocz - 1 - E->parallel.me_loc[3]; i > 0; i--)
	{
		from_proc = E->parallel.me + i;
		j = from_proc % E->parallel.nprocz;
		temp = temp + RG[j][1];
		if(j == E->parallel.nprocz - 1)
			E->data.T_adi0 = RG[j][2];
	}

	for(j = 1; j <= noz; j++)
		Tadi[j] = Tadi[j] + temp;


	E->data.T_adi1 = 0;
	if(E->parallel.me_loc[3] == 0)
		E->data.T_adi1 = Tadi[1];
	else
	{
		for(i = E->parallel.me_loc[3]; i > 0; i--)
		{
			from_proc = E->parallel.me - i;
			j = from_proc % E->parallel.nprocz;
			if(j < E->parallel.me_loc[3])
				E->data.T_adi1 += RG[j][1];
		}
		E->data.T_adi1 += Tadi[1];
	}

	return;
}

/* ================================================== */
double sum_across_depth(struct All_variables *E, double temp1)
{
	int j, d;
	double temp2;
	static int been_here = 0;
	static int *processors, nproc;

	static MPI_Comm world, vert_p;
	static MPI_Group world_g, vert_g;

	if(been_here == 0)
	{
		processors = (int *)malloc((E->parallel.nprocz + 2) * sizeof(int));

		for(j = 0; j < E->parallel.nprocz; j++)
		{
			d = j + E->parallel.me_loc[1] * E->parallel.nprocz + E->parallel.me_loc[2] * E->parallel.nprocxz;
			processors[j] = d;
		}

		nproc = E->parallel.nprocz;
		world = MPI_COMM_WORLD;
		MPI_Comm_group(world, &world_g);
		MPI_Group_incl(world_g, nproc, processors, &vert_g);
		MPI_Comm_create(world, vert_g, &vert_p);

		been_here++;
	}

	MPI_Allreduce(&temp1, &temp2, 1, MPI_DOUBLE, MPI_SUM, vert_p);

	return (temp2);
}
