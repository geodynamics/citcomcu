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

/* Set up the finite element problem to suit: returns with all memory */
/* allocated, temperature, viscosity, node locations and how to use */
/* them all established. 8.29.92 or 29.8.92 depending on your nationality*/

#include <signal.h>
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#ifdef USE_GGRD
#include "hc.h"
#endif
int Emergency_stop;


void read_instructions(struct All_variables *E, int argc, char **argv)
{
	double start_time;

	//int *temp, i;

	/* =====================================================
	 * Global interruption handling routine defined once here
	 * =====================================================  */

	if(E->parallel.me == 0)
		start_time = CPU_time0();
	Emergency_stop = 0;
	signal(SIGINT, interruption);
	signal(SIGTERM, interruption);

	E->control.PID = get_process_identifier();


	/* ==================================================
	 * Initialize from the command line 
	 * from startup files. (See Parsing.c).
	 * ==================================================  */

	setup_parser(E, argv[1]);

	if(E->parallel.me == 0)
		fprintf(stderr, "ok1\n");

	global_default_values(E);
	if(E->parallel.me == 0)
		fprintf(stderr, "ok2\n");
	read_initial_settings(E);
	if(E->parallel.me == 0)
		fprintf(stderr, "ok3\n");
	(E->problem_derived_values) (E);	/* call this before global_derived_  */
	global_derived_values(E);

	if(E->parallel.me == 0)
		fprintf(stderr, "ok4\n");
	parallel_domain_decomp1(E);

	if(E->parallel.me == 0)
		fprintf(stderr, "ok5\n");
	allocate_common_vars(E);

	if(E->parallel.me == 0)
		fprintf(stderr, "ok6\n");
	(E->problem_allocate_vars) (E);
	(E->solver_allocate_vars) (E);

	if(E->parallel.me == 0)
		fprintf(stderr, "ok6a\n");
	construct_ien(E);
	if(E->parallel.me == 0)
		fprintf(stderr, "ok9\n");
	construct_masks(E);			/* order is important here */
	if(E->parallel.me == 0)
		fprintf(stderr, "ok10\n");
	construct_id(E);
	construct_lm(E);
	if(E->parallel.me == 0)
		fprintf(stderr, "ok11\n");
	construct_sub_element(E);
	if(E->parallel.me == 0)
		fprintf(stderr, "ok12\n");

	node_locations(E);

	if(E->parallel.me == 0)
		fprintf(stderr, "ok7a\n");
	construct_mat_group(E);

	(E->problem_boundary_conds) (E);
	check_bc_consistency(E);

	if(E->parallel.me == 0)
		fprintf(stderr, "ok13\n");

	parallel_shuffle_ele_and_id(E);

	if(E->parallel.me == 0)
		fprintf(stderr, "ok14\n");
	parallel_communication_routs(E);

	construct_shape_functions(E);
	if(E->parallel.me == 0)
		fprintf(stderr, "ok15\n");
	mass_matrix(E);

	if(E->parallel.me == 0)
		fprintf(stderr, "ok16\n");
	(E->problem_initial_fields) (E);	/* temperature/chemistry/melting etc */

	if(E->parallel.me == 0)
		fprintf(stderr, "ok17\n");
	common_initial_fields(E);	/* velocity/pressure/viscosity (viscosity must be done LAST) */
	if(E->parallel.me == 0)
		fprintf(stderr, "ok18\n");

	shutdown_parser(E);

/*	if (E->parallel.me==0)
		fprintf (stderr,"done instruction\n");
	parallel_process_termination();
	exit(8);
*/

	return;
}


/* ===================================
   Functions which set up details 
   common to all problems follow ...
   ===================================  */

void allocate_common_vars(struct All_variables *E)
{
	//int nox, noy, noz, i, j, l, nno_l, npno_l, nozl, nnov_l, nxyz;
	int nox, noy, noz, i, j, l, nxyz;

	E->mesh.fnodal_malloc_size = (E->lmesh.nno + 2) * sizeof(float);
	E->mesh.dnodal_malloc_size = (E->lmesh.nno + 2) * sizeof(double);
	E->mesh.feqn_malloc_size = (E->mesh.nsd * E->lmesh.nno + 2) * sizeof(float);
	E->mesh.deqn_malloc_size = (E->mesh.nsd * E->lmesh.nno + 2) * sizeof(double);

	E->P = (double *)malloc((E->lmesh.npno + 1) * sizeof(double));
	E->S = (double *)malloc((E->lmesh.npno + 1) * sizeof(double));

	E->Have.f = (float *)malloc((E->mesh.noz + 1) * sizeof(float));
	E->Have.F = (float *)malloc((E->mesh.noz + 1) * sizeof(float));
	E->Have.T = (float *)malloc((E->mesh.noz + 1) * sizeof(float));
	E->Have.C = (float *)malloc((E->mesh.noz + 1) * sizeof(float));
	E->Have.Vi = (float *)malloc((E->mesh.noz + 1) * sizeof(float));
	E->Have.Rho = (float *)malloc((E->mesh.noz + 1) * sizeof(float));
	E->Have.vrms = (float *)malloc((E->mesh.noz + 1) * sizeof(float));
	E->Have.Tadi = (float *)malloc((E->mesh.noz + 1) * sizeof(float));

	E->segment.Tz = (float *)malloc((E->lmesh.noz + 1) * sizeof(float));
	E->segment.Vx = (float *)malloc((E->lmesh.noz + 1) * sizeof(float));

	E->F = (double *)malloc((E->mesh.nsd * E->lmesh.nnov + 1) * sizeof(double));
	E->U = (double *)malloc((E->lmesh.neq + 2) * sizeof(double));
	E->T = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	E->C = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	E->CE = (float *)malloc((E->lmesh.nel + 1) * sizeof(float));
	E->buoyancy = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	E->NP = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	E->heatflux = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	E->heatflux_adv = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	E->edot = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));

	E->heating_adi = (float *)malloc((E->lmesh.nel + 1) * sizeof(float));
	E->heating_visc = (float *)malloc((E->lmesh.nel + 1) * sizeof(float));
	E->heating_latent = (float *)malloc((E->lmesh.nel + 1) * sizeof(float));

	E->Fas670 = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	E->Fas670_b = (float *)malloc((E->lmesh.nsf + 1) * sizeof(float));
	E->Fas410 = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	E->Fas410_b = (float *)malloc((E->lmesh.nsf + 1) * sizeof(float));


	E->diffusivity = (float *)malloc((E->lmesh.noz + 1) * sizeof(float));
	E->expansivity = (float *)malloc((E->lmesh.noz + 1) * sizeof(float));


	for(i = 1; i <= E->mesh.nsd; i++)
	{
		E->TB[i] = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
		E->V[i] = (float *)malloc((E->lmesh.nnov + 1) * sizeof(float));
		E->VB[i] = (float *)malloc((E->lmesh.nnov + 1) * sizeof(float));
	}

	E->stress = (float *)malloc((E->lmesh.nsf * 12 + 12) * sizeof(float));
	E->slice.tpg = (float *)malloc((E->lmesh.nsf + 2) * sizeof(float));
	E->slice.tpgb = (float *)malloc((E->lmesh.nsf + 2) * sizeof(float));
	E->slice.vline = (float *)malloc((E->lmesh.nsf + 2) * sizeof(float));
	E->slice.vlinek = (float *)malloc((E->lmesh.nsf + 2) * sizeof(float));
	E->slice.shflux = (float *)malloc((E->lmesh.nsf + 2) * sizeof(float));
	E->slice.bhflux = (float *)malloc((E->lmesh.nsf + 2) * sizeof(float));
	E->slice.cen_mflux = (float *)malloc((E->lmesh.nsf + 2) * sizeof(float));
	E->slice.vxsurf[1] = (float *)malloc((E->lmesh.nsf + 2) * sizeof(float));
	E->slice.vxsurf[2] = (float *)malloc((E->lmesh.nsf + 2) * sizeof(float));

	E->mat = (int *)malloc((E->lmesh.nel + 2) * sizeof(int));

	E->XP[1] = (double *)malloc((E->lmesh.nox + 1) * sizeof(double));
	E->XP[2] = (double *)malloc((E->lmesh.noy + 1) * sizeof(double));
	E->XP[3] = (double *)malloc((E->lmesh.noz + 1) * sizeof(double));
	E->lmesh.rnoz = 50 * E->lmesh.elz + 1;
	E->XRG[3] = (double *)malloc((E->lmesh.rnoz + 1) * sizeof(double));
	E->RG[3] = (int *)malloc((E->lmesh.rnoz + 1) * sizeof(int));

	/* set up memory for different grids  */
	for(i = E->mesh.levmin; i <= E->mesh.levmax; i++)
	{
		for(j = 1; j <= E->mesh.nsd; j++)
			E->XX[i][j] = (float *)malloc((E->lmesh.NNO[i] + 1) * sizeof(float));
		if(E->control.Rsphere)
			for(j = 1; j <= E->mesh.nsd; j++)
				E->SXX[i][j] = (float *)malloc((E->lmesh.NNO[i] + 1) * sizeof(float));

		E->MASS[i] = (float *)malloc((E->lmesh.NNO[i] + 1) * sizeof(float));

		E->ECO[i] = (struct COORD *)malloc((E->lmesh.NNO[i] + 2) * sizeof(struct COORD));
		E->IEN[i] = (struct IEN *)malloc((E->lmesh.NNO[i] + 2) * sizeof(struct IEN));
		E->ID[i] = (struct ID *)malloc((E->lmesh.NNO[i] + 2) * sizeof(struct ID));
		E->GNX[i] = (struct Shape_function_dx *)malloc((E->lmesh.NEL[i] + 2) * sizeof(struct Shape_function_dx));
		E->GDA[i] = (struct Shape_function_dA *)malloc((E->lmesh.NEL[i] + 2) * sizeof(struct Shape_function_dA));
		E->EL[i] = (struct SUBEL *)malloc((E->lmesh.NEL[i] + 2) * sizeof(struct SUBEL));
		E->LMD[i] = (struct LM *)malloc((E->lmesh.NEL[i] + 2) * sizeof(struct LM));

		E->EVI[i] = (float *)malloc((E->lmesh.NEL[i] + 2) * vpoints[E->mesh.nsd] * sizeof(float));

		E->TW[i] = (float *)malloc((E->lmesh.NNO[i] + 2) * sizeof(float));
		E->VI[i] = (float *)malloc((E->lmesh.NNO[i] + 2) * sizeof(float));
		E->NODE[i] = (unsigned int *)malloc((E->lmesh.NNO[i] + 2) * sizeof(unsigned int));
		E->TWW[i] = (struct FNODE *)malloc((E->lmesh.NEL[i] + 2) * sizeof(struct FNODE));

		E->NEI[i].nels = (int *)malloc((E->lmesh.NNO[i] + 2) * sizeof(int));
		E->NEI[i].lnode = (int *)malloc((E->lmesh.NNO[i] + 2) * enodes[E->mesh.nsd] * sizeof(int));
		E->NEI[i].element = (int *)malloc((E->lmesh.NNO[i] + 2) * enodes[E->mesh.nsd] * sizeof(int));

		E->elt_del[i] = (struct EG *)malloc((E->lmesh.NEL[i] + 1) * sizeof(struct EG));

		E->BI[i] = (double *)malloc((E->lmesh.NEQ[i] + 2) * sizeof(double));
		E->BPI[i] = (double *)malloc((E->lmesh.NPNO[i] + 1) * sizeof(double));
		E->control.B_is_good[i] = 0;
	}

	E->temp = (double *)malloc((E->lmesh.NEQ[E->mesh.levmax] + 2) * sizeof(double));
	E->Element = (unsigned int *)malloc((E->lmesh.nel + 2) * sizeof(unsigned int));


	for(i = E->mesh.levmin; i <= E->mesh.levmax; i++)
	{
		nox = E->lmesh.NOX[i];
		noy = E->lmesh.NOY[i];
		noz = E->lmesh.NOZ[i];
		if(E->mesh.nsd == 2)
		{
			nxyz = max(nox, noz);
		}
		else if(E->mesh.nsd == 3)
		{
			nxyz = max(nox * noz, nox * noy);
			nxyz = max(nxyz, noz * noy);
		}

		E->parallel.IDD[i] = (int *)malloc((E->lmesh.NEQ[i] + 2) * sizeof(int));
		E->parallel.NODE[i] = (struct BOUND *)malloc((nxyz + 2) * sizeof(struct BOUND));
		E->parallel.IDPASS[i] = (struct BOUND *)malloc((10) * sizeof(struct BOUND));
		E->parallel.EXCHANGE_NODE[i] = (struct PASS *)malloc((nxyz + 2) * sizeof(struct PASS));
		E->parallel.EXCHANGE_ID[i] = (struct PASS *)malloc((nxyz * E->mesh.nsd + 3) * sizeof(struct PASS));

		if(i == E->mesh.levmax)
		{
			E->sien = (struct SIEN *)malloc((nxyz + 2) * sizeof(struct SIEN));
			E->surf_element = (int *)malloc((nxyz + 2) * sizeof(int));
			E->surf_node = (int *)malloc((E->lmesh.nsf + 2) * sizeof(int));
		}

	}



	for(l = E->mesh.levmin; l <= E->mesh.levmax; l++)
	{
		for(i = 1; i <= E->lmesh.NNO[l]; i++)
		{
			E->NODE[l][i] = (INTX | INTY | INTZ);	/* and any others ... */
			E->VI[l][i] = 1.0;
			E->TW[l][i] = 0.0;
		}

		for(i = 0; i < E->lmesh.NEQ[l]; i++)
		{
			E->BI[l][i] = 0.0;
			E->parallel.IDD[l][i] = 0;
		}
	}

	for(i = 1; i <= E->lmesh.nnov; i++)
		for(j = 1; j <= E->mesh.nsd; j++)
			E->V[j][i] = E->VB[j][i] = 0.0;

	for(i = 1; i <= E->lmesh.nno; i++)
		for(j = 1; j <= E->mesh.nsd; j++)
			E->TB[j][i] = 0.0;

	for(i = 1; i <= E->lmesh.npno; i++)
		E->P[i] = 0.0;

	for(i = 1; i <= E->lmesh.nel; i++)
	{
		E->mat[i] = 1;
		E->heating_visc[i] = 0;
		E->heating_adi[i] = 0;
		E->heating_latent[i] = 1;
	}

	for(i = 1; i <= E->lmesh.noz; i++)
	{
		E->diffusivity[i] = E->expansivity[i] = 1.0;
		E->Have.Tadi[i] = .0;
	}



	for(i = 1; i <= E->lmesh.nno; i++)
		E->T[i] = E->buoyancy[i] = 0.0;
	for(i = 0; i <= E->lmesh.neq + 1; i++)
		E->U[i] = 0.0;

	set_up_nonmg_aliases(E);

	return;
}


/*  =========================================================  */


void interruption(int signal_number)
{
	if(Emergency_stop++)
		exit(0);
	fprintf(stderr, "Cleaning up before exit\n");
	return;
}


void global_default_values(struct All_variables *E)
{
	//FILE *fp;

	/* FIRST: values which are not changed routinely by the user */

	E->control.v_steps_low = 10;
	E->control.v_steps_upper = 1;
	E->control.max_res_red_each_p_mg = 1.0e-3;
	E->control.accuracy = 1.0e-6;
	E->control.vaccuracy = 1.0e-8;
	E->control.true_vcycle = 0;
	E->control.depth_dominated = 0;
	E->control.eqn_zigzag = 0;
	E->control.verbose = 0;		/* debugging/profiles */
	E->control.stokes = 0;
	E->control.restart = 0;

	/* SECOND: values for which an obvious default setting is useful */

	E->control.ORTHO = 1;		/* for orthogonal meshes by default */
	E->control.ORTHOZ = 1;		/* for orthogonal meshes by default */


	E->control.KERNEL = 0;
	E->control.CONVECTION = 0;
	E->control.SLAB = 0;
	E->control.CART2D = 0;
	E->control.CART3D = 0;
	E->control.Rsphere = 0;
	E->control.CART2pt5D = 0;
	E->control.AXI = 0;
	E->control.CONJ_GRAD = 0;
	E->control.NMULTIGRID = 0;
	E->control.EMULTIGRID = 0;
	E->control.COMPRESS = 1;
	E->control.augmented_Lagr = 0;
	E->control.augmented = 0.0;

	/* Default: all optional modules set to `off' */
	E->control.MELTING_MODULE = 0;
	E->control.CHEMISTRY_MODULE = 0;

	E->control.composition = 0;

	E->control.GRID_TYPE = 1;
	E->mesh.hwidth[1] = E->mesh.hwidth[2] = E->mesh.hwidth[3] = 1.0;	/* divide by this one ! */
	E->mesh.magnitude[1] = E->mesh.magnitude[2] = E->mesh.magnitude[3] = 0.0;
	E->mesh.offset[1] = E->mesh.offset[2] = E->mesh.offset[3] = 0.0;

	E->parallel.automa = 0;
	E->parallel.nprocx = 1;
	E->parallel.nprocz = 1;
	E->parallel.nprocy = 1;

	E->mesh.levmax = 0;
	E->mesh.levmin = 0;
	E->mesh.nox = 1;
	E->mesh.nxs = 1;
	E->lmesh.nox = 1;
	E->lmesh.nxs = 1;
	E->mesh.noz = 1;
	E->mesh.nzs = 1;
	E->lmesh.noz = 1;
	E->lmesh.nzs = 1;
	E->mesh.noy = 1;
	E->mesh.nys = 1;
	E->lmesh.noy = 1;
	E->lmesh.nys = 1;

	E->viscosity.guess = 0;
	sprintf(E->viscosity.old_file, "initialize");

	E->control.precondition = 0;	/* for larger visc contrasts turn this back on  */
	E->control.vprecondition = 1;

	E->mesh.toptbc = 1;			/* fixed t */
	E->mesh.bottbc = 1;
	E->mesh.topvbc = 0;			/* stress */
	E->mesh.botvbc = 0;
	E->mesh.sidevbc = 0;
	E->mesh.periodic_x = 0;		/* reflection is default */
	E->mesh.periodic_y = 0;
	E->control.VBXtopval = 0.0;
	E->control.VBYtopval = 0.0;
	E->control.VBXbotval = 0.0;
	E->control.VBYbotval = 0.0;

	E->data.layer_km = 2800.0;	/* Earth, whole mantle defaults */
	E->monitor.length_scale = 2870000.0;
	E->data.grav_acc = 9.81;
	E->data.therm_exp = 3.28e-5;
	E->data.Cp = 1200.0;
	E->data.therm_diff = 8.0e-7;
	E->data.therm_cond = 3.168;
	E->data.density = 3340.0;
	E->data.res_density = 3295.0;	/* density when X = ... */
	E->data.res_density_X = 0.3;
	E->data.melt_density = 2800.0;
	E->data.permeability = 3.0e-10;
	E->data.density_above = 1030.0;	/* sea water */
	E->data.gas_const = 8.3;
	E->data.surf_heat_flux = 4.4e-2;
	E->data.grav_const = 6.673e-11;
	E->data.surf_temp = 0.0;
	E->data.disptn_number = 0.0;
	E->data.youngs_mod = 1.0e11;
	E->data.Te = 0.0;
	E->data.T_sol0 = 1373.0;	/* Dave's values 1991 (for the earth) */
	E->data.Tsurf = 273.0;
	E->data.dTsol_dz = 3.4e-3;
	E->data.dTsol_dF = 440.0;
	E->data.dT_dz = 0.48e-3;
	E->data.delta_S = 250.0;
	E->data.ref_temperature = 2 * 1350.0;	/* fixed temperature ... delta T */

	E->control.Ra_670 = 0.0;
	E->control.Ra_410 = 0.0;

	/* THIRD: you forgot and then went home, let's see if we can help out */

	sprintf(E->control.data_file, "citcom.tmp.%d", getpid());

	E->control.NASSEMBLE = 0;

	E->mesh.layer[1] = E->mesh.layer[2] = E->mesh.layer[3] = 1.0;
	E->monitor.elapsed_time = 0.0;

	return;
}


void global_derived_values(struct All_variables *E)
{
	int d, lx, lz, ly, i, nox, noz, noy;
	char logfile[100];
	FILE *fp;

	/* As early as possible, set up the log file to 
	 * record information about the progress of the 
	 * program as it runs 
	 */
	
#ifdef USE_GZDIR
	sprintf(logfile, "%s/log%d", E->control.data_file, E->parallel.me);
#else
	sprintf(logfile, "%s.log%d", E->control.data_file, E->parallel.me);
#endif
	if((fp = fopen(logfile, "w")) == NULL)
		E->fp = stdout;
	else
		E->fp = fp;

	//fprintf(E->fp, "visc_f %g %g\n", E->data.visc_factor, E->control.Q0);


	if(E->control.NMULTIGRID || E->control.EMULTIGRID)
	{
		E->mesh.levmax = E->mesh.levels - 1;
		E->mesh.nox = E->mesh.mgunitx * (int)pow(2.0, ((double)E->mesh.levmax)) + 1;
		E->mesh.noy = E->mesh.mgunity * (int)pow(2.0, ((double)E->mesh.levmax)) + 1;
		E->mesh.noz = E->mesh.mgunitz * (int)pow(2.0, ((double)E->mesh.levmax)) + 1;
	}

	E->mesh.nnx[1] = E->mesh.nox;
	E->mesh.nnx[3] = E->mesh.noz;
	E->mesh.nnx[2] = E->mesh.noy;
	E->mesh.nex[1] = E->mesh.elx = E->mesh.nox - 1;
	E->mesh.nex[3] = E->mesh.elz = E->mesh.noz - 1;
	E->mesh.nex[2] = E->mesh.ely = max(E->mesh.noy - 1, 1);
	E->mesh.nno = E->mesh.nel = 1;
	E->mesh.nmx = max(E->mesh.nox, max(E->mesh.noz, E->mesh.noy));
	for(d = 1; d <= E->mesh.nsd; d++)
	{
		E->mesh.nno *= E->mesh.nnx[d];
		E->mesh.nel *= E->mesh.nex[d];
	}

	E->mesh.nlno = E->mesh.nno;
	if(E->mesh.periodic_x)
	{
		E->mesh.nlno *= E->mesh.nox - 1;
		E->mesh.nlno /= E->mesh.nox;
	}

	if(E->mesh.periodic_y)
	{
		E->mesh.nlno *= E->mesh.noy - 1;
		E->mesh.nlno /= E->mesh.noy;
	}

	E->mesh.nnov = E->mesh.nno;
	E->mesh.neq = E->mesh.nnov * E->mesh.nsd;

	E->mesh.npno = E->mesh.nel;
	E->mesh.nsf = E->mesh.nox * E->mesh.noy;

	for(i = E->mesh.levmax; i >= E->mesh.levmin; i--)	/* set up dimensions for different grids  */
	{
		if(E->control.NMULTIGRID || E->control.EMULTIGRID)
		{
			nox = E->mesh.mgunitx * (int)pow(2.0, (double)i) + 1;
			noz = E->mesh.mgunitz * (int)pow(2.0, (double)i) + 1;
			noy = E->mesh.mgunity * (int)pow(2.0, (double)i) + 1;
		}
		else
		{
			nox = E->mesh.nox;
			noz = E->mesh.noz;
			noy = E->mesh.noy;
		}

		E->mesh.ELX[i] = nox - 1;
		E->mesh.ELZ[i] = noz - 1;
		E->mesh.ELY[i] = max(noy - 1, 1);
		E->mesh.NNO[i] = nox * noz * noy;
		E->mesh.NLNO[i] = nox * noz * noy;
		E->mesh.NEL[i] = (nox - 1) * (noz - 1) * max((noy - 1), 1);
		E->mesh.NPNO[i] = E->mesh.NEL[i];
		E->mesh.NOX[i] = nox;
		E->mesh.NOZ[i] = noz;
		E->mesh.NOY[i] = noy;
		E->mesh.NMX[i] = max(nox, max(noz, noy));

		if(E->mesh.periodic_x)
		{
			E->mesh.NLNO[i] *= E->mesh.NOX[i] - 1;
			E->mesh.NLNO[i] /= E->mesh.NOX[i];
		}

		if(E->mesh.periodic_y)
		{
			E->mesh.NLNO[i] *= E->mesh.NOY[i] - 1;
			E->mesh.NLNO[i] /= E->mesh.NOY[i];
		}

		E->mesh.NNOV[i] = E->mesh.NNO[i];
		E->mesh.NEQ[i] = E->mesh.nsd * E->mesh.NNOV[i];

		lx--;
		lz--;
		ly--;
	}

	if(E->control.print_convergence && E->parallel.me == 0)
		fprintf(stderr, "Problem has %d x %d x %d nodes\n", E->mesh.nox, E->mesh.noz, E->mesh.noy);



	return;
}


void read_initial_settings(struct All_variables *E)
{
	char tmp_string[100];

	/* first the problem type (defines subsequent behaviour) */
	int m;

	m = E->parallel.me;

	input_string("Problem", E->control.PROBLEM_TYPE, NULL, m);
	if(strcmp(E->control.PROBLEM_TYPE, "convection") == 0)
	{
		E->control.CONVECTION = 1;
		set_convection_defaults(E);
	}

	else if(strcmp(E->control.PROBLEM_TYPE, "convection-chemical") == 0)
	{
		E->control.CONVECTION = 1;
		E->control.CHEMISTRY_MODULE = 1;
		set_convection_defaults(E);
	}

	else
	{
		fprintf(E->fp, "Unable to determine problem type, assuming convection ... \n");
		E->control.CONVECTION = 1;
		set_convection_defaults(E);
	}

	input_string("Geometry", E->control.GEOMETRY, NULL, m);
	if(strcmp(E->control.GEOMETRY, "cart2d") == 0)
	{
		E->control.CART2D = 1;
		set_2dc_defaults(E);
	}
	else if(strcmp(E->control.GEOMETRY, "axi") == 0)
	{
		E->control.AXI = 1;
	}
	else if(strcmp(E->control.GEOMETRY, "cart2pt5d") == 0)
	{
		E->control.CART2pt5D = 1;
		set_2pt5dc_defaults(E);
	}
	else if(strcmp(E->control.GEOMETRY, "cart3d") == 0)
	{
		E->control.CART3D = 1;
		set_3dc_defaults(E);
	}
	else if(strcmp(E->control.GEOMETRY, "Rsphere") == 0)
	{
		E->control.Rsphere = 1;
		set_3ds_defaults(E);
	}
	else
	{
		fprintf(E->fp, "Unable to determine geometry, assuming cartesian 2d ... \n");
		E->control.CART2D = 1;
		set_2dc_defaults(E);
	}

	input_string("Solver", E->control.SOLVER_TYPE, NULL, m);
	if(strcmp(E->control.SOLVER_TYPE, "cgrad") == 0)
	{
		E->control.CONJ_GRAD = 1;
		set_cg_defaults(E);
	}
	else if(strcmp(E->control.SOLVER_TYPE, "multigrid") == 0)
	{
		E->control.NMULTIGRID = 1;
		set_mg_defaults(E);
	}
	else if(strcmp(E->control.SOLVER_TYPE, "multigrid-el") == 0)
	{
		E->control.EMULTIGRID = 1;
		set_mg_defaults(E);
	}
	else
	{
		if(E->parallel.me == 0)
			fprintf(stderr, "Unable to determine how to solve, specify Solver=VALID_OPTION \n");
		exit(0);
	}


	/* admin */

	/* Information on which files to print, which variables of the flow to calculate and print.
	 * Default is no information recorded (apart from special things for given applications.
	 */

	input_string("datatypes", E->control.which_data_files, "", m);
	input_string("averages", E->control.which_horiz_averages, "", m);
	input_string("timelog", E->control.which_running_data, "", m);
	input_string("observables", E->control.which_observable_data, "", m);

	input_string("datafile", E->control.data_file, "initialize", m);
	input_string("process_command", E->control.output_written_external_command, "", m);

#ifdef USE_GZDIR
	input_boolean("gzdir",&(E->control.gzdir),"on",m);
#endif
	input_string("use_scratch", tmp_string, "local", m);
	if(strcmp(tmp_string, "local") == 0)
		strcpy(E->control.data_file2, E->control.data_file);
	else
		sprintf(E->control.data_file2, "/scratch_%s/%s/%s", E->parallel.machinename, tmp_string, E->control.data_file);

	input_boolean("AVS", &(E->control.AVS), "off", m);
	input_boolean("CONMAN", &(E->control.CONMAN), "off", m);

	if(E->control.NMULTIGRID || E->control.EMULTIGRID)
	{
		input_int("mgunitx", &(E->mesh.mgunitx), "1", m);
		input_int("mgunitz", &(E->mesh.mgunitz), "1", m);
		input_int("mgunity", &(E->mesh.mgunity), "1", m);
		input_int("mgunitxl", &(E->lmesh.mgunitx), "1", m);
		input_int("mgunitzl", &(E->lmesh.mgunitz), "1", m);
		input_int("mgunityl", &(E->lmesh.mgunity), "1", m);
		input_int("levels", &(E->mesh.levels), "0", m);
	}

	input_int("restart", &(E->control.restart), "0", m);

	input_boolean("regular_grid", &(E->control.ORTHOZ), "off", m);

	input_int("stokes_flow_only", &(E->control.stokes), "0", m);

	input_boolean("node_assemble", &(E->control.NASSEMBLE), "off", m);
	/* general mesh structure */

	input_boolean("parallel_auto", &(E->parallel.automa), "off", m);
	if(E->parallel.automa == 0)
	{
		input_int("nprocx", &(E->parallel.nprocx), "1", m);
		input_int("nprocz", &(E->parallel.nprocz), "1", m);
		input_int("nprocy", &(E->parallel.nprocy), "1", m);
	}

	input_boolean("verbose", &(E->control.verbose), "off", m);
	input_boolean("see_convergence", &(E->control.print_convergence), "off", m);
	input_boolean("COMPRESS", &(E->control.COMPRESS), "on", m);
	input_float("sobtol", &(E->control.sob_tolerance), "0.0001", m);
	input_float("T_interior_max", &(E->monitor.T_interior_max), "1.5", m);

	input_int("obs_maxlongk", &(E->slice.maxlong), "100,1", m);
	input_int("obs_minlongk", &(E->slice.minlong), "1,1", m);

	/* for layers    */
	E->viscosity.zlm = 1.0;
	E->viscosity.z410 = 1.0;
	E->viscosity.zlith = 0.0;

	if(E->control.CART3D)
	{
		input_float("z_lmantle", &(E->viscosity.zlm), "1.0", m);
		input_float("z_410", &(E->viscosity.z410), "1.0", m);
		input_float("z_lith", &(E->viscosity.zlith), "0.0", m);
	}
	else if(E->control.Rsphere)
	{
		input_float("r_lmantle", &(E->viscosity.zlm), "1.0", m);
		input_float("r_410", &(E->viscosity.z410), "1.0", m);
		input_float("r_lith", &(E->viscosity.zlith), "0.0", m);
	}


#ifdef USE_GGRD
	/* ggrd control */
	ggrd_init_master(&(E->control.ggrd));
	input_boolean("ggrd_tinit",&(E->control.ggrd.use_temp),"off", m);
	input_double("ggrd_tinit_scale",&(E->control.ggrd.temp.scale),"1.0", m);
	input_boolean("ggrd_scale_with_prem",&(E->control.ggrd.temp.scale_with_prem),"off", m);
	input_string("ggrd_tinit_gfile",E->control.ggrd.temp.gfile,"", m);
	input_string("ggrd_tinit_dfile",E->control.ggrd.temp.dfile,"", m);
	input_double("ggrd_tinit_offset",&(E->control.ggrd.temp.offset),"0.0", m);
	/* comp */
	input_boolean("ggrd_cinit",&(E->control.ggrd.use_comp),"off", m);
	input_double("ggrd_cinit_scale",&(E->control.ggrd.comp.scale),"1.0", m);
	input_string("ggrd_cinit_gfile",E->control.ggrd.comp.gfile,"", m);
	input_string("ggrd_cinit_dfile",E->control.ggrd.comp.dfile,"", m);
	input_double("ggrd_cinit_offset",&(E->control.ggrd.comp.offset),"0.0", m);
	/* slab slice handling */
	input_boolean("slab_slice",&(E->control.slab_slice),"off", m);
	input_float("slab_theta_bound",&(E->control.slab_theta_bound),"1.0", m);
	
#endif

	E->control.transT670 = 1500;
	E->control.transT410 = 1500;

	input_float("Ra_670", &(E->control.Ra_670), "0.0", m);
	input_float("clapeyron670", &(E->control.clapeyron670), "0.0", m);
	input_float("transT670", &(E->control.transT670), "0.0", m);
	input_float("width670", &(E->control.width670), "0.0", m);

	input_float("Ra_410", &(E->control.Ra_410), "0.0", m);
	input_float("clapeyron410", &(E->control.clapeyron410), "0.0", m);
	input_float("transT410", &(E->control.transT410), "0.0", m);
	input_float("width410", &(E->control.width410), "0.0", m);

	input_int("ll_max", &(E->sphere.llmax), "1", m);
	input_int("nlong", &(E->sphere.noy), "1", m);
	input_int("nlati", &(E->sphere.nox), "1", m);


	input_int("topvbc", &(E->mesh.topvbc), "0", m);
	input_int("botvbc", &(E->mesh.botvbc), "0", m);
	input_int("sidevbc", &(E->mesh.sidevbc), "0", m);

	input_boolean("periodicx", &(E->mesh.periodic_x), "off", m);
	input_boolean("periodicy", &(E->mesh.periodic_y), "off", m);
	input_boolean("depthdominated", &(E->control.depth_dominated), "off", m);
	input_boolean("eqnzigzag", &(E->control.eqn_zigzag), "off", m);
	input_boolean("eqnviscosity", &(E->control.eqn_viscosity), "off", m);

	input_float("topvbxval", &(E->control.VBXtopval), "0.0", m);
	input_float("botvbxval", &(E->control.VBXbotval), "0.0", m);
	input_float("topvbyval", &(E->control.VBYtopval), "0.0", m);
	input_float("botvbyval", &(E->control.VBYbotval), "0.0", m);

	input_int("toptbc", &(E->mesh.toptbc), "1", m);
	input_int("bottbc", &(E->mesh.bottbc), "1", m);
	input_float("toptbcval", &(E->control.TBCtopval), "0.0", m);
	input_float("bottbcval", &(E->control.TBCbotval), "1.0", m);

	input_float("plate_velocity", &(E->control.plate_vel), "0.0", m);
	input_float("plate_age", &(E->control.plate_age), "0.0", m);
	input_float("plume_radius", &(E->segment.plume_radius), "0.0", m);
	input_float("plume_DT", &(E->segment.plume_DT), "0.0", m);
	input_float("plume_x", &(E->segment.plume_coord[1]), "0.0", m);
	input_float("plume_y", &(E->segment.plume_coord[2]), "0.0", m);
	input_float("plume_z", &(E->segment.plume_coord[3]), "0.0", m);

	/* check input temp / comp ranges */
	input_boolean("check_t_irange", &(E->control.check_t_irange), "on", m);
	input_boolean("check_c_irange", &(E->control.check_c_irange), "on", m);

	input_string("gridxfile", E->mesh.gridfile[1], " ", m);
	input_string("gridyfile", E->mesh.gridfile[2], " ", m);
	input_string("gridzfile", E->mesh.gridfile[3], " ", m);


	input_float("dimenx", &(E->mesh.layer[1]), "nodefault", m);
	input_float("dimeny", &(E->mesh.layer[2]), "nodefault", m);
	input_float("dimenz", &(E->mesh.layer[3]), "nodefault", m);

	input_float("radius_inner", &(E->sphere.ri), "nodefault", m);
	input_float("radius_outer", &(E->sphere.ro), "nodefault", m);
	input_float("theta_north", &(E->sphere.ti), "nodefault", m);
	input_float("theta_south", &(E->sphere.to), "nodefault", m);
	input_float("fi_west", &(E->sphere.fi), "nodefault", m);
	input_float("fi_east", &(E->sphere.fo), "nodefault", m);
	E->sphere.corner[1][1] = E->sphere.ti;
	E->sphere.corner[2][1] = E->sphere.to;
	E->sphere.corner[1][2] = E->sphere.fi;
	E->sphere.corner[2][2] = E->sphere.fo;
	E->sphere.corner[1][3] = E->sphere.ri;
	E->sphere.corner[2][3] = E->sphere.ro;

	input_int("nodex", &(E->mesh.nox), "nodefault,1,nomax", m);
	input_int("nodez", &(E->mesh.noz), "nodefault,1,nomax", m);
	input_int("nodey", &(E->mesh.noy), "1,1,nomax", m);
	input_boolean("aug_lagr", &(E->control.augmented_Lagr), "off", m);
	input_double("aug_number", &(E->control.augmented), "0.0", m);

	input_boolean("sdepv_print_convergence",&(E->control.sdepv_print_convergence),"off", m);
	//input_float("tole_compressibility", &(E->control.tole_comp), "0.0", m);
	input_boolean("orthogonal", &(E->control.ORTHO), "on", m);
	input_boolean("crust", &(E->control.crust), "off", m);
	input_float("crust_width", &(E->crust.width), "0.0", m);

	input_int("storage_spacing", &(E->control.record_every), "10", m);
	input_int("storage_always_before", &(E->control.record_all_until), "5", m);

	input_boolean("precond", &(E->control.precondition), "off", m);
	input_boolean("vprecond", &(E->control.vprecondition), "on", m);
	input_int("mg_cycle", &(E->control.mg_cycle), "2,0,nomax", m);
	input_int("down_heavy", &(E->control.down_heavy), "1,0,nomax", m);
	input_int("up_heavy", &(E->control.up_heavy), "1,0,nomax", m);
	input_double("accuracy", &(E->control.accuracy), "1.0e-4,0.0,1.0", m);
	input_int("viterations", &(E->control.max_vel_iterations), "250,0,nomax", m);


	input_int("vhighstep", &(E->control.v_steps_high), "1,0,nomax", m);
	input_int("vlowstep", &(E->control.v_steps_low), "250,0,nomax", m);
	input_int("vupperstep", &(E->control.v_steps_upper), "1,0,nomax", m);
	input_int("piterations", &(E->control.p_iterations), "100,0,nomax", m);
	input_int("maxsamevisc", &(E->control.max_same_visc), "25,0,nomax", m);

	/* data section */
	E->data.therm_exp_factor = 1.0;
	E->data.therm_diff_factor = 1.0;
	E->data.visc_factor = 1.0;

	input_float("ReferenceT", &(E->data.ref_temperature), "2600.0", m);
	input_float("Q0", &(E->control.Q0), "0.0", m);

//  input_float("layerd",&(E->data.layer_km),"2870000.0", m);

	if(E->control.CART3D)
		input_float("layerd", &(E->monitor.length_scale), "2870000.0", m);
	else if(E->control.Rsphere)
		input_float("radius", &(E->monitor.length_scale), "6370000.0", m);


	

	input_float("gravacc", &(E->data.grav_acc), "9.81", m);
	input_float("thermexp", &(E->data.therm_exp), "3.28e-5", m);
	input_float("thermexp_factor", &(E->data.therm_exp_factor), "1.0", m);
	input_float("visc_factor", &(E->data.visc_factor), "1.0", m);
	input_float("thermdiff_factor", &(E->data.therm_diff_factor), "1.0", m);
	input_float("cp", &(E->data.Cp), "1200.0", m);
	input_float("thermdiff", &(E->data.therm_diff), "8.0e-7", m);
	input_float("thermcond", &(E->data.therm_cond), "3.168", m);
	input_float("density", &(E->data.density), "3340.0", m);
	input_float("mdensity", &(E->data.melt_density), "2800.0", m);
	input_float("wdensity", &(E->data.density_above), "1030.0", m);
	input_float("rdensity", &(E->data.res_density), "3295.0", m);
	input_float("heatflux", &(E->data.surf_heat_flux), "4.4e-2", m);
	input_float("refvisc", &(E->data.ref_viscosity), "nodefault", m);
	input_float("meltvisc", &(E->data.melt_viscosity), "nodefault", m);
	input_float("surf_temp", &(E->data.surf_temp), "0.0", m);
	input_float("dissipation_number", &(E->data.disptn_number), "0.0", m);
	input_float("youngs", &(E->data.youngs_mod), "1.0e11", m);
	input_float("Te", &(E->data.Te), "0.0", m);
	input_float("Tsol0", &(E->data.T_sol0), "1373.0", m);
	input_float("dTsoldz", &(E->data.dTsol_dz), "3.4e-3", m);
	input_float("dTsoldF", &(E->data.dTsol_dF), "440.0", m);
	input_float("dTdz", &(E->data.dT_dz), "0.48e-3", m);
	input_float("deltaS", &(E->data.delta_S), "250.0", m);
	input_float("gasconst", &(E->data.gas_const), "8.3",m);	/* not much cause to change these !
	input_float("gravconst", &(E->data.grav_const), "6.673e-11", m);
	input_float("permeability", &(E->data.permeability), "3.0e-10", m);

	/* scaling */
	E->monitor.time_scale = pow(E->monitor.length_scale, 2.0) /E->data.therm_diff;
	/* million years */
	E->monitor.time_scale_ma = E->monitor.time_scale/(3600.0 * 24.0 * 365.25 * 1.0e6);
	
	E->monitor.velo_scale = E->data.therm_diff / (E->monitor.length_scale);
	E->monitor.tau_scale = (E->data.ref_viscosity/E->monitor.time_scale); /* scaling stress */


	(E->problem_settings) (E);

	viscosity_parameters(E);

	return;
}

void check_bc_consistency(struct All_variables *E)
{
	int i, lev;

	for(i = 1; i <= E->lmesh.nno; i++)
	{
		if((E->node[i] & VBX) && (E->node[i] & SBX))
			printf("Inconsistent x velocity bc at %d\n", i);
		if((E->node[i] & VBZ) && (E->node[i] & SBZ))
			printf("Inconsistent z velocity bc at %d\n", i);
		if((E->node[i] & VBY) && (E->node[i] & SBY))
			printf("Inconsistent y velocity bc at %d\n", i);
		if((E->node[i] & TBX) && (E->node[i] & FBX))
			printf("Inconsistent x temperature bc at %d\n", i);
		if((E->node[i] & TBZ) && (E->node[i] & FBZ))
			printf("Inconsistent z temperature bc at %d\n", i);
		if((E->node[i] & TBY) && (E->node[i] & FBY))
			printf("Inconsistent y temperature bc at %d\n", i);
	}

	for(lev = E->mesh.levmin; lev <= E->mesh.levmax; lev++)
		for(i = 1; i <= E->lmesh.NNO[lev]; i++)
		{
			if((E->NODE[lev][i] & VBX) && (E->NODE[lev][i] & SBX))
				printf("Inconsistent x velocity bc at %d,%d\n", lev, i);
			if((E->NODE[lev][i] & VBZ) && (E->NODE[lev][i] & SBZ))
				printf("Inconsistent z velocity bc at %d,%d\n", lev, i);
			if((E->NODE[lev][i] & VBY) && (E->NODE[lev][i] & SBY))
				printf("Inconsistent y velocity bc at %d,%d\n", lev, i);
			/* Tbc's not applicable below top level */ }

	return;
}

void set_up_nonmg_aliases(struct All_variables *E)
{								/* Aliases for functions only interested in the highest mg level */
	int i;

	E->eco = E->ECO[E->mesh.levmax];
	E->ien = E->IEN[E->mesh.levmax];
	E->id = E->ID[E->mesh.levmax];
	E->lm = E->LMD[E->mesh.levmax];
	E->Vi = E->VI[E->mesh.levmax];
	E->EVi = E->EVI[E->mesh.levmax];
	E->node = E->NODE[E->mesh.levmax];
	E->tw = E->TW[E->mesh.levmax];
	E->Mass = E->MASS[E->mesh.levmax];
	E->gDA = E->GDA[E->mesh.levmax];
	E->gNX = E->GNX[E->mesh.levmax];

	for(i = 1; i <= E->mesh.nsd; i++)
		E->X[i] = E->XX[E->mesh.levmax][i];
	if(E->control.Rsphere)
		for(i = 1; i <= E->mesh.nsd; i++)
			E->SX[i] = E->SXX[E->mesh.levmax][i];

	return;
}

void report(struct All_variables *E, char *string)
{
	if(E->control.verbose && E->parallel.me == 0)
	{
		fprintf(stderr, "%s\n", string);
		fflush(stderr);
	}
	return;
}

void record(struct All_variables *E, char *string)
{
	if(E->control.verbose)
	{
		fprintf(E->fp, "%s\n", string);
		fflush(E->fp);
	}

	return;
}



/* =============================================================
   Initialize values which are not problem dependent.
   NOTE: viscosity may be a function of all previous
   input fields (temperature, pressure, velocity, chemistry) and 
   so is always to be done last.
   ============================================================= */


void common_initial_fields(struct All_variables *E)
{
	report(E, "Initialize pressure field");
	initial_pressure(E);
	report(E, "Initialize velocity field");
	initial_velocity(E);
	report(E, "Initialize viscosity field");
	get_viscosity_option(E);

	return;
}

/* ========================================== */

void initial_pressure(struct All_variables *E)
{
	//int i, node, ii;
	int i;

	for(i = 1; i <= E->lmesh.npno; i++)
		E->P[i] = 0.0;

	return;
}

void initial_velocity(struct All_variables *E)
{
	//int i, node, ii;
	int i;

	for(i = 1; i <= E->lmesh.nnov; i++)
	{
		E->V[1][i] = 0.0;
		E->V[2][i] = 0.0;
		E->V[3][i] = 0.0;
	}

	return;
}
