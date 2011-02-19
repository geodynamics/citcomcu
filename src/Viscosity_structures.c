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
#ifdef CITCOM_ALLOW_ANISOTROPIC_VISC
#include "anisotropic_viscosity.h"
#endif

static void visc_from_B(struct All_variables *, float *, float *, int );
static void visc_from_C(struct All_variables *, float *, float *, int );
void *safe_malloc (size_t );

void viscosity_parameters(struct All_variables *E)
{
	int i, l;
	float temp;
	int m = E->parallel.me ;
	/* default values .... */
	E->viscosity.update_allowed = 0;
	E->viscosity.SDEPV = E->viscosity.TDEPV = E->viscosity.BDEPV = 
	  E->viscosity.CDEPV = 0;
	E->viscosity.EXPX = 0;
	E->viscosity.SMOOTH = 0;

	input_boolean("plasticity_dimensional",&(E->viscosity.plasticity_dimensional),"on",m);

	for(i = 0; i < CITCOM_CU_VISC_MAXLAYER; i++)
	{
		E->viscosity.N0[i] = 1.0;
		E->viscosity.T[i] = 0.0;
		E->viscosity.Z[i] = 0.0;
		E->viscosity.E[i] = 0.0;

		
		/* 
		   byerlee or plasticity law 
		*/
		if(E->viscosity.plasticity_dimensional){
		  /* for byerlee, dimensional values are used */
		  E->viscosity.abyerlee[i]=2.1e4; /* t_y = (a * z[m] + b) *p 
						     whereas a is rho * g * 0.6(from friction)
						     and b is 60MPa cohesion
						  */
		  E->viscosity.bbyerlee[i]=6e7;
		  
		  E->viscosity.lbyerlee[i]=.1;
		}else{
		  /* for plasticity: t_y = min (a + (1-x) * b, l) */
		  E->viscosity.abyerlee[i]=0.0;
		  E->viscosity.bbyerlee[i]=1.0;
		  
		  E->viscosity.lbyerlee[i]=1e20;
		}
		
		/* comp dep visc */
		E->viscosity.pre_comp[i] = 1.0;


		E->viscosity.zbase_layer[i] = -999; /* not assigned by default */
	} /* end layer loop */

	/* read in information */
	input_int("rheol", &(E->viscosity.RHEOL), "essential",m);
	input_int("num_mat", &(E->viscosity.num_mat), "1",m);


	/* four layers    */
	E->viscosity.zlm = 1.0;
	E->viscosity.z410 = 1.0;
	E->viscosity.zlith = 0.0;
	if(E->control.CART3D)	/* defaults could be betters */
	{
		input_float("z_lmantle", &(E->viscosity.zlm), "1.0", m);
		input_float("z_410", &(E->viscosity.z410), "1.0", m);
		input_float("z_lith", &(E->viscosity.zlith), "0.0", m);
		input_float_vector("z_layer",E->viscosity.num_mat,(E->viscosity.zbase_layer),m);
	}
	else if(E->control.Rsphere)
	{
		input_float("r_lmantle", &(E->viscosity.zlm), "1.0", m);
		input_float("r_410", &(E->viscosity.z410), "1.0", m);
		input_float("r_lith", &(E->viscosity.zlith), "0.0", m);
		input_float_vector("r_layer",E->viscosity.num_mat,(E->viscosity.zbase_layer),m);
	}

	if((fabs(E->viscosity.zbase_layer[0]+999) < 1e-5) &&
	   (fabs(E->viscosity.zbase_layer[1]+999) < 1e-5)) {
	  /* no z_layer input found */	  
	  if(E->viscosity.num_mat != 4)
            myerror("either use z_layer for non dim layer depths, or set num_mat to four",E);
	  E->viscosity.zbase_layer[0] = E->viscosity.zlith;
	  E->viscosity.zbase_layer[1] = E->viscosity.z410;
	  E->viscosity.zbase_layer[2] = E->viscosity.zlm;
	  if(E->control.Rsphere)
	    E->viscosity.zbase_layer[3] = 0.55;
	  else
	    E->viscosity.zbase_layer[3] = 0.;
	}

	input_float_vector("viscT", E->viscosity.num_mat, (E->viscosity.T),m);	/* redundant */
	input_float_vector("viscT1", E->viscosity.num_mat, (E->viscosity.T),m);
	input_float_vector("viscZ", E->viscosity.num_mat, (E->viscosity.Z),m);
	input_float_vector("viscE", E->viscosity.num_mat, (E->viscosity.E),m);
	input_float_vector("visc0", E->viscosity.num_mat, (E->viscosity.N0),m);	/* redundant */
	input_float_vector("viscN0", E->viscosity.num_mat, (E->viscosity.N0),m);

	input_boolean("TDEPV", &(E->viscosity.TDEPV), "on",m);
	input_boolean("SDEPV", &(E->viscosity.SDEPV), "off",m);
	input_boolean("BDEPV",&(E->viscosity.BDEPV),"off",m);
	input_boolean("CDEPV",&(E->viscosity.CDEPV),"off",m);

	/* plasticity offset viscosity */
	input_float("plasticity_viscosity_offset", &(E->viscosity.plasticity_viscosity_offset),"0.0",m);

	/* 
	   byerlee 
	*/
	input_float_vector("abyerlee",E->viscosity.num_mat,
			   (E->viscosity.abyerlee),m);
	input_float_vector("bbyerlee",E->viscosity.num_mat,
			   (E->viscosity.bbyerlee),m);
	input_float_vector("lbyerlee",E->viscosity.num_mat,
			   (E->viscosity.lbyerlee),m);
	
	
	/* 1: transition 0: min/max transition */
	input_boolean("plasticity_trans",&(E->viscosity.plasticity_trans),"on",m);
	/* SRW */
	input_boolean("psrw",&(E->viscosity.psrw),"off",m);
	/* 
	   
	*/
	input_int("allow_anisotropic_viscosity",&(E->viscosity.allow_anisotropic_viscosity),"0",m);
#ifndef CITCOM_ALLOW_ANISOTROPIC_VISC 
	if(E->viscosity.allow_anisotropic_viscosity){ /* error */
	  fprintf(stderr,"error: allow_anisotropic_viscosity is not zero, but code not compiled with CITCOM_ALLOW_ANISOTROPIC_VISC\n");
	  parallel_process_termination();
	}
#else
	if(E->viscosity.allow_anisotropic_viscosity){ /* read additional
							 parameters for
							 anisotropic
							 viscosity */
	  input_int("anisotropic_init",&(E->viscosity.anisotropic_init),"0",m); /* 0: isotropic
										   1: random
										   2: read in director orientation
										      and log10(eta_s/eta) 
										   3: align with velocity, use ani_vis2_factor for eta_s/eta
										   4: align with ISA, use ani_vis2_factor for eta_s/eta
										   5: align mixed depending on deformation state, use ani_vis2_factor for eta_s/eta
										   6: use radially aligned director and taper eta_s/eta from base (1) to top of layer (ani_vis2_factor)
										   
										*/
	  input_string("anisotropic_init_dir",(E->viscosity.anisotropic_init_dir),"",m); /* directory
											    for
											    ggrd
											    type
											    init */
	  input_int("anivisc_layer",&(E->viscosity.anivisc_layer),"1",m); /* >0: assign to layers on top of anivisc_layer
									     <0: assign to layer = anivisc_layer
									  */
	  if((E->viscosity.anisotropic_init == 6) && (E->viscosity.anivisc_layer >= 0))
	    myerror("anisotropic init mode 6 requires selection of layer where anisotropy applies",E);
	  
	  input_boolean("anivisc_start_from_iso",
			&(E->viscosity.anivisc_start_from_iso),"on",m); /* start
									   from
									   isotropic
									   solution? */
	  if(!E->viscosity.anivisc_start_from_iso)
	    if(E->viscosity.anisotropic_init == 3){
	      if(E->parallel.me == 0)fprintf(stderr,"WARNING: overriding anisotropic first step for ani init mode\n");
	      E->viscosity.anivisc_start_from_iso = TRUE;
	    }
	  /* ratio between weak and strong direction */
	  input_double("ani_vis2_factor",&(E->viscosity.ani_vis2_factor),"1.0",m);


	}
#endif



	/* composition factors */
	input_float_vector("pre_comp",2,(E->viscosity.pre_comp),m);



	input_float("sdepv_misfit", &(E->viscosity.sdepv_misfit), "0.001",m);
	input_float_vector("sdepv_expt", E->viscosity.num_mat, (E->viscosity.sdepv_expt),m);
	input_float_vector("sdepv_trns", E->viscosity.num_mat, (E->viscosity.sdepv_trns),m);

	/* iteration damping for alpha < 1 */
	input_float("sdepv_iter_damp", &(E->viscosity.sdepv_iter_damp), "1.0",m);

	input_boolean("TDEPV_AVE", &(E->viscosity.TDEPV_AVE), "off",m);
	input_boolean("VFREEZE", &(E->viscosity.FREEZE), "off",m);
	input_boolean("VMAX", &(E->viscosity.MAX), "off",m);
	input_boolean("VMIN", &(E->viscosity.MIN), "off",m);
	input_boolean("VISC_UPDATE", &(E->viscosity.update_allowed), "on",m);

	input_float("freeze_thresh", &(E->viscosity.freeze_thresh), "0.0",m);
	input_float("freeze_value", &(E->viscosity.freeze_value), "1.0",m);
	input_float("visc_max", &(E->viscosity.max_value), "nodefault",m);
	input_float("visc_min", &(E->viscosity.min_value), "nodefault",m);

	input_boolean("VISC_GUESS", &(E->viscosity.guess), "off",m);
	input_string("visc_old_file", E->viscosity.old_file, " ",m);

	return;
}

void get_viscosity_option(struct All_variables *E)
{
	/* general, essential default */
  int m = E->parallel.me ;

  input_string("Viscosity", E->viscosity.STRUCTURE, NULL,m);	/* Which form of viscosity */
  
  input_boolean("VISC_EQUIVDD", &(E->viscosity.EQUIVDD), "off",m);	/* Whether to average it */
  input_int("equivdd_opt", &(E->viscosity.equivddopt), "1",m);
  input_int("equivdd_x", &(E->viscosity.proflocx), "1",m);
  input_int("equivdd_y", &(E->viscosity.proflocy), "1",m);
  
  input_int("update_every_steps", &(E->control.KERNEL), "1",m);
  
  input_boolean("VISC_SMOOTH", &(E->viscosity.SMOOTH), "off",m);
  input_int("visc_smooth_cycles", &(E->viscosity.smooth_cycles), "0",m);
  
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

#ifdef CITCOM_ALLOW_ANISOTROPIC_VISC
	if(E->viscosity.allow_anisotropic_viscosity){
	  if(!E->viscosity.anisotropic_viscosity_init)
	    set_anisotropic_viscosity_at_element_level(E,1);
	  else
	    set_anisotropic_viscosity_at_element_level(E,0);
	}
#endif
	

	if(E->viscosity.TDEPV)
		visc_from_T(E, visc, evisc, propogate);
	else
		visc_from_mat(E, visc, evisc);
#ifdef USE_GGRD
	/* pre-factor applies here */
	if(E->control.ggrd.mat_control != 0){
	  ggrd_read_mat_from_file(E);
	  for(i = 1; i <= E->lmesh.nel; i++){
	    for(j = 1; j <= vpoints[E->mesh.nsd]; j++){
	      evisc[(i - 1) * vpoints[E->mesh.nsd] + j] *= E->VIP[i];
	    }
	  }
	}
#endif	

	if(E->viscosity.CDEPV)
	  visc_from_C(E, visc, evisc, propogate);

	if(E->viscosity.SDEPV)
		visc_from_S(E, visc, evisc, propogate);

	if(E->viscosity.BDEPV)
	  visc_from_B(E, visc, evisc, propogate);


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

#ifdef USE_GZDIR
	/* this is much preferred over v_to_nodes */
	visc_from_gint_to_nodes(E,evisc,visc,E->mesh.levmax);
#ifdef CITCOM_ALLOW_ANISOTROPIC_VISC /* allow for anisotropy */
	if(E->viscosity.allow_anisotropic_viscosity){
	  visc_from_gint_to_nodes(E,E->EVI2[E->mesh.levmax], E->VI2[E->mesh.levmax],E->mesh.levmax);
	  visc_from_gint_to_nodes(E,E->EVIn1[E->mesh.levmax], E->VIn1[E->mesh.levmax],E->mesh.levmax);
	  visc_from_gint_to_nodes(E,E->EVIn2[E->mesh.levmax], E->VIn2[E->mesh.levmax],E->mesh.levmax);
	  visc_from_gint_to_nodes(E,E->EVIn3[E->mesh.levmax], E->VIn3[E->mesh.levmax],E->mesh.levmax);
	  normalize_director_at_nodes(E,E->VIn1[E->mesh.levmax],E->VIn2[E->mesh.levmax],E->VIn3[E->mesh.levmax],E->mesh.levmax);
	}
#endif
#endif

	/* v_to_nodes(E,evisc,visc,E->mesh.levmax);  */
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

#ifdef CITCOM_ALLOW_ANISOTROPIC_VISC /* allow for anisotropy */
		if(E->viscosity.allow_anisotropic_viscosity){
		  p_to_centres(E, E->EVI2[E->mesh.levmax], ViscCentre, E->mesh.levmax);p_to_nodes(E, ViscCentre, E->EVI2[E->mesh.levmax], E->mesh.levmax);
		  p_to_centres(E, E->EVIn1[E->mesh.levmax], ViscCentre, E->mesh.levmax);p_to_nodes(E, ViscCentre, E->EVIn1[E->mesh.levmax], E->mesh.levmax);
		  p_to_centres(E, E->EVIn2[E->mesh.levmax], ViscCentre, E->mesh.levmax);p_to_nodes(E, ViscCentre, E->EVIn2[E->mesh.levmax], E->mesh.levmax);
		  p_to_centres(E, E->EVIn3[E->mesh.levmax], ViscCentre, E->mesh.levmax);p_to_nodes(E, ViscCentre, E->EVIn3[E->mesh.levmax], E->mesh.levmax);
		  normalize_director_at_gint(E,E->EVIn1[E->mesh.levmax],E->EVIn2[E->mesh.levmax],E->EVIn3[E->mesh.levmax],E->mesh.levmax);
		}
#endif

		
	}

	free((void *)ViscCentre);

	return;
}

void visc_from_mat(struct All_variables *E, float *Eta, float *EEta)
{
	//int i, j, k, l, z, jj, kk;
	int i, jj,l;
	for(i = 1; i <= E->lmesh.nel; i++){
	  l = E->mat[i] - 1;
	  for(jj = 1; jj <= vpoints[E->mesh.nsd]; jj++){
	    EEta[(i - 1) * vpoints[E->mesh.nsd] + jj] = E->viscosity.N0[l];
	  }
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

	  switch(E->viscosity.RHEOL){ 
	  case 0:
	    /* eta = eta0 exp(E * (1-T)) */
	    for(i = 1; i <= nel; i++)
	      {
		l = E->mat[i] - 1; /* moved -1 up here */
		e = (i - 1) % E->lmesh.elz + 1;
		
		tempa = E->viscosity.N0[l];
		
		for(kk = 1; kk <= ends; kk++)
		  TT[kk] = E->T[E->ien[i].node[kk]];
		
		for(jj = 1; jj <= vpts; jj++)
		  {
		    temp = 1.0e-32;
		    for(kk = 1; kk <= ends; kk++)
		      {
			temp += max(zero, TT[kk]) * E->N.vpt[GNVINDEX(kk, jj)];;
		      }
		    EEta[(i - 1) * vpts + jj] = tempa * exp(E->viscosity.E[l] * (one - temp ));
		  }
	      }
	    break;
	  case 1:
	    /* 
	       eta = eta0 * exp(E/(T+T0))
	    */
	    for(i = 1; i <= nel; i++)
	      {
		l = E->mat[i] - 1 ;
		tempa = E->viscosity.N0[l];
		
		for(kk = 1; kk <= ends; kk++)
		  TT[kk] = E->T[E->ien[i].node[kk]];
		
		for(jj = 1; jj <= vpts; jj++)
		  {
		    temp = 1.0e-32;
		    for(kk = 1; kk <= ends; kk++)
		      {
			temp += max(zero, TT[kk]) * E->N.vpt[GNVINDEX(kk, jj)];;
		      }
		    EEta[(i - 1) * vpts + jj] = tempa * exp((E->viscosity.E[l]) / (temp + E->viscosity.T[l]));
		  }
	      }
	  case 2:
	    /* 
	       eta = eta0 * exp(E + (1-z)*Z0/(T+T0))
	    */
	    for(i = 1; i <= nel; i++)
	      {
		l = E->mat[i] - 1 ;

		tempa = E->viscosity.N0[l];
		
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
		    EEta[(i - 1) * vpts + jj] = tempa * exp((E->viscosity.E[l] + (1 - zz) * E->viscosity.Z[l]) /(temp + E->viscosity.T[l]) );
		  }
	      }
	    break;
	  case 3:
	    /* 
	       eta = eta0 exp(E * (T0 - T)) 

	       (this is like rheol=0 with T0=1, but I left rheol==0
	       unchanged for backward compatibility

	    */
	    for(i = 1; i <= nel; i++)
	      {
		l = E->mat[i] - 1;
		e = (i - 1) % E->lmesh.elz + 1;
		
		tempa = E->viscosity.N0[l];
		
		for(kk = 1; kk <= ends; kk++)
		  TT[kk] = E->T[E->ien[i].node[kk]];
		
		for(jj = 1; jj <= vpts; jj++)
		  {
		    temp = 1.0e-32;
		    for(kk = 1; kk <= ends; kk++)
		      {
			temp += max(zero, TT[kk]) * E->N.vpt[GNVINDEX(kk, jj)];;
		      }
		    EEta[(i - 1) * vpts + jj] = tempa * exp(E->viscosity.E[l] * (E->viscosity.T[l]  - temp ));
		  }
	      }
	    break;
	  case 4:
	    /* 
	       eta = eta0 * exp(E(T_c - T) + (1-z) * Z ) 
	    */
	    for(i = 1; i <= nel; i++) /* loop over all elements */
	      {
		l = E->mat[i] - 1; /* material element, determines [,,,,]  */
		tempa = E->viscosity.N0[l];	/* prefactor */
		
		for(kk = 1; kk <= ends; kk++) /* temperature at nodes of element */
		  TT[kk] = E->T[E->ien[i].node[kk]];
		
		for(jj = 1; jj <= vpts; jj++) /* loop through integration points of element */
		  {
		    zz = 0;
		    temp = 1.0e-32;
		    for(kk = 1; kk <= ends; kk++) /* loop through points in element */
		      {
			temp += max(zero, TT[kk]) * E->N.vpt[GNVINDEX(kk, jj)];;
			zz += Xtmp[3][E->ien[i].node[kk]] * E->N.vpt[GNVINDEX(kk, jj)]; /* depth */
		      }
		    /* now we have the temperature at this integration point */
		    EEta[(i - 1) * vpts + jj] = 
		      tempa * exp(E->viscosity.E[l] * (E->viscosity.T[l] - temp) + (1 - zz) * E->viscosity.Z[l] );
		  }
	      }
	    break;
	  case 10:
	    /* 
	       eta = eta0 * exp(E/(T+T0) + Z) for testing purposes
	    */
	    for(i = 1; i <= nel; i++)
	      {
		l = E->mat[i] - 1 ;
		tempa = E->viscosity.N0[l];
		
		for(kk = 1; kk <= ends; kk++)
		  TT[kk] = E->T[E->ien[i].node[kk]];
		
		for(jj = 1; jj <= vpts; jj++)
		  {
		    temp = 1.0e-32;
		    for(kk = 1; kk <= ends; kk++)
		      {
			temp += max(zero, TT[kk]) * E->N.vpt[GNVINDEX(kk, jj)];;
		      }
		    EEta[(i - 1) * vpts + jj] = tempa * exp((E->viscosity.E[l]) / (temp + E->viscosity.T[l]) + E->viscosity.Z[l]) ;
		  }
	      }
	    break;
	  case 11:
            /* 
               eta = eta0 * exp(E/(T+T0) - E/(0.5+T0))
            */
            for(i = 1; i <= nel; i++)
              {
                l = E->mat[i] - 1 ;
                tempa = E->viscosity.N0[l];

                for(kk = 1; kk <= ends; kk++)
                  TT[kk] = E->T[E->ien[i].node[kk]];

                for(jj = 1; jj <= vpts; jj++)
                  {
                    temp = 1.0e-32;
                    for(kk = 1; kk <= ends; kk++)
                      {
                        temp += max(zero, TT[kk]) * E->N.vpt[GNVINDEX(kk, jj)];;
                      }
                    EEta[(i - 1) * vpts + jj] = tempa * exp((E->viscosity.E[l]) / (temp + E->viscosity.T[l]) - (E->viscosity.E[l]) / (0.5 + E->viscosity.T[l]));
                  }
              }
            break;
	  default:
	    myerror("RHEOL option undefined",E);
	    break;
	  } /* end swith */

	  visits++;
	  
	}



	/* 	fprintf(E->fp,"aaa\n"); */
	/* 	for(i=1;i<=nel;i++) */
	/* 		fprintf(E->fp,"%d %d %g\n",i,E->mat[i],EEta[(i-1)*vpts+1]); */
	

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


void strain_rate_2_inv(struct All_variables *E, float *EEDOT, 
		       int SQRT)
{
  double edot[4][4], dudx[4][4], rtf[4][9];
  float VV[4][9], Vxyz[9][9];
  
  //int e, i, j, p, q, n, nel, k;
  int e, i, j, p, q, n, nel;
  
  const int dims = E->mesh.nsd;
  const int ends = enodes[dims];
  const int lev = E->mesh.levmax;
  const int ppts = ppoints[dims];
  
  nel = E->lmesh.nel;
  
  
  if(E->control.Rsphere){
    
    for(e = 1; e <= nel; e++){

      get_rtf(E, e, 1, rtf, lev); /* pressure */
      
      for(i = 1; i <= ends; i++){
	n = E->ien[e].node[i];
	VV[1][i] = E->V[1][n];
	VV[2][i] = E->V[2][n];
	VV[3][i] = E->V[3][n];
      }
      
      for(j = 1; j <= ppts; j++){
	Vxyz[1][j] = 0.0;
	Vxyz[2][j] = 0.0;
	Vxyz[3][j] = 0.0;
	Vxyz[4][j] = 0.0;
	Vxyz[5][j] = 0.0;
	Vxyz[6][j] = 0.0;
      }

      for(j = 1; j <= ppts; j++){ /* only makes "sense" for ppts = 1 */
	for(i = 1; i <= ends; i++){
	  Vxyz[1][j] += (VV[1][i] * E->gNX[e].ppt[GNPXINDEX(0, i, j)] + VV[3][i] * E->N.ppt[GNPINDEX(i, j)]) * rtf[3][j]; /* tt */
	  Vxyz[2][j] += ((VV[2][i] * E->gNX[e].ppt[GNPXINDEX(1, i, j)] + VV[1][i] * E->N.ppt[GNPINDEX(i, j)] * cos(rtf[1][j])) / sin(rtf[1][j]) + VV[3][i] * E->N.ppt[GNPINDEX(i, j)]) * rtf[3][j]; /* pp */
	  Vxyz[3][j] += VV[3][i] * E->gNX[e].ppt[GNPXINDEX(2, i, j)]; /* rr */
	  
	  Vxyz[4][j] += ((VV[1][i] * E->gNX[e].ppt[GNPXINDEX(1, i, j)] - VV[2][i] * E->N.ppt[GNPINDEX(i, j)] * cos(rtf[1][j])) / sin(rtf[1][j]) + VV[2][i] * E->gNX[e].ppt[GNPXINDEX(0, i, j)]) * rtf[3][j]; /* tp */
	  Vxyz[5][j] += VV[1][i] * E->gNX[e].ppt[GNPXINDEX(2, i, j)] + rtf[3][j] * (VV[3][i] * E->gNX[e].ppt[GNPXINDEX(0, i, j)] - VV[1][i] * E->N.ppt[GNPINDEX(i, j)]); /* tr */
	  Vxyz[6][j] += VV[2][i] * E->gNX[e].ppt[GNPXINDEX(2, i, j)] + rtf[3][j] * (VV[3][i] * E->gNX[e].ppt[GNPXINDEX(1, i, j)] / sin(rtf[1][j]) - VV[2][i] * E->N.ppt[GNPINDEX(i, j)]); /* pr */
	}
	edot[1][1] = 2.0 * Vxyz[1][j]; /* this should be a summation */
	edot[2][2] = 2.0 * Vxyz[2][j];
	edot[3][3] = 2.0 * Vxyz[3][j];
	edot[1][2] = Vxyz[4][j];
	edot[1][3] = Vxyz[5][j];
	edot[2][3] = Vxyz[6][j];
      }
      /* /\* eps: rr, rt, rp, tt, tp, pp *\/ */
      /* if(e < 5)fprintf(stderr,"1: %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n",edot[3][3]/2,edot[1][3]/2,edot[2][3]/2,edot[1][1]/2,edot[1][2]/2,edot[2][2]/2); */
      EEDOT[e] = edot[1][1] * edot[1][1] + edot[1][2] * edot[1][2] * 2.0 + edot[2][2] * edot[2][2] + edot[2][3] * edot[2][3] * 2.0 + edot[3][3] * edot[3][3] + edot[1][3] * edot[1][3] * 2.0;
      
    }

  }else if(E->control.CART3D){
    for(e = 1; e <= nel; e++){
      for(i = 1; i <= ends; i++){
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

      /* if(e < 5)fprintf(stderr,"1: %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n",edot[1][1]/2,edot[1][2]/2,edot[1][3]/2,edot[2][2]/2,edot[2][3]/2,edot[3][3]/2); */

      
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

/* compute second invariant from a strain-rate tensor in 0,...2 format

 */
double second_invariant_from_3x3(double e[3][3])
{
  return(sqrt(0.5*
	      (e[0][0] * e[0][0] + 
	       e[0][1] * e[0][1] * 2.0 + 
	       e[1][1] * e[1][1] + 
	       e[1][2] * e[1][2] * 2.0 + 
	       e[2][2] * e[2][2] + 
	       e[0][2] * e[0][2] * 2.0)));
}

/* 

   calculate velocity gradient matrix for all elements

   l is stored as (e-1)*9 
   
*/
void calc_vgm_matrix(struct All_variables *E, double *l,double *evel)
{
  double rtf[4][9];
  double VV[4][9],vgm[3][3];
  int e,j,i,n,p,q,eoff;
  double d[6];
  static struct CC Cc;
  static struct CCX Ccx;
  const int dims = E->mesh.nsd;
  const int ppts = ppoints[dims];
  const int ends = enodes[dims];
  const int lev = E->mesh.levmax;
  const int nel = E->lmesh.nel;

  for(eoff=0,e=1; e <= nel; e++, eoff += 9) {
    if(E->control.Rsphere){	/* need rtf for spherical */
      get_rtf(E, e, 1, rtf, lev); /* pressure points */
      if((e-1)%E->lmesh.elz==0)
	construct_c3x3matrix_el(E,e,&Cc,&Ccx,lev,1);
    }
    for(i = 1; i <= ends; i++){	/* velocity at element nodes */
      n = E->ien[e].node[i];
      VV[1][i] = E->V[1][n];
      VV[2][i] = E->V[2][n];
      VV[3][i] = E->V[3][n];
    }
    get_vgm_p(VV,&(E->N),&(E->GNX[lev][e]),&Cc, &Ccx,rtf,
	      E->mesh.nsd,ppts,ends,(E->control.Rsphere),vgm,evel);
    get_9vec_from_3x3((l+eoff),vgm);
  }
}

/* 

   get velocity gradient matrix at element, and also compute the
   average velocity in this element
   

*/

void get_vgm_p(double VV[4][9],struct Shape_function *N,
	       struct Shape_function_dx *GNx,
	       struct CC *cc, struct CCX *ccx, double rtf[4][9],
	       int dims,int ppts, int ends, int spherical,
	       double l[3][3], double v[3])
{

  int i,k,j,a;
  double ra[9], si[9], ct[9];
  const double one = 1.0;
  const double two = 2.0;
  double vgm[3][3];
  double shp, cc1, cc2, cc3,d_t,d_r,d_p,up,ur,ut;
  /* init L matrix */
  for(i=0;i < 3;i++){
    v[i] = 0.0;
    for(j=0;j < 3; j++)
      l[i][j] = 0.0;
  }
  /* mean velocity at pressure integration point */
  for(a=1;a <= ends;a++){
    v[0] += N->ppt[GNPINDEX(a, 1)] * VV[1][a];
    v[1] += N->ppt[GNPINDEX(a, 1)] * VV[2][a];
    v[2] += N->ppt[GNPINDEX(a, 1)] * VV[3][a];
  }
  if(spherical){
    for(k = 1; k <= ppts; k++){
      ra[k] = rtf[3][k];	      /* 1/r */
      si[k] = one / sin(rtf[1][k]); /* 1/sin(t) */
      ct[k] = cos(rtf[1][k]) * si[k]; /* cot(t) */
    }
    for(a = 1; a <= ends; a++){
      for(k = 1; k <= ppts; k++){
	d_t = GNx->ppt[GNPXINDEX(0, a, k)]; /* d_t */
	d_p = GNx->ppt[GNPXINDEX(1, a, k)]; /* d_p */
	d_r = GNx->ppt[GNPXINDEX(2, a, k)]; /* d_r */
	shp = N->ppt[GNPINDEX(a, k)];
	for(i = 1; i <= dims; i++){
	  ut = cc->ppt[BPINDEX(1, i, a, k)]; /* ut */
	  up = cc->ppt[BPINDEX(2, i, a, k)]; /* up */
	  ur = cc->ppt[BPINDEX(3, i, a, k)]; /* ur */
	  
	  /* velocity gradient matrix is transpose of grad v, using Citcom sort t, p, r
	
	     | d_t(vt) d_p(vt) d_r(vt) |
	     | d_t(vp) d_p(vp) d_r(vp) |
	     | d_t(vr) d_p(vr) d_r(vr) |

	  */

	  /* d_t vt = 1/r (d_t vt + vr) */
	  vgm[0][0] =  ((d_t * ut + shp * ccx->ppt[BPXINDEX(1, i, 1, a, k)]) + 
			shp * ur) * ra[k];
	  /* d_p vt = 1/r (1/sin(t) d_p vt -vp/tan(t)) */
	  vgm[0][1] =  ((d_p * ut + shp * ccx->ppt[BPXINDEX(1, i, 2, a, k)]) * si[k] - 
			shp * up * ct[k]) * ra[k];
	  /* d_r vt = d_r v_t */
	  vgm[0][2] = d_r * ut;
	  /* d_t vp = 1/r d_t v_p*/
	  vgm[1][0] = (d_t * up + shp * ccx->ppt[BPXINDEX(2, i, 1, a, k)]) * ra[k];
	  /* d_p vp = 1/r((d_p vp)/sin(t) + vt/tan(t) + vr) */
	  vgm[1][1] = ((d_p * up + shp * ccx->ppt[BPXINDEX(2, i, 2, a, k)]) * si[k] + 
		       shp * ut * ct[k] + shp * ur) * ra[k];
	  /* d_r vp = d_r v_p */
	  vgm[1][2] =  d_r * up;
	  /* d_t vr = 1/r(d_t vr - vt) */
	  vgm[2][0] = ((d_t * ur + shp * ccx->ppt[BPXINDEX(3, i, 1, a, k)]) -
		       shp * ut) * ra[k];
	  /* d_p vr =  1/r(1/sin(t) d_p vr - vp) */
	  vgm[2][1] = (( d_p * ur + shp * ccx->ppt[BPXINDEX(3,i, 2,a,k)] ) * si[k] -
		       shp * up ) * ra[k];
	  /* d_r vr = d_r vr */
	  vgm[2][2] = d_r * ur;


	  l[0][0] += vgm[0][0] * VV[i][a];
	  l[0][1] += vgm[0][1] * VV[i][a];
	  l[0][2] += vgm[0][2] * VV[1][a];
	  
	  l[1][0] += vgm[1][0] * VV[i][a];
	  l[1][1] += vgm[1][1] * VV[i][a];
	  l[1][2] += vgm[1][2] * VV[i][a];
	  
	  l[2][0] += vgm[2][0] * VV[i][a];
	  l[2][1] += vgm[2][1] * VV[i][a];
	  l[2][2] += vgm[2][2] * VV[i][a];
	  
	}
      }
    }
  }else{		
    /* cartesian */
    for(k = 1; k <= ppts; k++){
      for(a = 1; a <= ends; a++){
	/* velocity gradient matrix is transpose of grad v
	
	     | d_x(vx) d_y(vx) d_z(vx) |
	     | d_x(vy) d_y(vy) d_z(vy) |
	     | d_x(vz) d_y(vz) d_z(vz) |
	*/
	l[0][0] += GNx->ppt[GNPXINDEX(0, a, k)] * VV[1][a]; /* other contributions are zero */
	l[0][1] += GNx->ppt[GNPXINDEX(1, a, k)] * VV[1][a];
	l[0][2] += GNx->ppt[GNPXINDEX(2, a, k)] * VV[1][a];

	l[1][0] += GNx->ppt[GNPXINDEX(0, a, k)] * VV[2][a];
	l[1][1] += GNx->ppt[GNPXINDEX(1, a, k)] * VV[2][a];
	l[1][2] += GNx->ppt[GNPXINDEX(2, a, k)] * VV[2][a];

	l[2][0] += GNx->ppt[GNPXINDEX(0, a, k)] * VV[3][a];
	l[2][1] += GNx->ppt[GNPXINDEX(1, a, k)] * VV[3][a];
	l[2][2] += GNx->ppt[GNPXINDEX(2, a, k)] * VV[3][a];

      }
    }
  }
  if(ppts != 1){
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	l[i][j] /= (float)ppts;
  }

}




/* 

   given a 3x3 velocity gradient matrix l, compute a d[3][3]
   strain-rate matrix

*/

void calc_strain_from_vgm(double l[3][3], double d[3][3])
{
  int i,j;
  for(i=0;i < 3;i++)
    for(j=0;j < 3;j++)
      d[i][j] = 0.5 * (l[i][j] + l[j][i]);
}
void calc_strain_from_vgm9(double *l9, double d[3][3])
{
  double l[3][3];
  get_3x3_from_9vec(l, l9);
  calc_strain_from_vgm(l, d);
}
/* 

   given a 3x3 velocity gradient matrix l, compute a rotation matrix

*/

void calc_rot_from_vgm(double l[3][3], double r[3][3])
{
  int i,j;
  for(i=0;i < 3;i++)
    for(j=0;j < 3;j++)
      r[i][j] = 0.5 * (l[i][j] - l[j][i]);
}

/* 

   different version of strain-rate computation, obtain strain-rate
   matrix at all elements, storing it in edot as (e-1)*6 upper
   triangle format


*/
void calc_strain_rate_matrix(struct All_variables *E, 
			     double *edot)
{
  double rtf[4][9];
  double VV[4][9], Vxyz[7][9];
  int e,j,i,n,p,q,eoff;
  static struct CC Cc;
  static struct CCX Ccx;
  double ba[9][4][9][7];
  const int dims = E->mesh.nsd;
  const int ppts = ppoints[dims];
  const int ends = enodes[dims];
  const int lev = E->mesh.levmax;
  const int nel = E->lmesh.nel;
  //fprintf(stderr,"\n");
  for(eoff=0,e=1; e <= nel; e++,eoff+=6) {
    if(E->control.Rsphere){	/* need rtf for spherical */
      get_rtf(E, e, 1, rtf, lev); /* pressure */
      if((e-1)%E->lmesh.elz==0)
	construct_c3x3matrix_el(E,e,&Cc,&Ccx,lev,1);
    }
    for(i = 1; i <= ends; i++){	/* velocity at element nodes */
      n = E->ien[e].node[i];
      VV[1][i] = E->V[1][n];
      VV[2][i] = E->V[2][n];
      VV[3][i] = E->V[3][n];
    }
    for(j = 1; j <= ppts; j++){
      Vxyz[1][j] = 0.0;
      Vxyz[2][j] = 0.0;
      Vxyz[3][j] = 0.0;
      Vxyz[4][j] = 0.0;
      Vxyz[5][j] = 0.0;
      Vxyz[6][j] = 0.0;
    }
    get_ba_p(&(E->N),&(E->GNX[lev][e]),&Cc, &Ccx,rtf,
	     E->mesh.nsd,ppts,ends,(E->control.Rsphere),ba);
    for(j=1;j <= ppts;j++)
      for(p=1;p <= 6;p++)
	for(i=1;i <= ends;i++)
	  for(q=1;q <= dims;q++) {
	    Vxyz[p][j] += ba[i][q][j][p] * VV[q][i];
	  }
    edot[eoff+0] = edot[eoff+1] = edot[eoff+2] =
      edot[eoff+3] = edot[eoff+4] = edot[eoff+5] = 0.0;

    for(j=1; j <= ppts; j++) {
      edot[eoff+0] += Vxyz[1][j];   /* e_xx = e_tt */
      edot[eoff+1] += Vxyz[4][j]/2; /* e_xy = e_tp */
      edot[eoff+2] += Vxyz[5][j]/2; /* e_xz = e_tr */
      edot[eoff+3] += Vxyz[2][j];   /* e_yy = e_pp */
      edot[eoff+4] += Vxyz[6][j]/2; /* e_yz = e_pr */
      edot[eoff+5] += Vxyz[3][j];   /* e_zz = e_rr */
    }
    if(ppts != 1){
      for(i=0;i<6;i++)
	edot[eoff+i] /= (float)ppts;
    }

    /* if(e < 5){ */
    /*   if(E->control.Rsphere){ */
    /* 	/\* rr rt rp tt tp pp *\/ */
    /* 	fprintf(stderr,"2: %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e \n",edot[eoff+5],edot[eoff+2],edot[eoff+4],edot[eoff+0],edot[eoff+1],edot[eoff+3]); */
    /*   }else{ */
    /* 	fprintf(stderr,"2: %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e \n",edot[eoff+0],edot[eoff+1],edot[eoff+2],edot[eoff+3],edot[eoff+4],edot[eoff+5]); */
    /*   } */
    /* } */
  }
}


int layers(struct All_variables *E, float x3)
{
	int llayers = 0;
	int i,ncheck;
	ncheck = E->viscosity.num_mat-1;
	for(i=0;i < ncheck;i++)
	  if(x3 >= E->viscosity.zbase_layer[i]){
	    llayers = i+1;
	    break;
	  }
	if(!llayers)
	  llayers = E->viscosity.num_mat;
	    
  
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

/* 

limit the viscosity by a plastic-type yielding with depth dependent
yield stress (aka byerlee law)



psrw = true

		
a strain-rate weakening rheology applies
based on a steady state stress-strain
relationship following, e.g., Tackley 1998

\eta_ett = (\sigma_y^2 \eta_0)/(\sigma_y^2 + \eta_0^2 (2 \eps_II^2))

where sigma_y is defined as above

where the 2\eps_II arises because our \eps_II has the 1/2 factor in it


*/
//#define DEBUG
static void visc_from_B(struct All_variables *E, float *Eta, float *EEta, int propogate)
{
  static int visited = 0;
  float scale,stress_magnitude,depth,exponent1,eta_old,eta_old2,eta_new;
  float *eedot;
  float zzz,zz[9];
  float tau,tau2,ettby,ettnew;
  int m,l,z,jj,kk,i;
  static float ndz_to_m;
#ifdef DEBUG
  FILE *out;
#endif
  const int vpts = vpoints[E->mesh.nsd];
  const int nel = E->lmesh.nel;
  const int ends = enodes[E->mesh.nsd];
  
  eedot = (float *) safe_malloc((2+nel)*sizeof(float));
 
#ifdef DEBUG
  out=fopen("tmp.visc","w");
#endif
  if (!visited)   {
    /* 
       scaling from nod-dim radius (0...1) to m 

       (only used for dimensional version)

    */
    ndz_to_m = E->monitor.length_scale;

    /*  */
    if(E->parallel.me==0){	/* control output */
      for(l=1;l <= E->viscosity.num_mat;l++) {
	fprintf(stderr,"Plasticity: %d/%d: a=%g b=%g p=%g\n",
		l,E->viscosity.num_mat,
		E->viscosity.abyerlee[l-1],
		E->viscosity.bbyerlee[l-1],
		E->viscosity.lbyerlee[l-1]);
      }
      fprintf(stderr,"\tdim: %i trans: %i offset: %g\n",
	      E->viscosity.plasticity_dimensional,
	      E->viscosity.plasticity_trans,
	      E->viscosity.plasticity_viscosity_offset);
      fprintf(stderr,"\tpsrw: %i\n",E->viscosity.psrw);
    }
    /* 
       get strain rates for all elements 
    */
    for(i=1;i<=nel;i++)
      eedot[i] = 1.0; 

  }else{
    if(E->viscosity.psrw)
      strain_rate_2_inv(E,eedot,0);
    else
      strain_rate_2_inv(E,eedot,1);
  }
  if(E->viscosity.psrw){
    /* strain-rate weakening */
    for(i=1;i <= nel;i++)   {	
      l = E->mat[i] - 1;	/* material of element */
      for(kk=1;kk <= ends;kk++){/* loop through integration points*/
	if(E->control.Rsphere)
	  zz[kk] = (1.0 - E->SX[3][E->ien[i].node[kk]]); 
	else
	  zz[kk] = (1.0 - E->X[3][E->ien[i].node[kk]]); 
	if(E->viscosity.plasticity_dimensional)
	  zz[kk] *= ndz_to_m;	/* scale to meters */
      }
      for(jj=1;jj <= vpts;jj++) { 
	zzz = 0.0;
	for(kk=1;kk <= ends;kk++)   {
	  zzz += zz[kk] * E->N.vpt[GNVINDEX(kk,jj)];
	}
	if(E->viscosity.plasticity_dimensional){
	  tau = (E->viscosity.abyerlee[l] * zzz + 
		 E->viscosity.bbyerlee[l]) * E->viscosity.lbyerlee[l];
	  tau /= E->monitor.tau_scale;
	}else{
	  tau = E->viscosity.abyerlee[l] * zzz + E->viscosity.bbyerlee[l];
	  tau = min(tau,  E->viscosity.lbyerlee[l]);
	}
	if((visited > 1) && (tau < 1e15)){
	  tau2 = tau * tau;
	  eta_old = EEta[ (i-1)*vpts + jj ] ;
	  eta_old2 = eta_old * eta_old;
	  eta_new = (tau2 * eta_old)/(tau2 + 2.0 * eta_old2 * eedot[i]);
	  EEta[ (i-1)*vpts + jj ] = ettnew;
	  //if(E->parallel.me==0)fprintf(stderr,"tau: %11.4e eII: %11.4e eta_old: %11.4e eta_new: %11.4e\n",tau, eedot[i],eta_old,eta_new);
	}
      }
    }
  }else{

    /* regular plasticity */

    
    for(i=1;i <= nel;i++)   {	
      /* 
	 loop through all elements 
      */
      l = E->mat[i] - 1;	/* material of element */
      for(kk=1;kk <= ends;kk++){/* loop through integration points*/
	/* depth in meters */
	if(E->control.Rsphere)
	  zz[kk] = (1.0 - E->SX[3][E->ien[i].node[kk]]); 
	else
	  zz[kk] = (1.0 - E->X[3][E->ien[i].node[kk]]); 
	if(E->viscosity.plasticity_dimensional)
	  zz[kk] *= ndz_to_m;	/* scale to meters */
      }
      for(jj=1;jj <= vpts;jj++) { 
	/* loop over nodes in element */
	zzz = 0.0;
	for(kk=1;kk <= ends;kk++)   {
	  /* 
	     depth  [m] 
	  */
	  zzz += zz[kk] * E->N.vpt[GNVINDEX(kk,jj)];
	}
	if(E->viscosity.plasticity_dimensional){
	  /* byerlee type */
	  
	  /* 
	     yield stress in [Pa] 
	  */
	  tau=(E->viscosity.abyerlee[l] * zzz + E->viscosity.bbyerlee[l]) * E->viscosity.lbyerlee[l];
	  /* 
	     scaled stress 
	  */
	  tau /= E->monitor.tau_scale;
	}else{
	  
	  tau = E->viscosity.abyerlee[l] * zzz + E->viscosity.bbyerlee[l];
	  
	  tau = min(tau,  E->viscosity.lbyerlee[l]);
	}
	/* 
	   
	`byerlee viscosity' : tau = 2 eps eta, this is non-dim
	plus some offset as in Stein et al. 
	
	*/
	ettby = tau/(2.0 * (eedot[i]+1e-7)) + E->viscosity.plasticity_viscosity_offset;
	/* 
	   
	decide on the plasticity transition
	
	
	*/
	if(E->viscosity.plasticity_trans){
	  /* 
	     eta = 1./(1./eta(k)+1./ettby)  
	  */
	  
	  ettnew = 1.0/(1.0/EEta[ (i-1)*vpts + jj ] + 1.0/ettby);
	  //fprintf(stderr,"a: %g %g %g\n",EEta[ (i-1)*vpts + jj ],ettby,ettnew);
	}else{
	  /* 
	     min(\eta_p, \eta_visc )
	  */
	  ettnew = min(EEta[ (i-1)*vpts + jj ],ettby);
	  //fprintf(stderr,"m: %g %g %g\n",EEta[ (i-1)*vpts + jj ],ettby,ettnew);
	}
#ifdef DEBUG
	/* output format 
	   
	z[km] tau[MPa] eps[s^-1] eta_b[Pas] eta_T[Pas] eta_c[Pas]
	
	*/
	if(visited)
	  fprintf(out,"%10.2f %17.4e %17.4e %17.4e %17.4e %17.4e\n", 
		  zzz/1e3,tau*E->monitor.tau_scale/1e6, 
		  eedot[i]/E->monitor.time_scale, 
		  ettby*E->data.ref_viscosity, 
		  EEta[ (i-1)*vpts + jj ]*E->data.ref_viscosity, 
		  ettnew*E->data.ref_viscosity); 
#endif
	//      if(visited)
	//	fprintf(stderr,"%11.4e %11.4e %11.4e %11.4e\n",ettnew,EEta[ (i-1)*vpts + jj ] ,ettby,eedot[i]);
	EEta[ (i-1)*vpts + jj ] = ettnew;
      }
    } /* end regular plasticity */
  }
#ifdef DEBUG
  fclose(out);
#endif
  visited = 1;
  free ((void *)eedot);
  return;  
}
/* 

multiply with composition factor

*/
static void visc_from_C(struct All_variables *E, float *Eta, float *EEta, int propogate)
{
  float comp,comp_fac,CC[9],tcomp;
  double vmean,cc_loc;
  int m,l,z,jj,kk,i;
  static int visited = 0;
  static double logv[2];
  const int vpts = vpoints[E->mesh.nsd];
  const int nel = E->lmesh.nel;
  const int ends = enodes[E->mesh.nsd];
  if(!visited){
    /* log of the material viscosities */
    logv[0] = log(E->viscosity.pre_comp[0]);
    logv[1] = log(E->viscosity.pre_comp[1]);
  }
  for(i = 1; i <= nel; i++)
    {
      /* determine composition of each of the nodes of the element */
      for(kk = 1; kk <= ends; kk++){
	CC[kk] = E->C[E->ien[i].node[kk]];
	if(E->control.check_c_irange){
	  if(CC[kk] < 0)CC[kk]=0.0;
	  if(CC[kk] > 1)CC[kk]=1.0;
	}
      }
      for(jj = 1; jj <= vpts; jj++)
	{
	  /* compute mean composition  */
	  cc_loc = 0.0;
	  for(kk = 1; kk <= ends; kk++){/* the vpt takes care of averaging */
	    cc_loc += CC[kk] * E->N.vpt[GNVINDEX(kk, jj)];
	  }
	  /* geometric mean of viscosity */
	  vmean = exp(cc_loc  * logv[1] + (1.0-cc_loc) * logv[0]);
	  EEta[ (i-1)*vpts + jj ] *= vmean;
	} /* end jj loop */
    }
  visited++;
}

