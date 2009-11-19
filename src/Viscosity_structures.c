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

	for(i = 0; i < 40; i++)
	{
		E->viscosity.N0[i] = 1.0;
		E->viscosity.T[i] = 0.0;
		E->viscosity.Z[i] = 0.0;
		E->viscosity.E[i] = 0.0;
		E->viscosity.T0[i] = 0.0;

		
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
		
	} /* end layer loop */

	/* read in information */
	input_int("rheol", &(E->viscosity.RHEOL), "essential",m);
	input_int("num_mat", &(E->viscosity.num_mat), "1",m);


	input_float_vector("viscT", E->viscosity.num_mat, (E->viscosity.T),m);	/* redundant */
	input_float_vector("viscT1", E->viscosity.num_mat, (E->viscosity.T),m);
	input_float_vector("viscZ", E->viscosity.num_mat, (E->viscosity.Z),m);
	input_float_vector("viscE", E->viscosity.num_mat, (E->viscosity.E),m);
	input_float_vector("viscT0", E->viscosity.num_mat, (E->viscosity.T0),m);
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


	
	if(E->viscosity.TDEPV)
		visc_from_T(E, visc, evisc, propogate);
	else
		visc_from_mat(E, visc, evisc);

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

	/* this is the old logic, llayer = 4 would never get assigned */
	/* 	if(x3 >= E->viscosity.zlith) */
	/* 		llayers = 1; */
	/* 	else if((x3 < E->viscosity.zlith) && (x3 >= E->viscosity.zlm)) */
	/* 		llayers = 2; */
	/* 	else if(x3 < E->viscosity.zlm) */
	/* 		llayers = 3; */
	
	if(x3 >= E->viscosity.zlith) /* above 410 */
	  llayers = 1;
	else if(x3 >= E->viscosity.z410) /* above 410 */
	  llayers = 2;
	else if(x3 >= E->viscosity.zlm) /* above 660 */
	  llayers = 3;
	else /* lower mantle */
	  llayers = 4;
  
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
	  //if(E->parallel.me==0)fprintf(stderr,"tau: %11g eII: %11g eta_old: %11g eta_new: %11g\n",tau, eedot[i],eta_old,eta_new);
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
	//	fprintf(stderr,"%11g %11g %11g %11g\n",ettnew,EEta[ (i-1)*vpts + jj ] ,ettby,eedot[i]);
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
	if(CC[kk] < 0)CC[kk]=0.0;
	if(CC[kk] > 1)CC[kk]=1.0;
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
