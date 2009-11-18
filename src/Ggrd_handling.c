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

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include <stdlib.h>				/* for "system" command */
#include <string.h>
#include <mpi.h>
#include "hc.h"
#include "function_defs.h"

/* ===============================
   Initialization of fields .....
   =============================== */

void convection_initial_temperature_ggrd(struct All_variables *E)
{
  int ll, mm, i, j, k, p, node, ii;
  double temp, temp1, temp2, temp3, base, radius, radius2;
  unsigned int type;

  FILE *fp;
  double  x1, y1, z1, con;

  int noz2, nfz, in1, in2, in3, instance, nox, noy, noz,nxs,setflag;
  char input_s[200], output_file[255];
  float weight, para1, plate_velocity, delta_temp, age;

  char *char_dummy="";
  
  /* twb additions */
  double rho_prem;
  char pfile[1000];
  double t1,f1,r1,tgrad,tbot,tadd,tz,tmean;

  /* for dealing with several processors */
  MPI_Status mpi_stat;
  int mpi_rc, mpi_tag=1;  
  int mpi_inmsg, mpi_success_message = 1;
	

  const int dims = E->mesh.nsd;

  noy = E->lmesh.noy;
  noz = E->lmesh.noz;
  nox = E->lmesh.nox;


  setflag = 0;


  for(i=1;i<=E->lmesh.nno;i++) {
    /* init as zeros */
    E->T[i] = E->C[i] = 0.0;
  }
  if(!E->control.restart)
    {				/* 

				regular init

				*/

      if(E->control.ggrd.use_temp){ /* read T init from grid */
	if(E->parallel.me==0)  
	  fprintf(stderr,"convection_initial_temperature: using GMT grd files for temperatures\n");
	/* 
	       
	read in tempeatures/density from GMT grd files
	    
	    
	*/
	if(E->parallel.me > 0){
	  /* 
	     wait for the previous processor 
	  */
	  mpi_rc = MPI_Recv(&mpi_inmsg, 1, MPI_INT, (E->parallel.me-1), mpi_tag, 
			    MPI_COMM_WORLD, &mpi_stat);
	}
	if(E->parallel.me == 0)
	  fprintf(stderr,"Init: initializing PREM and ggrd files\n");
	if(E->control.ggrd.temp.scale_with_prem){
	  if(!E->control.Rsphere)
	    myerror("cannot use PREM with Cartesian setting",E);
	  /* initialize PREM */
	  if(prem_read_model(PREM_MODEL_FILE,&E->control.ggrd.temp.prem, 
			     (E->parallel.me==0)))
	    myerror("prem initialization",E);
	}
	/* 
	   initialize the GMT grid files 
	*/
	if(E->control.slab_slice){ /* only slice of x - depth */

	  if(ggrd_grdtrack_init_general(FALSE,E->control.ggrd.temp.gfile, 
					char_dummy,"",E->control.ggrd.temp.d,
					(E->parallel.me==0),
					FALSE))
	    myerror("grdtrack x-z init error",E);
	}else{			/* 3-D */
	  if(ggrd_grdtrack_init_general(TRUE,E->control.ggrd.temp.gfile, 
					E->control.ggrd.temp.dfile,"",
					E->control.ggrd.temp.d,(E->parallel.me==0),
					FALSE))
	    myerror("grdtrack 3-D init error",E);
	}


	
	if(E->parallel.me <  E->parallel.nproc-1){
	  /* 
	     tell the next processor to go ahead 
	  */
	  mpi_rc = MPI_Send(&mpi_success_message, 1, MPI_INT, 
			    (E->parallel.me+1), mpi_tag, MPI_COMM_WORLD);
	}else{
	  fprintf(stderr,"convection_initial_temperature: last processor done with grd init\n");
	}
	/* 
	       
	end MPI synchornization  part
	    
	*/	
	/* 
	       
	interpolate densities to temperature given PREM variations
	    
	*/
	if(E->mesh.bottbc != 0){
	  /* bottom has specified temperature */
	  tbot =  E->control.TBCbotval;
	}else{
	  /* 
	     bottom has specified heat flux
	     start with unity bottom temperature
	  */
	  tbot = 1.0;
	}
      	tmean = (tbot + E->control.TBCtopval)/2.0 +  E->control.ggrd.temp.offset;
	if(E->parallel.me == 0)
	  fprintf(stderr,"tinit: top: %g bot: %g mean: %g scale: %g PREM %i\n",
		  E->control.TBCtopval,tbot,tmean,E->control.ggrd.temp.scale,
		  E->control.ggrd.temp.scale_with_prem);
	for(i=1;i<=noy;i++)  
	  for(j=1;j<=nox;j++) 
	    for(k=1;k<=noz;k++)  {
	      node=k+(j-1)*noz+(i-1)*nox*noz; /* offset */
	      if(E->control.slab_slice){
		if(E->control.Rsphere) {
		  if(E->SX[1][node] <= E->control.slab_theta_bound)
		    /* spherical interpolation */
		    ggrd_grdtrack_interpolate_xy((double)E->SX[2][node] * ONEEIGHTYOVERPI,
						 (double)E->SX[1][node],
						 E->control.ggrd.temp.d,&tadd,
						 FALSE);
		  else{
		    if(E->SX[3][node] == E->segment.zzlayer[E->segment.zlayers-1])
		      tadd = E->control.TBCtopval;
		    else
		      tadd = 1.0;
		  }
		}else{		/* cartesian interpolation */
		  if(E->X[2][node] <= E->control.slab_theta_bound)
		    ggrd_grdtrack_interpolate_xy((double)E->X[1][node],
						 (double)E->X[3][node],
						 E->control.ggrd.temp.d,&tadd,
						 FALSE);
		  else{
		    if(E->X[3][node] == E->segment.zzlayer[E->segment.zlayers-1])
		      tadd = E->control.TBCtopval;
		    else
		      tadd = 1.0;
		  }
		}
	      }else{
		/* 
		   
		3-D

		*/
		if(E->control.Rsphere) /* spherical interpolation */
		  ggrd_grdtrack_interpolate_rtp((double)E->SX[3][node],
						(double)E->SX[1][node],
						(double)E->SX[2][node],
						E->control.ggrd.temp.d,&tadd,
						FALSE,FALSE);
		else{		/* cartesian interpolation */
		  ggrd_grdtrack_interpolate_xyz((double)E->X[1][node],
						(double)E->X[2][node],
						(double)E->X[3][node],
						E->control.ggrd.temp.d,&tadd,
						FALSE);
		}
	      }
	      if(E->control.ggrd.temp.scale_with_prem){ /* only works for spherical! */
		/* 
		   get the PREM density at r for additional scaling  
		*/
		prem_get_rho(&rho_prem,(double)E->SX[3][node],&E->control.ggrd.temp.prem);
		if(rho_prem < 3200.0)
		  rho_prem = 3200.0; /* we don't want the viscosity of water */
		/* 
		   assign temperature 
		*/
		E->T[node] = tmean + tadd * E->control.ggrd.temp.scale * 
		  rho_prem / E->data.density;
	      }else{
		/* no PREM scaling */
		E->T[node] = tmean + tadd * E->control.ggrd.temp.scale;
		
		//fprintf(stderr,"z: %11g mean: %11g tadd: %11g scale: %11g T: %11g\n", E->X[3][node],tmean, tadd ,E->control.ggrd_tinit_scale,	E->T[node] );

	      }
	      /* 
		 if we're on the surface or bottom, reassign T and
		 temperature boundary condition if toptbc or bottbsc == 2
	      */
	      if(E->control.Rsphere){
		if((E->mesh.toptbc == 2) && (E->SX[3][node] == E->segment.zzlayer[E->segment.zlayers-1])){
		  E->TB[1][node] = E->TB[2][node] = E->TB[3][node] = E->T[node];
		}
		if((E->mesh.bottbc == 2) && (E->SX[3][node] == E->segment.zzlayer[0])){
		  E->TB[1][node] = E->TB[2][node] = E->TB[3][node] = E->T[node];
		}
	      }else{
		if((E->mesh.toptbc == 2) && (E->X[3][node] == E->segment.zzlayer[E->segment.zlayers-1]))
		  E->TB[1][node] = E->TB[2][node] = E->TB[3][node] = E->T[node];
		if((E->mesh.bottbc == 2) && (E->X[3][node] == E->segment.zzlayer[0]))
		  E->TB[1][node] = E->TB[2][node] = E->TB[3][node] = E->T[node];
	      }
	      /* 
		 boundary condition flags 
	      */
	      if(!setflag)
		E->node[node] = E->node[node] | (INTX | INTZ | INTY);
	    } /* end node loop */
	setflag=1;
	/* free the structure, not needed anymore */
	ggrd_grdtrack_free_gstruc(E->control.ggrd.temp.d);
	/* 
	   end temperature from GMT grd init
	*/
	/* end grid init branch */
      }	else {



	/* 
	   
	add a linear profile and potentially perturbations to the temperature
	
	*/
	if(E->control.CART3D)
	  {
	    
	  for(i = 1; i <= noy; i++)
	    for(j = 1; j <= nox; j++)
	      for(k = 1; k <= noz; k++)
		{
		  node = k + (j - 1) * noz + (i - 1) * noz * nox;
		  x1 = E->X[1][node];
		  z1 = E->X[3][node];
		  y1 = E->X[2][node];
		  /* linear gradient */
		  if(E->mesh.bottbc){
		    tz = z1* (E->control.TBCtopval - E->control.TBCbotval) + 
		      E->control.TBCbotval;
		  }else{
		    tz = z1* (E->control.TBCtopval - 1) + 1;
		  }
		  E->T[node] += tz;

		  if(E->convection.random_t_init){
		    /* random init */
		    E->T[node] += (-0.5 + drand48()) * E->convection.perturb_mag[0];
		  }else{
		    
		    for(p=0;p < E->convection.number_of_perturbations;p++){ /* perturbations */
		      
		      E->T[node] += E->convection.perturb_mag[p] * 
			sin(M_PI * (1.0 - z1)) * cos(E->convection.perturb_k[p] * M_PI * x1) * 
			((E->mesh.nsd != 3) ? 1.0 : cos(E->convection.perturb_k[p] * M_PI * y1));
		    }
		  }
		  if(!setflag)
		    E->node[node] = E->node[node] | (INTX | INTZ | INTY);
		}
	  setflag=1;
	}
      else if(E->control.Rsphere)
	{

	  con = (E->mesh.noz - 1) / (E->sphere.ro - E->sphere.ri);
	  noz2 = (E->mesh.noz - 1) / 2 + 1;

	  for(i = 1; i <= noy; i++)
	    for(j = 1; j <= nox; j++)
	      for(k = 1; k <= noz; k++)
		{
		  ii = k + E->lmesh.nzs - 1;
		  node = k + (j - 1) * noz + (i - 1) * noz * nox;
		  x1 = E->SX[1][node];
		  z1 = E->SX[3][node];
		  y1 = E->SX[2][node];
		  /* linear */
		  tz = (z1 -  E->sphere.ri)/(E->sphere.ro - E->sphere.ri) * 
		    (E->control.TBCtopval-E->control.TBCbotval) + E->control.TBCbotval;
		  E->T[node] += tz;
		  if(E->convection.random_t_init){
		    /* random init */
		    E->T[node] += (-0.5 + drand48())*E->convection.perturb_mag[0];
		  }else{
		    for(p=0;p < E->convection.number_of_perturbations;p++){ /* perturbations */
		      mm = E->convection.perturb_mm[p];
		      ll = E->convection.perturb_ll[p];
		      con = E->convection.perturb_mag[p];
		      
		      E->T[node] += E->convection.perturb_mag[p] * modified_plgndr_a(ll, mm, x1) * cos(mm * y1) * 
			sin(M_PI * (z1 - E->sphere.ri) / (E->sphere.ro - E->sphere.ri));
		    }
		  }
		  if(!setflag)
		    E->node[node] = E->node[node] | (INTX | INTZ | INTY);
		}
	  setflag =1;
	}						// end for if else if of geometry
	
      }	/* end perturnbation branch  */
      /* 

      now deal with composition

      */
      if(E->control.ggrd.use_comp){ /* read composition init from grid */
	if(!E->control.composition)
	  myerror("comp grd init but no composition control set!?",E);
	if(E->parallel.me==0)  
	  fprintf(stderr,"convection_initial_temperature: using GMT grd files for composition\n");
	if(E->parallel.me > 0){
	  mpi_rc = MPI_Recv(&mpi_inmsg, 1, MPI_INT, (E->parallel.me-1), mpi_tag, 
			    MPI_COMM_WORLD, &mpi_stat);
	}
	if(E->control.slab_slice){ 
	  /* x - z slice */
	  if(ggrd_grdtrack_init_general(FALSE,E->control.ggrd.comp.gfile, 
					char_dummy,"",
					E->control.ggrd.comp.d,
					/* (E->parallel.me==0)*/
					FALSE,FALSE))
	    myerror("grdtrack init error",E);
	}else{
	  /* 3-D  */
	  if(ggrd_grdtrack_init_general(TRUE,E->control.ggrd.comp.gfile, 
					E->control.ggrd.comp.dfile,"",
					E->control.ggrd.comp.d,
					/* (E->parallel.me==0)*/
					FALSE,FALSE))
	    myerror("grdtrack init error",E);
	}
	if(E->parallel.me <  E->parallel.nproc-1){
	  /* 
	     tell the next processor to go ahead 
	  */
	  mpi_rc = MPI_Send(&mpi_success_message, 1, MPI_INT, 
			    (E->parallel.me+1), mpi_tag, MPI_COMM_WORLD);
	}else{
	  fprintf(stderr,"convection_initial_temperature: last processor done with grd init\n");
	}
	for(i=1;i<=noy;i++)  
	  for(j=1;j<=nox;j++) 
	    for(k=1;k<=noz;k++)  {
	      node=k+(j-1)*noz+(i-1)*nox*noz; /* offset */
	      if(E->control.slab_slice){
		/* slab */
		if(E->control.Rsphere) {
		  if(E->SX[1][node] <= E->control.slab_theta_bound)
		    /* spherical interpolation */
		    ggrd_grdtrack_interpolate_xy((double)E->SX[2][node] * ONEEIGHTYOVERPI,
						 (double)E->SX[1][node],
						 E->control.ggrd.comp.d,&tadd,
						 FALSE);
		  else
		    tadd = 0.0;
		}else{		/* cartesian interpolation */
		  if(E->X[2][node] <= E->control.slab_theta_bound){
		    ggrd_grdtrack_interpolate_xy((double)E->X[1][node],
						 (double)E->X[3][node],
						 E->control.ggrd.comp.d,&tadd,
						 FALSE);
		  }else{
		    tadd = 0.0;
		  }
		}
	      }else{
		/* 3-D */
		if(E->control.Rsphere) /* spherical interpolation */
		  ggrd_grdtrack_interpolate_rtp((double)E->SX[3][node],
						(double)E->SX[1][node],
						(double)E->SX[2][node],
						E->control.ggrd.comp.d,&tadd,
						FALSE,FALSE);
		else		/* cartesian interpolation */
		  ggrd_grdtrack_interpolate_xyz((double)E->X[1][node],
						(double)E->X[2][node],
						(double)E->X[3][node],
						E->control.ggrd.comp.d,&tadd,
						FALSE);
	      }

	      E->C[node] =  E->control.ggrd.comp.offset + tadd *  
		E->control.ggrd.comp.scale;
	    }
	/* free the structure, not needed anymore */
	ggrd_grdtrack_free_gstruc(E->control.ggrd.comp.d);
      }	/* end grid cinit */
      /* 
	 
      check T and C range 
      
      */
      if(E->control.check_t_irange)
	for(i=1;i<=E->lmesh.nno;i++){
	  if(E->T[i]>1)E->T[i]=1;
	  if(E->T[i]<0)E->T[i]=0;
	}
      if(E->control.check_c_irange)
	for(i=1;i<=E->lmesh.nno;i++){
	  if(E->C[i]>1)E->C[i]=1;
	  if(E->C[i]<0)E->C[i]=0;
	}
      
      if(E->control.composition)
	convection_initial_markers(E,1);
    }							// end for restart==0

  else if(E->control.restart)
    {
#ifdef USE_GZDIR
      /* restart part */
      if(E->control.gzdir)
	process_restart_tc_gzdir(E);
      else
#endif
	process_restart_tc(E);

    }

  /* make sure temps are assigned as TBC */
  //assign_T_as_TBC(E,1e-4);


  temperatures_conform_bcs(E);


  thermal_buoyancy(E);

  return;
} 

