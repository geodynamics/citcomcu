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
#include "prototypes.h"

/* ===============================
   Initialization of fields .....
   =============================== */

void convection_initial_temperature_ggrd(struct All_variables *E)
{
  int ll, mm, i, j, k, p, node, ii;
  double temp, temp1, temp2, temp3, base, radius, radius2;
  unsigned int type;

  FILE *fp;
  double  x1, y1, z1, con, top_t, bot_t;

  int noz2, nfz, in1, in2, in3, instance, nox, noy, noz,nxs,setflag;
  char input_s[200], output_file[255];
  float weight, para1, plate_velocity, delta_temp, age;

  char *char_dummy="";
  
  /* twb additions */
  double rho_prem;
  char pfile[1000];
  double t1,f1,r1,tgrad,tadd,tz,tmean;

  /* for dealing with several processors */
  MPI_Status mpi_stat;
  int mpi_rc, mpi_tag=1;  
  int mpi_inmsg, mpi_success_message = 1;
	

  const int dims = E->mesh.nsd;

  noy = E->lmesh.noy;
  noz = E->lmesh.noz;
  nox = E->lmesh.nox;


  setflag = 0;

  /* top and bottom temperatures for initial assign (only use if
     temperatures are set, else defaults */
  bot_t = (E->mesh.bottbc) ? E->control.TBCbotval : 1.0;
  top_t = (E->mesh.toptbc) ? E->control.TBCtopval : 1.0;

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

      	tmean = (top_t + bot_t)/2.0 +  E->control.ggrd.temp.offset;
	if(E->parallel.me == 0)
	  fprintf(stderr,"tinit: top: %g bot: %g mean: %g scale: %g PREM %i\n",
		  top_t,bot_t,tmean,E->control.ggrd.temp.scale,
		  E->control.ggrd.temp.scale_with_prem);
	for(i=1;i<=noy;i++)  
	  for(j=1;j<=nox;j++) 
	    for(k=1;k<=noz;k++)  {
	      node=k+(j-1)*noz+(i-1)*nox*noz; /* offset */
	      if(E->control.slab_slice){
		/* 

		slab slice input 

		*/
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
		   
		3-D grid model input

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
		  /* linear gradient assuming unity box height */
		  tz = z1 * (top_t - bot_t) + bot_t;
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
		  tz = (z1 -  E->sphere.ri)/(E->sphere.ro - E->sphere.ri) * (top_t-bot_t) + bot_t;
		    
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

#ifdef CITCOM_ALLOW_ANISOTROPIC_VISC
/*


read in anisotropic viscosity from a directory which holds


vis2.grd for the viscosity factors, read in log10(eta_S/eta)

nr.grd, nt.grd, np.grd for the directors

*/
void ggrd_read_anivisc_from_file(struct All_variables *E)
{
  MPI_Status mpi_stat;
  int mpi_rc;
  int mpi_inmsg, mpi_success_message = 1;
  int el,i,j,k,l,inode,i1,i2,elxlz,elxlylz,ind,nel;
  int llayer,nox,noy,noz,level,lselect,idim,elx,ely,elz;
  char gmt_string[10],char_dummy;
  double vis2,log_vis,ntheta,nphi,nr;
  float rout[3],xloc[3];
  double cvec[3];
  float base[9];
  char tfilename[1000];
  static ggrd_boolean shift_to_pos_lon = FALSE;
  const int dims=E->mesh.nsd;
  const int ends = enodes[dims];
  FILE *in;
  struct ggrd_gt *vis2_grd,*ntheta_grd,*nphi_grd,*nr_grd;
  const int vpts = vpoints[E->mesh.nsd];

  
  nox=E->mesh.nox;noy=E->mesh.noy;noz=E->mesh.noz;
  elx=E->lmesh.elx;elz=E->lmesh.elz;ely=E->lmesh.ely;
  elxlz = elx * elz;
  elxlylz = elxlz * ely;


  if(E->viscosity.allow_anisotropic_viscosity == 0)
    myerror("ggrd_read_anivisc_from_file: called, but allow_anisotropic_viscosity is FALSE?!",E);
  
  /* isotropic default */
  for(i=E->mesh.levmin;i <= E->mesh.levmax;i++){
    nel  = E->lmesh.NEL[i];
    for(k=1;k <= nel;k++){
      for(l=1;l <= vpts;l++){ /* assign to all integration points */
	ind = (k-1)*vpts + l;
	E->EVI2[i][ind] = 0.0;
	E->EVIn1[i][ind] = 1.0; E->EVIn2[i][ind] = E->EVIn3[i][ind] = 0.0;
      }
    }
  }
  sprintf(gmt_string,"");	/* regional */

  /*
    
  initialization steps
  
  */
  if(E->parallel.me > 0)	/* wait for previous processor */
    mpi_rc = MPI_Recv(&mpi_inmsg, 1, MPI_INT, (E->parallel.me-1),
		      0, MPI_COMM_WORLD, &mpi_stat);
  /*
    read in the material file(s)
  */
  vis2_grd =   (struct  ggrd_gt *)calloc(1,sizeof(struct ggrd_gt));
  nphi_grd =   (struct  ggrd_gt *)calloc(1,sizeof(struct ggrd_gt));
  nr_grd =     (struct  ggrd_gt *)calloc(1,sizeof(struct ggrd_gt));
  ntheta_grd = (struct  ggrd_gt *)calloc(1,sizeof(struct ggrd_gt));


  if(E->parallel.me==0)
    if(E->viscosity.anivisc_layer > 0)
      fprintf(stderr,"ggrd_read_anivisc_from_file: initializing, assigning to all elements above %g km, input is %s\n",
	      E->monitor.length_scale*E->viscosity.zbase_layer[E->viscosity.anivisc_layer - 1],
	      ("regional"));
    else
      fprintf(stderr,"ggrd_read_anivisc_from_file: initializing, assigning to all elements between  %g and %g km, input is %s\n",
	      E->monitor.length_scale*((E->viscosity.anivisc_layer<-1)?(E->viscosity.zbase_layer[-E->viscosity.anivisc_layer - 2]):(0)),
	      E->monitor.length_scale*E->viscosity.zbase_layer[-E->viscosity.anivisc_layer - 1],
	      ("regional"));

  /* 
     read viscosity ratio, and east/north direction of normal azimuth 
  */
  /* viscosity factor */
  sprintf(tfilename,"%s/vis2.grd",E->viscosity.anisotropic_init_dir);
  if(ggrd_grdtrack_init_general(FALSE,tfilename,"",gmt_string,
				vis2_grd,(E->parallel.me == 0),FALSE))
    myerror("ggrd init error",E);
  /* n_r */
  if(E->control.CART3D)
    sprintf(tfilename,"%s/nz.grd",E->viscosity.anisotropic_init_dir);
  else
    sprintf(tfilename,"%s/nr.grd",E->viscosity.anisotropic_init_dir);

  if(ggrd_grdtrack_init_general(FALSE,tfilename,"",gmt_string,
				nr_grd,(E->parallel.me == 0),FALSE))
    myerror("ggrd init error",E);
  if(E->control.CART3D)
    /* n_theta */
    sprintf(tfilename,"%s/ny.grd",E->viscosity.anisotropic_init_dir);
  else
    sprintf(tfilename,"%s/nt.grd",E->viscosity.anisotropic_init_dir);

  if(ggrd_grdtrack_init_general(FALSE,tfilename,"",gmt_string,
				ntheta_grd,(E->parallel.me == 0),FALSE))
    myerror("ggrd init error",E);
  if(E->control.CART3D)
  /* n_phi */
    sprintf(tfilename,"%s/nx.grd",E->viscosity.anisotropic_init_dir);
  else
    sprintf(tfilename,"%s/np.grd",E->viscosity.anisotropic_init_dir);

  if(ggrd_grdtrack_init_general(FALSE,tfilename,"",gmt_string,
				nphi_grd,(E->parallel.me == 0),FALSE))
    myerror("ggrd init error",E);

  /* done */
  if(E->parallel.me <  E->parallel.nproc-1){ /* tell the next proc to go ahead */
    mpi_rc = MPI_Send(&mpi_success_message, 1,
		      MPI_INT, (E->parallel.me+1), 0, MPI_COMM_WORLD);
  }else{
    fprintf(stderr,"ggrd_read_anivisc_from_file: last processor done with ggrd init\n");
    fprintf(stderr,"ggrd_read_anivisc_from_file: WARNING: assuming a regular grid geometry\n");
  }
  /*

  loop through all elements and assign

  */
  for (j=1;j <= elz;j++)  {	/* this assumes a regular grid sorted as in (1)!!! */
    if(((E->viscosity.anivisc_layer > 0)&&(E->mat[j] <=   E->viscosity.anivisc_layer))||
       ((E->viscosity.anivisc_layer < 0)&&(E->mat[j] ==  -E->viscosity.anivisc_layer))){
      /* within top layers */
      for (k=1;k <= ely;k++){
	for (i=1;i <= elx;i++)   {
	  /* eq.(1) */
	  el = j + (i-1) * elz + (k-1)*elxlz;
	  /*
	    find average coordinates
	  */
	  if(E->control.CART3D){
	    rout[0]=rout[1]=rout[2]=0.0;
	    for(inode=1;inode <= 4;inode++){
	      ind = E->ien[el].node[inode];
	      rout[0] += E->X[1][ind];
	      rout[1] += E->X[2][ind];
	      rout[2] += E->X[3][ind];
	    }
	    rout[0]/=4.;rout[1]/=4.;rout[2]/=4.;
	    if(!ggrd_grdtrack_interpolate_xy((double)rout[0],(double)rout[1],vis2_grd,&log_vis,FALSE)){
	      fprintf(stderr,"ggrd_read_anivisc_from_file: interpolation error at x: %g y: %g\n",
		      rout[0],rout[1]);
	      parallel_process_termination();
	    }
	    /* read in cartesian vector */
	    if(!ggrd_grdtrack_interpolate_xy((double)rout[0],(double)rout[1],nphi_grd,cvec,FALSE)){
	      fprintf(stderr,"ggrd_read_anivisc_from_file: interpolation error at x: %g y: %g\n",
		      rout[0],rout[1]);
	      parallel_process_termination();
	    }
	    if(!ggrd_grdtrack_interpolate_xy((double)rout[0],(double)rout[1],ntheta_grd,(cvec+1),FALSE)){
	      fprintf(stderr,"ggrd_read_anivisc_from_file: interpolation error at x: %g y: %g\n",
		      rout[0],rout[1]);
	      parallel_process_termination();
	    }
	    if(!ggrd_grdtrack_interpolate_xy((double)rout[0],(double)rout[1],nr_grd,(cvec+2),FALSE)){
	      fprintf(stderr,"ggrd_read_anivisc_from_file: interpolation error at x: %g y: %g\n",
		      rout[0],rout[1]);
	      parallel_process_termination();
	    }
	    normalize_vec3d(cvec,(cvec+1),(cvec+2));
	  }else{
	    /* spherical */
	    xloc[0] = xloc[1] = xloc[2] = 0.0;
	    for(inode=1;inode <= 4;inode++){
	      ind = E->ien[el].node[inode];
	      rtp2xyz((float)E->SX[3][ind],(float)E->SX[1][ind],
		      (float)E->SX[2][ind],rout);
	      xloc[0] += rout[0];
	      xloc[1] += rout[1];
	      xloc[2] += rout[2];
	    }
	    xloc[0]/=4.;xloc[1]/=4.;xloc[2]/=4.;
	    xyz2rtp(xloc[0],xloc[1],xloc[2],rout); /* convert to spherical */
	    
	    /* vis2 */
	    if(!ggrd_grdtrack_interpolate_tp(rout[1],rout[2],
					     vis2_grd,&log_vis,FALSE,shift_to_pos_lon)){
	      fprintf(stderr,"ggrd_read_anivisc_from_file: interpolation error at lon: %g lat: %g\n",
		      rout[2]*180/M_PI,90-rout[1]*180/M_PI);
	      parallel_process_termination();
	    }
	    /* nr */
	    if(!ggrd_grdtrack_interpolate_tp(rout[1],rout[2],
					     nr_grd,&nr,FALSE,shift_to_pos_lon)){
	      fprintf(stderr,"ggrd_read_anivisc_from_file: interpolation error at lon: %g lat: %g\n",
		      rout[2]*180/M_PI,90-rout[1]*180/M_PI);
	      parallel_process_termination();
	    }
	    /* ntheta */
	    if(!ggrd_grdtrack_interpolate_tp(rout[1],rout[2],
					     ntheta_grd,&ntheta,FALSE,shift_to_pos_lon)){
	      fprintf(stderr,"ggrd_read_anivisc_from_file: interpolation error at lon: %g lat: %g\n",
		      rout[2]*180/M_PI,90-rout[1]*180/M_PI);
	      parallel_process_termination();
	    }
	    /* nphi */
	    if(!ggrd_grdtrack_interpolate_tp(rout[1],rout[2],
					     nphi_grd,&nphi,FALSE,shift_to_pos_lon)){
	      fprintf(stderr,"ggrd_read_anivisc_from_file: interpolation error at lon: %g lat: %g\n",
		      rout[2]*180/M_PI,90-rout[1]*180/M_PI);
	      parallel_process_termination();
	    }
	    normalize_vec3d(&nr,&ntheta,&nphi);
	    calc_cbase_at_tp(rout[1],rout[2],base); /* convert from
						       spherical
						       coordinates to
						       Cartesian */
	    convert_pvec_to_cvec(nr,ntheta,nphi,base,xloc);
	    cvec[0]=xloc[0];cvec[1]=xloc[1];cvec[2]=xloc[2];
	  }
	  /* transform */
	  vis2 = 1.0 - pow(10.0,log_vis);
	  for(l=1;l <= vpts;l++){ /* assign to all integration points */
	    ind = (el-1)*vpts + l;
	    E->EVI2[E->mesh.levmax][ind]  =   vis2;
	    E->EVIn1[E->mesh.levmax][ind]  = cvec[0];
	    E->EVIn2[E->mesh.levmax][ind]  = cvec[1];
	    E->EVIn3[E->mesh.levmax][ind]  = cvec[2];
	  }
	}
      }
    }	/* end insize lith */
  }	/* end elz loop */


  ggrd_grdtrack_free_gstruc(vis2_grd);
  ggrd_grdtrack_free_gstruc(nr_grd);
  ggrd_grdtrack_free_gstruc(ntheta_grd);
  ggrd_grdtrack_free_gstruc(nphi_grd);
  
  
}


#endif	/* for ANISOTROPIC */

