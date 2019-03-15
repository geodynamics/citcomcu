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




 minor changes by twb to change output format

*/

/* Routine to process the output of the finite element cycles 
   and to turn them into a coherent suite  files  */


#include <fcntl.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>				/* for "system" command */
#include <zlib.h>



#include "element_definitions.h"
#include "global_defs.h"

#ifdef CITCOM_ALLOW_ANISOTROPIC_VISC
#include "anisotropic_viscosity.h"
#endif

gzFile safe_gzopen(char *,char *);
FILE *safe_fopen(char *,char *);
void calc_cbase_at_tp(float , float , float *);
void convert_pvec_to_cvec(float ,float , float , float *,float *);
void convert_cvec_to_pvec(float ,float , float , float *,float *);
void rtp2xyz(float r, float , float , float *);


void output_velo_related_gzdir(E, file_number)
     struct All_variables *E;
     int file_number;
{
  int el, els, i, j, k, ii, m, node, fd,snode;
  int nox, noz, noy, nfx, nfz, nfy1, nfy2, size1, size2,offset;
  static int vtkout = 1;
  char output_file[1000];
  float cvec[3],locx[3],*fdummy,x[4];
  double rtf[4][9];
  double VV[4][9],lgrad[3][3],isa[3],evel[3];
  static struct CC Cc;
  static struct CCX Ccx;

  static float *SV, *EV;
  static int been_here = 0;
  float vs;
  int vtk_comp_out;

  float *SZZ, *SXX, *SYY, *SXY, *SXZ, *SZY;
  
  const int dims = E->mesh.nsd;
  const int ends = enodes[dims];
  const int vpts = vpoints[dims];
  const int nno = E->mesh.nno;
  const int nel = E->lmesh.nel;
  const int lev = E->mesh.levmax;
  const int ppts = ppoints[dims];
  float *e2,*e2n, Pressure;
  gzFile gzout;

  if(E->control.composition)
    vtk_comp_out = 1;
  else
    vtk_comp_out = 0;
  
  /* make a directory */
  if(E->parallel.me == 0){
    fprintf(stderr,"Output_gzdir: output directory");
    //sprintf(output_file,"if [ ! -s %s/%d ];then mkdir -p %s/%d;fi 2> /dev/null",
    sprintf(output_file,"if [ ! -s %s/%d/ ];then mkdir -p %s/%d/;fi ",
	    E->control.data_file2,file_number,E->control.data_file2,file_number);
    fprintf(stderr," %s...",output_file);
    system(output_file);
    fprintf(stderr,"mkdir done\n");
  }
  /* and wait for the other jobs */
  parallel_process_sync();
  if(E->debug && (E->parallel.me==0))
    fprintf(stderr,"output_velo_related_gzdir: directory done\n");

  if(been_here == 0){
    /* 
       
    initial nodal coordinate output 
    
    */
    if(E->control.COMPRESS){
      /* compressed I/O */
      /* 
	 all nodes 
      */
      sprintf(output_file, "%s/coord.%d.gz", E->control.data_file2, E->parallel.me);
      gzout = safe_gzopen(output_file, "w");
      gzprintf(gzout, "%6d\n", E->lmesh.nno);
      if(!E->control.Rsphere)
	for(i = 1; i <= E->lmesh.nno; i++){
	  gzprintf(gzout, "%13.6e %13.6e %13.6e\n", E->X[1][i], E->X[2][i], E->X[3][i]);
	}
      else
	for(i = 1; i <= E->lmesh.nno; i++)
	  gzprintf(gzout, "%13.6e %13.6e %13.6e\n", E->SX[1][i], E->SX[2][i], E->SX[3][i]);
      gzclose(gzout);
      /* 
	 surface nodes 
      */
      sprintf(output_file, "%s/scoord.%d.gz", E->control.data_file2, E->parallel.me);
      gzout = safe_gzopen(output_file, "w");
      gzprintf(gzout, "%6d\n", E->lmesh.nsf);
      for(i = 1; i <= E->lmesh.nsf; i++){
	node = E->surf_node[i];
	if(!E->control.Rsphere)	/* leave in redundancy for now */
	  gzprintf(gzout, "%13.6e %13.6e %13.6e\n", E->X[1][node], E->X[2][node], E->X[3][node]);
	else
	  gzprintf(gzout, "%13.6e %13.6e %13.6e\n", E->SX[1][node], E->SX[2][node], E->SX[3][node]);
      }
      gzclose(gzout);
    }else{			/* regular I/O */
      /* 
	 
      all nodes

      */
      sprintf(output_file, "%s/coord.%d", E->control.data_file2, E->parallel.me);
      E->filed[13] = safe_fopen(output_file, "w");
      fprintf(E->filed[13], "%6d\n", E->lmesh.nno);
      if(!E->control.Rsphere)
	for(i = 1; i <= E->lmesh.nno; i++)
	  fprintf(E->filed[13], "%13.6e %13.6e %13.6e\n", E->X[1][i], E->X[2][i], E->X[3][i]);
      else
	for(i = 1; i <= E->lmesh.nno; i++)
	  fprintf(E->filed[13], "%13.6e %13.6e %13.6e\n", E->SX[1][i], E->SX[2][i], E->SX[3][i]);
      fclose(E->filed[13]);
      if((E->parallel.me_loc[3] == E->parallel.nprocz - 1) || /* top */ 
	 (E->parallel.me_loc[3] == 0)){ /* bottom */
	/* 
	   surface nodes 
	*/
	sprintf(output_file, "%s/scoord.%d", E->control.data_file2, E->parallel.me);
	E->filed[13] = safe_fopen(output_file, "w");
	fprintf(E->filed[13], "%6d\n", E->lmesh.nsf);
	for(i = 1; i <= E->lmesh.nsf; i++){
	  node = E->surf_node[snode];
	  if(!E->control.Rsphere)	/* leave in redundancy for now */
	    fprintf(E->filed[13], "%13.6e %13.6e %13.6e\n", E->X[1][node], E->X[2][node], E->X[3][node]);
	  else
	    fprintf(E->filed[13], "%13.6e %13.6e %13.6e\n", E->SX[1][node], E->SX[2][node], E->SX[3][node]);
	}
	fclose(E->filed[13]);
      }
    } 
    if(vtkout){
      /* 
	 vtk output, always compressed
      */
      /* 
	   node coordinates 
      */
      
      
      sprintf(output_file,"%s/vtk_ecor.%d.gz",E->control.data_file2,E->parallel.me);
      gzout = safe_gzopen(output_file,"w");
      if(E->control.Rsphere){	/* spherical */
	//if(E->parallel.me == 0)fprintf(stderr, "converting spherical to cartesian  vtk coords to %s \n",output_file);
	for(i=1;i <= E->lmesh.nno;i++) {
	  rtp2xyz((float)E->SX[3][i],(float)E->SX[1][i],(float)E->SX[2][i],locx);
	  gzprintf(gzout,"%9.6f %9.6f %9.6f\n",locx[0],locx[1],locx[2]);
	}
      }else{			/* cartesian */
	if(E->parallel.me == 0)
	  fprintf(stderr, "writing cartesian vtk coords to %s \n",output_file);
	for(i=1;i <= E->lmesh.nno;i++){
	  gzprintf(gzout,"%9.6f %9.6f %9.6f\n",E->X[1][i],E->X[2][i],E->X[3][i]);
	}
      }
      gzclose(gzout);
      
      /* 
	 connectivity 
      */
      offset = E->lmesh.nno * E->parallel.me - 1;
      sprintf(output_file,"%s/vtk_econ.%d.gz",E->control.data_file2,E->parallel.me);
      gzout = safe_gzopen(output_file,"w");
      for(i=1;i<= E->lmesh.nel;i++) {
	gzprintf(gzout,"%2i\t",ends);
	/* 
	   need to add offset according to the processor for global
	   node numbers
	*/
	gzprintf(gzout,"%6i %6i %6i %6i %6i %6i %6i %6i\n",
		 E->ien[i].node[1]+offset,E->ien[i].node[2]+offset,
		 E->ien[i].node[3]+offset,E->ien[i].node[4]+offset,
		 E->ien[i].node[5]+offset,E->ien[i].node[6]+offset,
		 E->ien[i].node[7]+offset,E->ien[i].node[8]+offset);
      }
      gzclose(gzout);
      //if(E->parallel.me == 0)fprintf(stderr, "writing element connectivity %s \n",output_file);
      /* end vtk branch */
    }
    /* end init loop */
    been_here++;
    if(E->debug && E->parallel.me==0)
      fprintf(stderr,"initial coord out done\n");
  
  }
  
  if(E->parallel.me < E->parallel.nprocz)
    {
      sprintf(output_file, "%s/%d/ave_r.%d.%d", E->control.data_file2, file_number, 
	      E->parallel.me, file_number);
      E->filed[10] = safe_fopen(output_file, "w");
      /* header line */
      fprintf(E->filed[10], "# %6d %6d %13.6e %13.6e %13.6e %12.4e %12.4e %13.6e %13.6e\n", 
	      E->lmesh.noz, E->advection.timesteps, E->monitor.elapsed_time, 
	      E->slice.Nut, E->slice.Nub, E->data.T_adi0, E->data.T_adi1, 
	      E->monitor.Sigma_interior, E->monitor.Sigma_max);
      /* 
	 added first column of radial coordinate TWB 
      */
      if(!E->control.Rsphere){	/* cartesian */
	for(i = 1; i <= E->lmesh.noz; i++){
	  fprintf(E->filed[10], "%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n", 
		  E->X[3][i],E->Have.T[i], E->Have.vrms[i], E->Have.Vi[i], E->Have.Rho[i], 
		  E->Have.F[i], E->Have.f[i], E->Have.C[i], E->Have.Tadi[i]);
	}
      }else{			/* spherical */
	for(i = 1; i <= E->lmesh.noz; i++){
	  fprintf(E->filed[10], "%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n", 
		  E->SX[3][i],E->Have.T[i], E->Have.vrms[i], E->Have.Vi[i], E->Have.Rho[i], 
		  E->Have.F[i], E->Have.f[i], E->Have.C[i], E->Have.Tadi[i]);
	}
      }
      fclose(E->filed[10]);
      if(E->debug && E->parallel.me==0)
	fprintf(stderr,"averages  done\n");

    }
  


  if((file_number % E->control.record_every) == 0)
    {
      if(vtkout){
	/* 
	   new vtk output
	
	*/
	/* temperature at nodes */
	sprintf(output_file,"%s/%d/t.%d.%d.gz",
		E->control.data_file2,file_number,E->parallel.me,file_number);
	gzout=safe_gzopen(output_file,"w");
	gzprintf(gzout,"%d %d %13.6e\n",
		 file_number,E->lmesh.nno,
		 E->monitor.elapsed_time);
	gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno); /* for backward compatibility */
	for(i=1;i<=E->lmesh.nno;i++)           
	  gzprintf(gzout,"%.6e\n",E->T[i]);
	gzclose(gzout);
	if(E->debug && E->parallel.me==0)
	  fprintf(stderr,"temp out done\n");

	/* velocities */
	if(E->control.Rsphere && (!E->sphere.vtk_base_init)){ 
	  force_report(E,"constructing basis function");
	  /* need to init spherical/cartesian conversion vectors */
	  E->sphere.vtk_base = (float *)safe_malloc(sizeof(float)*E->lmesh.nno*9);
	  for(k=0,i=1;i<=E->lmesh.nno;i++,k+=9)
	    calc_cbase_at_tp(E->SX[1][i],E->SX[2][i],(E->sphere.vtk_base+k));
	  E->sphere.vtk_base_init = 1;
	}
	/* write velocities to file */
	sprintf(output_file,"%s/%d/vtk_v.%d.%d.gz",
		E->control.data_file2,file_number, E->parallel.me,file_number);

	/* 
	 */
	gzout=safe_gzopen(output_file,"w");
	if(E->control.Rsphere){
	  for(k=0,i=1;i<=E->lmesh.nno;i++,k+=9) {
	    /* 
	       convert r,theta,phi vector to x,y,z at base location, output is in 
	       cm/yr
	    */
	    convert_pvec_to_cvec(E->V[3][i],E->V[1][i],E->V[2][i],(E->sphere.vtk_base+k),cvec);
	    /* output of cartesian vector */
	    gzprintf(gzout,"%.6e %.6e %.6e\n",cvec[0],cvec[1],cvec[2]);
	  }
	}else{			/* cartesian output in cm/yr */
	  for(i=1;i<=E->lmesh.nno;i++) {
	    gzprintf(gzout,"%.6e %.6e %.6e\n",E->V[1][i],E->V[2][i],E->V[3][i]);
	  }
	}
	gzclose(gzout);
	if(E->debug && E->parallel.me==0)
	  fprintf(stderr,"vel out done\n");

	if(E->control.ele_pressure_out){
	  if(E->control.Rsphere)
	    myerror("ele_pressure_out not set up for spherical geometry",E);
	  /*  */
	  sprintf(output_file,"%s/%d/pe.%d.%d.gz",E->control.data_file2,
		  file_number, E->parallel.me,file_number);
	  gzout = safe_gzopen(output_file,"w");
	  gzprintf(gzout,"%d %d %13.6e\n",file_number,E->lmesh.nno,E->monitor.elapsed_time);
	  for(i=1;i <= nel;i++){
	    /* get element centroid location */
	    x[1]=x[2]=x[3]=0;
	    for(j=1;j <= ends;j++){
	      node = E->ien[i].node[j]; /* global node number */
	      for(k=1;k <= 3;k++)
		x[k] += E->X[k][node];
	    }
	    for(k=1;k <= 3;k++)
	      x[k] /= (float)ends;
	    gzprintf(gzout,"%12.8f %12.8f %12.8f %.6e\n",x[1],x[2],x[3],E->P[i]);
	  } /* end element loop */
	  gzclose(gzout);
	}


	if(E->control.vtk_pressure_out){
	  /* 
	     pressure at nodes 
	  */
	  sprintf(output_file,"%s/%d/p.%d.%d.gz",E->control.data_file2,
		  file_number, E->parallel.me,file_number);
	  gzout = safe_gzopen(output_file,"w");
	  gzprintf(gzout,"%d %d %13.6e\n",file_number,E->lmesh.nno,E->monitor.elapsed_time);
	  gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	  for(i=1;i<=E->lmesh.nno;i++)          
	    gzprintf(gzout,"%.6e\n",E->NP[i]);
	  gzclose(gzout);
	}
	
        if(E->control.vtk_stress2_out){    
	  force_report(E,"stress output, large mem demand...");
	  /* begin Adam style stress output */
	  SXX = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	  SYY = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	  SXY = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	  SXZ = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	  SZY = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	  SZZ = (float *)malloc((E->lmesh.nno + 1) * sizeof(float));
	  /* compute the stresses */
	  get_stress(SXX,SXY,SXZ,SYY,SZY,SZZ,E);
	  /* don't use isotropic component */
	  remove_isotropic_component(SXX,SYY,SZZ,E->lmesh.nno);
	  /* output */
	  /* s_xx */
	  sprintf(output_file,"%s/%d/sXX.%d.%d.gz",E->control.data_file2,
		  file_number, E->parallel.me,file_number);
	  gzout = safe_gzopen(output_file,"w");
	  gzprintf(gzout,"%d %d %13.6e\n",file_number,E->lmesh.nno,
		   E->monitor.elapsed_time);
	  gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	  for(i=1;i<=E->lmesh.nno;i++)
	    gzprintf(gzout,"%.6e\n",SXX[i]);
	  gzclose(gzout);
	  /* s_zz */
	  sprintf(output_file,"%s/%d/sZZ.%d.%d.gz",E->control.data_file2,
		  file_number, E->parallel.me,file_number);
	  gzout = safe_gzopen(output_file,"w");
	  gzprintf(gzout,"%d %d %13.6e\n",file_number,E->lmesh.nno,
		   E->monitor.elapsed_time);
	  gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	  for(i=1;i<=E->lmesh.nno;i++)
	    gzprintf(gzout,"%.6e\n",SZZ[i]);
	  gzclose(gzout);
	  /* s_xz */
	  sprintf(output_file,"%s/%d/sXZ.%d.%d.gz",E->control.data_file2,
		  file_number, E->parallel.me,file_number);
	  gzout = safe_gzopen(output_file,"w");
	  gzprintf(gzout,"%d %d %13.6e\n",file_number,E->lmesh.nno,
		   E->monitor.elapsed_time);
	  gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	  for(i=1;i<=E->lmesh.nno;i++)
	    gzprintf(gzout,"%.6e\n",SXZ[i]);
	  gzclose(gzout);

	  if(E->control.vtk_stress_3D){
	    /* s_xy */
	    sprintf(output_file,"%s/%d/sXY.%d.%d.gz",
		    E->control.data_file2,file_number, E->parallel.me,file_number);
	    gzout = safe_gzopen(output_file,"w");
	    gzprintf(gzout,"%d %d %13.6e\n",file_number,
		     E->lmesh.nno,E->monitor.elapsed_time);
	    gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	    for(i=1;i<=E->lmesh.nno;i++)
	      gzprintf(gzout,"%.6e\n",SXY[i]);
	    gzclose(gzout);
	    
	    /* s_yy */
	    sprintf(output_file,"%s/%d/sYY.%d.%d.gz",
		    E->control.data_file2,file_number, E->parallel.me,file_number);
	    gzout = safe_gzopen(output_file,"w");
	    gzprintf(gzout,"%d %d %13.6e\n",file_number,
		     E->lmesh.nno,E->monitor.elapsed_time);
	    gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	    for(i=1;i<=E->lmesh.nno;i++)
	      gzprintf(gzout,"%.6e\n",SYY[i]);
	    gzclose(gzout);
	    /* s_yz */
	    sprintf(output_file,"%s/%d/sYZ.%d.%d.gz",
		    E->control.data_file2,file_number, E->parallel.me,file_number);
	    gzout = safe_gzopen(output_file,"w");
	    gzprintf(gzout,"%d %d %13.6e\n",file_number,
		     E->lmesh.nno,E->monitor.elapsed_time);
	    gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	    for(i=1;i<=E->lmesh.nno;i++)
	      gzprintf(gzout,"%.6e\n",SZY[i]);
	    gzclose(gzout);


	    

	  }


	  
	  free((void *)SXX);
	  free((void *)SYY);
	  free((void *)SXY);
	  free((void *)SXZ);
	  free((void *)SZY);
	  free((void *)SZZ);

	  /* end Adam stress addition */
        }

	if(E->control.vtk_e2_out){
	  /* 
	     second invariant 
	  */
	  force_report(E,"second invariant output...");
	  e2n = (float *)malloc((E->lmesh.nno + 10) * sizeof(float));
	  e2 = (float *)malloc((E->lmesh.nel + 1) * sizeof(float));
	  strain_rate_2_inv(E, e2, 1,0,fdummy); /* compute strain rate for every element */
	  e2_to_nodes(E, e2,e2n, E->mesh.levmax); /* project to nodes */
	  sprintf(output_file,"%s/%d/e2inv.%d.%d.gz",E->control.data_file2,
		  file_number, E->parallel.me,file_number);
	  gzout = safe_gzopen(output_file,"w");
	  gzprintf(gzout,"%d %d %13.6e\n",file_number,E->lmesh.nno,E->monitor.elapsed_time);
	  gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	  for(i=1;i<=E->lmesh.nno;i++)          
	    gzprintf(gzout,"%.6e\n",e2n[i]);
	  gzclose(gzout);
	  free(e2n);free(e2);
	}
	if(E->control.vtk_vgm_out){
	  /* 
	     
	     print velocity gradient matrix  and ISA

	  */
	  sprintf(output_file,"%s/%d/edot.%d.%d.gz",
		  E->control.data_file2,
		  file_number, E->parallel.me,file_number);
	  gzout = safe_gzopen(output_file,"w");
	  for(i=1;i <= E->lmesh.nel;i++){ /* element loop */
	    if(E->control.Rsphere){	/* need rtf for spherical */
	      get_rtf(E, i, 1, rtf, lev); /* pressure points */
	      if((i-1)%E->lmesh.elz==0)
		construct_c3x3matrix_el(E,i,&Cc,&Ccx,lev,1);
	    }
	    for(j = 1; j <= ends; j++){	/* velocity at element nodes */
	      node = E->ien[i].node[j];
	      VV[1][j] = E->V[1][node];
	      VV[2][j] = E->V[2][node];
	      VV[3][j] = E->V[3][node];
	    }
	    /* find mean location */
	    locx[0] = locx[1] = locx[2] = 0.0;
	    for(j = 1; j <= ends; j++){
	      node = E->ien[i].node[j];
	      if(E->control.Rsphere){
		rtp2xyz((float)E->SX[3][node],
			(float)E->SX[1][node],
			(float)E->SX[2][node],cvec);
		locx[0] += cvec[0];
		locx[1] += cvec[1];
		locx[2] += cvec[2];
	      }else{
		locx[0] += E->X[1][node];
		locx[1] += E->X[2][node];
		locx[2] += E->X[3][node];
	      }
	    }
	    locx[0] /= ends;locx[1] /= ends;locx[2] /= ends;
	    /* calculate velocity gradient matrix */
	    get_vgm_p(VV,&(E->N),&(E->GNX[lev][i]),&Cc, &Ccx,rtf,
		      E->mesh.nsd,ppts,ends,(E->control.Rsphere),lgrad,
		      evel);
#ifdef  CITCOM_ALLOW_ANISOTROPIC_VISC
	    /* calculate the ISA axis or whichever was used to align */
	    switch(E->viscosity.anisotropic_init){
	    case 3:
	      calc_isa_from_vgm(lgrad,evel,i,isa,E,CITCOM_ANIVISC_ALIGN_WITH_VEL);
	      break;
	    case 4:
	      calc_isa_from_vgm(lgrad,evel,i,isa,E,CITCOM_ANIVISC_ALIGN_WITH_ISA);
	      break;
	    case 5:
	      calc_isa_from_vgm(lgrad,evel,i,isa,E,CITCOM_ANIVISC_MIXED_ALIGN);
	      break;
	    default:
	      calc_isa_from_vgm(lgrad,evel,i,isa,E,CITCOM_ANIVISC_ALIGN_WITH_ISA);
	      break;
	    }
	    normalize_vec3d(evel,(evel+1),(evel+2));
#else
	    isa[0]=isa[1]=isa[2]=evel[0]=evel[1]=evel[2]=0.0;
#endif	    
	    /*  */
	    if(E->control.Rsphere){
	      /* print velocity gradient matrix, ISA, and normalized velocities */
	      xyz2rtp(locx[0],locx[1],locx[2],cvec);
	      /* r theta phi d[upper right half] (CitcomS spherical system) ISA[3] */
	      gzprintf(gzout,"%11g %11g %11g\t%13.6e %13.6e %13.6e %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e\t%13.6e %13.6e %13.6e\t%13.6e %13.6e %13.6e\n",
		       cvec[0],cvec[1],cvec[2],
		       lgrad[0][0],lgrad[0][1],lgrad[0][2],
		       lgrad[1][0],lgrad[1][1],lgrad[1][2],
		       lgrad[2][0],lgrad[2][1],lgrad[2][2],
		       isa[0],isa[1],isa[2],evel[0],evel[1],evel[2]);
	    }else{
	      gzprintf(gzout,"%11g %11g %11g\t%13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\t%13.6e %13.6e %13.6e\t%13.6e %13.6e %13.6e\n",
		       locx[0],locx[1],locx[2],
		       lgrad[0][0],lgrad[0][1],lgrad[0][2],
		       lgrad[1][0],lgrad[1][1],lgrad[1][2],
		       lgrad[2][0],lgrad[2][1],lgrad[2][2],
		       isa[0],isa[1],isa[2],evel[0],evel[1],evel[2]);
	    }
	  } /* end element loop */
	  gzclose(gzout);
	}
	
	if(E->control.vtk_viscosity_out){
	  /* 
	     
	     viscosity output
	  
	  */
	  /* write to file */
	  sprintf(output_file,"%s/%d/visc.%d.%d.gz",
		  E->control.data_file2,file_number, E->parallel.me,file_number);
	  gzout=safe_gzopen(output_file,"w");
	  gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	  for(i=1;i<=E->lmesh.nno;i++){
	    gzprintf(gzout,"%12.4e\n",E->VI[lev][i]);
	  }
	  gzclose(gzout);
#ifdef CITCOM_ALLOW_ANISOTROPIC_VISC
	  if(E->viscosity.allow_anisotropic_viscosity){
	    sprintf(output_file,"%s/%d/avisc.%d.%d.gz",
		    E->control.data_file2,file_number, E->parallel.me,file_number);
	    gzout=safe_gzopen(output_file,"w");
	    gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	    for(i=1;i<=E->lmesh.nno;i++){
	      gzprintf(gzout,"%10.4e %10.4e %10.4e %10.4e\n",
		       E->VI2[lev][i],
		       E->VIn1[lev][i],
		       E->VIn2[lev][i],
		       E->VIn3[lev][i]);
	    }
	    gzclose(gzout);


	  }
#endif
	}
	if(E->control.vtk_ddpart_out){
	  /* 
	     
	     diffusion/dislocation partitioning
	  
	  */
	  /* write to file */
	  sprintf(output_file,"%s/%d/ddpart.%d.%d.gz",
		  E->control.data_file2,file_number, E->parallel.me,file_number);
	  gzout=safe_gzopen(output_file,"w");
	  gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	  if(E->viscosity.sdepv_rheology == 3){
	    for(i=1;i<=E->lmesh.nno;i++)
	      gzprintf(gzout,"%12.4e\n",E->Vddpart[i]);
	  }else{
	    for(i=1;i<=E->lmesh.nno;i++)
	      gzprintf(gzout,"%12.4e\n",-8);
	  }
	  gzclose(gzout);
	}
	if(vtk_comp_out){
	  /* 
	     
	     composition output
	  
	  */

	  /* write to file */
	  sprintf(output_file,"%s/%d/c.%d.%d.gz",
		  E->control.data_file2,file_number, E->parallel.me,file_number);
	  gzout=safe_gzopen(output_file,"w");
	  gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	  for(i=1;i<=E->lmesh.nno;i++)           
	    gzprintf(gzout,"%12.4e\n",E->C[i]);
	  gzclose(gzout);

	  if(E->tracers_track_strain){
	    sprintf(output_file,"%s/%d/strain.%d.%d.gz",
		    E->control.data_file2,file_number, E->parallel.me,file_number);
	    gzout=safe_gzopen(output_file,"w");
	    gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	    for(i=1;i<=E->lmesh.nno;i++)           
	      gzprintf(gzout,"%12.4e\n",E->strain[i]);
	    gzclose(gzout);
	  }

	  if(E->tracers_add_flavors){
	    sprintf(output_file,"%s/%d/cf.%d.%d.gz",
		    E->control.data_file2,file_number, E->parallel.me,file_number);
	    gzout=safe_gzopen(output_file,"w");
	    gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	    for(i=1;i <= E->lmesh.nno;i++) {
	      for(j=0;j < E->tracers_add_flavors;j++)
		gzprintf(gzout,"%3i ",E->CF[j][i]);
	      gzprintf(gzout,"\n");
	    }
	    gzclose(gzout);
	  }

	}
	/* 
	   end vtk out 
	*/
      }else{
	/*

 
	   old output format 

	*/
	if(E->control.COMPRESS){
	  sprintf(output_file, "%s/temp.%d.%d.gz", E->control.data_file2, E->parallel.me, file_number);
	  gzout = safe_gzopen(output_file, "w");
	  gzprintf(gzout, "%6d %6d %13.6e\n", E->lmesh.nno, E->advection.timesteps, 
		  E->monitor.elapsed_time);
	  if((file_number % 2* E->control.record_every) == 0)
	    {
	      if(!E->control.composition)
		for(i = 1; i <= E->lmesh.nno; i++)
		  //        gzprintf(gzout,"%13.6e %12.4e %12.4e %13.6e %13.6e %13.6e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i],E->V[1][i],E->V[2][i],E->V[3][i]);
		  //       gzprintf(gzout,"%13.6e %12.4e %12.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i]);
		  gzprintf(gzout, "%13.6e %12.4e %12.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i]);
	      else
		for(i = 1; i <= E->lmesh.nno; i++)
		  //      gzprintf(gzout,"%13.6e %12.4e %12.4e %12.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i],E->C[i]);
		  gzprintf(gzout, "%13.6e %12.4e %12.4e %12.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i], E->C[i]);
	    }
	  else
	    {
	      if(!E->control.composition)
		for(i = 1; i <= E->lmesh.nno; i++)
		  //      gzprintf(gzout,"%13.6e %12.4e %12.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i]);
		  gzprintf(gzout, "%13.6e %12.4e %12.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i]);
	      else
		for(i = 1; i <= E->lmesh.nno; i++)
		  //      gzprintf(gzout,"%13.6e %12.4e %12.4e %12.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i],E->C[i]);
		  gzprintf(gzout, "%13.6e %12.4e %12.4e %12.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i], E->C[i]);
	    }
	  gzclose(gzout);

	}else{
	  sprintf(output_file, "%s/temp.%d.%d", E->control.data_file2, E->parallel.me, file_number);
	  E->filed[10] = safe_fopen(output_file, "w");
	  fprintf(E->filed[10], "%6d %6d %13.6e\n", E->lmesh.nno, E->advection.timesteps, 
		  E->monitor.elapsed_time);
	  if(file_number % (2 * E->control.record_every) == 0)
	    {
	      if(!E->control.composition)
		for(i = 1; i <= E->lmesh.nno; i++)
		  //        fprintf(E->filed[10],"%13.6e %12.4e %12.4e %13.6e %13.6e %13.6e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i],E->V[1][i],E->V[2][i],E->V[3][i]);
		  //       fprintf(E->filed[10],"%13.6e %12.4e %12.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i]);
		  fprintf(E->filed[10], "%13.6e %12.4e %12.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i]);
	      else
		for(i = 1; i <= E->lmesh.nno; i++)
		  //      fprintf(E->filed[10],"%13.6e %12.4e %12.4e %12.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i],E->C[i]);
		  fprintf(E->filed[10], "%13.6e %12.4e %12.4e %12.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i], E->C[i]);
	    }
	  else
	    {
	      if(!E->control.composition)
		for(i = 1; i <= E->lmesh.nno; i++)
		  //      fprintf(E->filed[10],"%13.6e %12.4e %12.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i]);
		  fprintf(E->filed[10], "%13.6e %12.4e %12.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i]);
	      else
		for(i = 1; i <= E->lmesh.nno; i++)
		  //      fprintf(E->filed[10],"%13.6e %12.4e %12.4e %12.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i],E->C[i]);
		  fprintf(E->filed[10], "%13.6e %12.4e %12.4e %12.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i], E->C[i]);
	    }
	  fclose(E->filed[10]);
	}

	/* 
	   
	end old field output

	*/
      }
      /* 
	 
      surface and bottom properties

      */
      if(E->parallel.me_loc[3] == E->parallel.nprocz - 1) /* top */
	{
	  if(!E->control.COMPRESS){
	    sprintf(output_file, "%s/%d/th_t.%d.%d", 
		    E->control.data_file2, file_number, E->parallel.me, file_number);
	    E->filed[11] = safe_fopen(output_file, "w");
	    fprintf(E->filed[11], "%6d %6d %13.6e %13.6e\n", E->lmesh.nsf, E->advection.timesteps, E->monitor.elapsed_time, E->slice.Nut);
	    for(i = 1; i <= E->lmesh.nsf; i++)
	      fprintf(E->filed[11], "%13.6e %13.6e %13.6e %13.6e\n", 
		      E->slice.tpg[i], E->slice.shflux[i], E->Fas410_b[i], E->Fas670_b[i]);
	    fclose(E->filed[11]);
	  }else{
	    sprintf(output_file, "%s/%d/th_t.%d.%d.gz", 
		    E->control.data_file2, file_number, E->parallel.me, file_number);
	    gzout = safe_gzopen(output_file, "w");
	    gzprintf(gzout, "%6d %6d %13.6e %13.6e\n", E->lmesh.nsf, E->advection.timesteps, E->monitor.elapsed_time, E->slice.Nut);
	    for(i = 1; i <= E->lmesh.nsf; i++)
	      gzprintf(gzout, "%13.6e %13.6e %13.6e %13.6e\n", 
		      E->slice.tpg[i], E->slice.shflux[i], E->Fas410_b[i], E->Fas670_b[i]);
	    gzclose(gzout);
	  }
	}
      if(E->parallel.me_loc[3] == 0) /* bottom */
	{
	  if(!E->control.COMPRESS){
	    sprintf(output_file, "%s/%d/th_b.%d.%d", 
		    E->control.data_file2, file_number, E->parallel.me, file_number);
	    E->filed[11] = safe_fopen(output_file, "w");
	    fprintf(E->filed[11], "%6d %6d %13.6e %13.6e\n", E->lmesh.nsf, E->advection.timesteps, E->monitor.elapsed_time, E->slice.Nub);
	    for(i = 1; i <= E->lmesh.nsf; i++)
	      fprintf(E->filed[11], "%13.6e %13.6e %13.6e %13.6e\n", 
		      E->slice.tpgb[i], E->slice.bhflux[i], E->Fas410_b[i], E->Fas670_b[i]);
	    fclose(E->filed[11]);
	  }else{
	    sprintf(output_file, "%s/%d/th_b.%d.%d.gz", 
		    E->control.data_file2, file_number, E->parallel.me, file_number);
	    gzout = safe_gzopen(output_file, "w");
	    gzprintf(gzout, "%6d %6d %13.6e %13.6e\n", 
		     E->lmesh.nsf, E->advection.timesteps, E->monitor.elapsed_time, E->slice.Nub);
	    for(i = 1; i <= E->lmesh.nsf; i++)
	      gzprintf(gzout, "%13.6e %13.6e %13.6e %13.6e\n", E->slice.tpgb[i], E->slice.bhflux[i], E->Fas410_b[i], E->Fas670_b[i]);
	    gzclose(gzout);
	  }
	}

      /* end output number control */
    }

  /* 

  tracer output only every 10th step

  */
  if((0) && (E->control.composition) && ((file_number % (10*E->control.record_every)) == 0))
    {
      /* output of markers */

      if(E->control.COMPRESS){	/* compressed output */
	sprintf(output_file, "%s/%d/traces.%d.gz", E->control.data_file2, file_number,E->parallel.me);
	gzout = safe_gzopen(output_file, "w");
	gzprintf(gzout, "%6d %6d %13.6e\n", E->advection.markers, E->advection.timesteps, E->monitor.elapsed_time);
	for(i = 1; i <= E->advection.markers; i++)
	  gzprintf(gzout, "%g %g %g %d %d\n", E->XMC[1][i], E->XMC[2][i], E->XMC[3][i], E->CElement[i], E->C12[i]);
	for(i = 1; i <= E->lmesh.nel; i++)
	  gzprintf(gzout, "%g\n", E->CE[i]);
	gzclose(gzout);
      }else{
	sprintf(output_file, "%s/%d/traces.%d", E->control.data_file2, file_number, E->parallel.me);
	E->filed[10] = safe_fopen(output_file, "w");
	fprintf(E->filed[10], "%6d %6d %13.6e\n", E->advection.markers, E->advection.timesteps, E->monitor.elapsed_time);
	for(i = 1; i <= E->advection.markers; i++)
	  fprintf(E->filed[10], "%g %g %g %d %d\n", E->XMC[1][i], E->XMC[2][i], E->XMC[3][i], E->CElement[i], E->C12[i]);
	for(i = 1; i <= E->lmesh.nel; i++)
	  fprintf(E->filed[10], "%g\n", E->CE[i]);
	fclose(E->filed[10]);
      }
    }

  
  return;
}

/* safe gzopen function */
gzFile safe_gzopen(char *name,char *mode)
{
  gzFile tmp;
  char m2[300];
  if((tmp=gzopen(name,mode))==NULL){	
    fprintf(stderr,"error: cannot gzopen file %s, exiting\n",
	    name);
    parallel_process_termination();
  }
  return (tmp);
}	

/* 

   restart 

   gzdir version, will not limit temperatures to [0;1] 
*/
void process_restart_tc_gzdir(struct All_variables *E)
{

  int node, i, j, k, p;
  FILE *fp;
  float temp1, temp2, *temp,cvec[3],pvec[3];
  char input_s[200], input_file[255];
  gzFile gzin;

  if(E->control.Rsphere && (!E->sphere.vtk_base_init)){ 
    /* need to init spherical/cartesian conversion vectors */
    E->sphere.vtk_base = (float *)safe_malloc(sizeof(float)*E->lmesh.nno * 9);
    for(k=0,i=1;i <= E->lmesh.nno;i++,k+=9)
      calc_cbase_at_tp(E->SX[1][i],E->SX[2][i],(E->sphere.vtk_base+k));
    E->sphere.vtk_base_init = 1;
  }
  
  temp = (float *)malloc((E->mesh.noz + 1) * sizeof(float));
  
  if(E->control.restart == 1 || E->control.restart == 3)
    {
      if(E->parallel.me == 0)
	fprintf(stderr,"restarting from %s timestep %i time %g\n",E->convection.old_T_file,E->monitor.solution_cycles,E->monitor.elapsed_time);
  
      /* 
	 temperatures 
      */
      sprintf(input_file,"%s/%d/t.%d.%d.gz",
	      E->convection.old_T_file,  
	      E->monitor.solution_cycles,E->parallel.me, E->monitor.solution_cycles);
      gzin = safe_gzopen(input_file, "r");
      if(gzgets (gzin,input_s, 200) == Z_NULL){
	fprintf(stderr,"read error\n");
	parallel_process_termination();
      }
      sscanf(input_s, "%i %i %f", &i, &j, &E->monitor.elapsed_time);
      if(j != E->lmesh.nno)
	myerror("mismatch of total node number upon restart",E);
      gzgets (gzin,input_s, 200);
      sscanf(input_s, "%i %i", &i, &j);
      /*  */
      for(node = 1; node <= E->lmesh.nno; node++)
	{
	  gzgets (gzin,input_s, 200);
	  sscanf(input_s, "%f", &E->T[node]);
	  //if(E->parallel.me==0)fprintf(stderr,"rt: %6i %11g\n",node,E->T[node]);
	  E->C[node] = 0;
	}
      gzclose(gzin);
      if(E->parallel.me == 0)fprintf(stderr,"restart T OK\n");
      /* 
	 
      velocities

      */
      sprintf(input_file,"%s/%d/vtk_v.%d.%d.gz",
	      E->convection.old_T_file,  
	      E->monitor.solution_cycles,E->parallel.me, E->monitor.solution_cycles);
      gzin = safe_gzopen(input_file, "r");
      if(E->control.Rsphere){
	for(node = 1,k=0; node <= E->lmesh.nno; node++,k+=9)
	  {
	    gzgets (gzin,input_s, 200);
	    sscanf(input_s, "%f %f %f", (cvec),(cvec+1),(cvec+2));
	    convert_cvec_to_pvec(cvec[0],cvec[1],cvec[2],(E->sphere.vtk_base+k),pvec);
	    E->V[3][i] = pvec[0];E->V[1][i]=pvec[1];E->V[2][i]=pvec[2];
	  }
      }else{			/* cartesian */
	for(i=1;i<=E->lmesh.nno;i++) {
	   gzgets (gzin,input_s, 200);
	   sscanf(input_s, "%f %f %f",&(E->V[1][i]),&(E->V[2][i]),&(E->V[3][i]));
	   //if(E->parallel.me==0)fprintf(stderr,"rv: %6i %11g %11g %11g\n",node,E->V[1][i],E->V[2][i],E->V[3][i]);
	}
      }
      gzclose(gzin);
      if(E->parallel.me == 0)fprintf(stderr,"restart V OK\n");
      if(E->control.vtk_pressure_out){
	/* pressures */
	sprintf(input_file,"%s/%d/p.%d.%d.gz",
		E->convection.old_T_file,  
		E->monitor.solution_cycles,E->parallel.me, E->monitor.solution_cycles);
	gzin = safe_gzopen(input_file, "r");
	if(gzgets (gzin,input_s, 200) == Z_NULL){
	  fprintf(stderr,"read error\n");
	  parallel_process_termination();
	}
	sscanf(input_s, "%i %i %f", &i, &j, &E->monitor.elapsed_time);
	if(j != E->lmesh.nno)
	  myerror("mismatch of total node number upon restart",E);
	gzgets (gzin,input_s, 200);
	sscanf(input_s, "%i %i", &i, &j);
	/*  */
	for(node = 1; node <= E->lmesh.nno; node++)
	  {
	    gzgets (gzin,input_s, 200);
	    sscanf(input_s, "%f", &E->NP[node]);
	  }
	gzclose(gzin);
	if(E->parallel.me == 0)fprintf(stderr,"restart P OK\n");
      }
    }
  else 
    {
      myerror("error: restart mode -1 or 2 not implemented for gzdir/vtkout\n",E);
    }
  
  // for composition
  
  if(E->control.composition && (E->control.restart == 1 || E->control.restart == 2)){
    convection_initial_markers(E,0);
    if(E->parallel.me==0)
      fprintf(stderr,"WARNING: not reusing composition, flavor, or strain information for restart\n");
  }else if(E->control.composition && (E->control.restart == 3)){ 
    /* composition */
    sprintf(input_file,"%s/%d/c.%d.%d.gz",
	    E->convection.old_T_file, E->monitor.solution_cycles,
	    E->parallel.me, E->monitor.solution_cycles);
    gzin = safe_gzopen(input_file, "r");
    gzgets (gzin,input_s, 200);
    sscanf(input_s, "%d %d ", &i, &j);
    for(node = 1; node <= E->lmesh.nno; node++){  
	gzgets (gzin,input_s, 200);
	sscanf(input_s, "%f", &E->C[node]);
    }
    gzclose(gzin);

    if(E->tracers_add_flavors)
      fprintf(stderr,"WARNING: not restoring flavor information!!!\n");

    if(E->tracers_track_strain){ 
      /* reinit the strain information */
      sprintf(input_file,"%s/%d/strain.%d.%d.gz",
	      E->convection.old_T_file, E->monitor.solution_cycles,
	      E->parallel.me, E->monitor.solution_cycles);
      gzin = safe_gzopen(input_file, "r");
      gzgets (gzin,input_s, 200);
      sscanf(input_s, "%d %d ", &i, &j);
      for(node = 1; node <= E->lmesh.nno; node++){  
	gzgets (gzin,input_s, 200);
	sscanf(input_s, "%f", &E->strain[node]);
      }
      gzclose(gzin);
    }
    convection_initial_markers1(E);
  }
  
  E->advection.timesteps = E->monitor.solution_cycles;
  

  return;
}


