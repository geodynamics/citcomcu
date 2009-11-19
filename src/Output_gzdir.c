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

static gzFile *safe_gzopen(char *,char *);
FILE *safe_fopen(char *,char *);
void *safe_malloc (size_t );
void calc_cbase_at_tp(float , float , float *);
void convert_pvec_to_cvec(float ,float , float , float *,
			  float *);
void rtp2xyz(float r, float , float , float *);


void output_velo_related_gzdir(E, file_number)
     struct All_variables *E;
     int file_number;
{
  int el, els, i, j, k, ii, m, node, fd,snode;
  int nox, noz, noy, nfx, nfz, nfy1, nfy2, size1, size2,offset;
  static int vtkout = 1;

  char output_file[1000];
  float cvec[3],locx[3];
  static float *SV, *EV;
  static int been_here = 0;
  float vs;
  static int vtk_base_init = 0;	/* for spherical */
  static int vtk_pressure_out = 1,
    vtk_viscosity_out = 1;
  int vtk_comp_out;



  const int dims = E->mesh.nsd;
  const int ends = enodes[dims];

  const int nno = E->mesh.nno;

  gzFile gzout;






  if(E->control.composition)
    vtk_comp_out = 1;
  else
    vtk_comp_out = 0;

  /* make a directory */
  if(E->parallel.me == 0){
    fprintf(stderr,"Output_gzdir: processing output\n");


    
    sprintf(output_file,"if [ ! -s %s/%d ];then mkdir -p %s/%d;fi 2> /dev/null",
	    E->control.data_file2,file_number,E->control.data_file2,file_number);
    system(output_file);
    fprintf(stderr,"making directory: %s\n",output_file);

  }
  /* and wait for the other jobs */
  parallel_process_sync();
  

  if(been_here == 0){
    /* 
       switch off composition, if not needed
    */
    if(!E->control.composition)
      vtk_comp_out = 0;

    /* 
       
    initial nodeal coordinate output 
    
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
	  gzprintf(gzout, "%.5e %.5e %.5e\n", E->X[1][i], E->X[2][i], E->X[3][i]);
	}
      else
	for(i = 1; i <= E->lmesh.nno; i++)
	  gzprintf(gzout, "%.5e %.5e %.5e\n", E->SX[1][i], E->SX[2][i], E->SX[3][i]);
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
	  gzprintf(gzout, "%.5e %.5e %.5e\n", E->X[1][node], E->X[2][node], E->X[3][node]);
	else
	  gzprintf(gzout, "%.5e %.5e %.5e\n", E->SX[1][node], E->SX[2][node], E->SX[3][node]);
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
	  fprintf(E->filed[13], "%.5e %.5e %.5e\n", E->X[1][i], E->X[2][i], E->X[3][i]);
      else
	for(i = 1; i <= E->lmesh.nno; i++)
	  fprintf(E->filed[13], "%.5e %.5e %.5e\n", E->SX[1][i], E->SX[2][i], E->SX[3][i]);
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
	    fprintf(E->filed[13], "%.5e %.5e %.5e\n", E->X[1][node], E->X[2][node], E->X[3][node]);
	  else
	    fprintf(E->filed[13], "%.5e %.5e %.5e\n", E->SX[1][node], E->SX[2][node], E->SX[3][node]);
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
	if(E->parallel.me == 0)
	  fprintf(stderr, "converting spherical to cartesian  vtk coords to %s \n",output_file);
	for(i=1;i <= E->lmesh.nno;i++) {
	  rtp2xyz((float)E->SX[3][i],(float)E->SX[1][i],(float)E->SX[2][i],locx);
	  gzprintf(gzout,"%9.6f %9.6f %9.6f\n",locx[0],locx[1],locx[2]);
	}
      }else{			/* cartesian */
	if(E->parallel.me == 0)
	  fprintf(stderr, "writing cartesian  vtk coords to %s \n",output_file);
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
      if(E->parallel.me == 0)
	fprintf(stderr, "writing element connectivity %s \n",output_file);
      /* end vtk branch */
    }
    /* end init loop */
    been_here++;
  }
  
  if(E->parallel.me < E->parallel.nprocz)
    {
      sprintf(output_file, "%s/%d/ave_r.%d.%d", E->control.data_file2, file_number, 
	      E->parallel.me, file_number);
      E->filed[10] = safe_fopen(output_file, "w");
      /* header line */
      fprintf(E->filed[10], "# %6d %6d %.5e %.5e %.5e %.4e %.4e %.5e %.5e\n", 
	      E->lmesh.noz, E->advection.timesteps, E->monitor.elapsed_time, 
	      E->slice.Nut, E->slice.Nub, E->data.T_adi0, E->data.T_adi1, 
	      E->monitor.Sigma_interior, E->monitor.Sigma_max);
      /* 
	 added first column of radial coordinate TWB 
      */
      if(!E->control.Rsphere){	/* cartesian */
	for(i = 1; i <= E->lmesh.noz; i++){
	  fprintf(E->filed[10], "%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n", 
		  E->X[3][i],E->Have.T[i], E->Have.vrms[i], E->Have.Vi[i], E->Have.Rho[i], 
		  E->Have.F[i], E->Have.f[i], E->Have.C[i], E->Have.Tadi[i]);
	}
      }else{			/* spherical */
	for(i = 1; i <= E->lmesh.noz; i++){
	  fprintf(E->filed[10], "%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n", 
		  E->SX[3][i],E->Have.T[i], E->Have.vrms[i], E->Have.Vi[i], E->Have.Rho[i], 
		  E->Have.F[i], E->Have.f[i], E->Have.C[i], E->Have.Tadi[i]);
	}
      }
      fclose(E->filed[10]);
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
	gzprintf(gzout,"%d %d %.5e\n",
		 file_number,E->lmesh.nno,
		 E->monitor.elapsed_time);
	gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno); /* for backward compatibility */
	for(i=1;i<=E->lmesh.nno;i++)           
	  gzprintf(gzout,"%.6e\n",E->T[i]);
	gzclose(gzout);
	/* velocities */
	if(E->control.Rsphere && (!vtk_base_init)){ 
	  /* 

	  STILL NEED TO CHECK THE CONVENTION OF ORDERING FOR SPHERICAL HERE!
	  */

	  /* need to init spherical/cartesian conversion vectors */
	  E->sphere.vtk_base = (float *)safe_malloc(sizeof(float)*E->lmesh.nno*9);
	  for(k=0,i=1;i<=E->lmesh.nno;i++,k+=9)
	    calc_cbase_at_tp(E->SX[1][i],E->SX[2][i],(E->sphere.vtk_base+k));
	  vtk_base_init = 1;
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
	  for(i=1;i<=E->lmesh.nno;i++) 
	    gzprintf(gzout,"%.6e %.6e %.6e\n",E->V[1][i],E->V[2][i],E->V[3][i]);
	}
	gzclose(gzout);
	if(vtk_pressure_out){
	  /* 
	     pressure at nodes 
	  */
	  sprintf(output_file,"%s/%d/p.%d.%d.gz",E->control.data_file2,
		  file_number, E->parallel.me,file_number);
	  gzout = safe_gzopen(output_file,"w");
	  gzprintf(gzout,"%d %d %.5e\n",file_number,E->lmesh.nno,E->monitor.elapsed_time);
	  gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	  for(i=1;i<=E->lmesh.nno;i++)          
	    gzprintf(gzout,"%.6e\n",E->NP[i]);
	  gzclose(gzout);
	}
	if(vtk_viscosity_out){
	  /* 
	     
	     viscosity output
	  
	  */
	  /* write to file */
	  sprintf(output_file,"%s/%d/visc.%d.%d.gz",
		  E->control.data_file2,file_number, E->parallel.me,file_number);
	  gzout=safe_gzopen(output_file,"w");
	  gzprintf(gzout,"%3d %7d\n",1,E->lmesh.nno);
	  for(i=1;i<=E->lmesh.nno;i++){
	    gzprintf(gzout,"%.4e\n",E->VI[E->mesh.levmax][i]);
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
	    gzprintf(gzout,"%.4e\n",E->C[i]);
	  gzclose(gzout);
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
	  gzprintf(gzout, "%6d %6d %.5e\n", E->lmesh.nno, E->advection.timesteps, 
		  E->monitor.elapsed_time);
	  if((file_number % 2* E->control.record_every) == 0)
	    {
	      if(E->control.composition == 0)
		for(i = 1; i <= E->lmesh.nno; i++)
		  //        gzprintf(gzout,"%.5e %.4e %.4e %.5e %.5e %.5e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i],E->V[1][i],E->V[2][i],E->V[3][i]);
		  //       gzprintf(gzout,"%.5e %.4e %.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i]);
		  gzprintf(gzout, "%.5e %.4e %.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i]);
	      else if(E->control.composition)
	    for(i = 1; i <= E->lmesh.nno; i++)
	      //      gzprintf(gzout,"%.5e %.4e %.4e %.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i],E->C[i]);
	      gzprintf(gzout, "%.5e %.4e %.4e %.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i], E->C[i]);
	    }
	  else
	    {
	      if(E->control.composition == 0)
		for(i = 1; i <= E->lmesh.nno; i++)
		  //      gzprintf(gzout,"%.5e %.4e %.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i]);
		  gzprintf(gzout, "%.5e %.4e %.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i]);
	      else if(E->control.composition)
		for(i = 1; i <= E->lmesh.nno; i++)
		  //      gzprintf(gzout,"%.5e %.4e %.4e %.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i],E->C[i]);
		  gzprintf(gzout, "%.5e %.4e %.4e %.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i], E->C[i]);
	    }
	  gzclose(gzout);

	}else{
	  sprintf(output_file, "%s/temp.%d.%d", E->control.data_file2, E->parallel.me, file_number);
	  E->filed[10] = safe_fopen(output_file, "w");
	  fprintf(E->filed[10], "%6d %6d %.5e\n", E->lmesh.nno, E->advection.timesteps, 
		  E->monitor.elapsed_time);
	  if(file_number % (2 * E->control.record_every) == 0)
	    {
	      if(E->control.composition == 0)
		for(i = 1; i <= E->lmesh.nno; i++)
		  //        fprintf(E->filed[10],"%.5e %.4e %.4e %.5e %.5e %.5e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i],E->V[1][i],E->V[2][i],E->V[3][i]);
		  //       fprintf(E->filed[10],"%.5e %.4e %.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i]);
		  fprintf(E->filed[10], "%.5e %.4e %.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i]);
	      else if(E->control.composition)
	    for(i = 1; i <= E->lmesh.nno; i++)
	      //      fprintf(E->filed[10],"%.5e %.4e %.4e %.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i],E->C[i]);
	      fprintf(E->filed[10], "%.5e %.4e %.4e %.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i], E->C[i]);
	    }
	  else
	    {
	      if(E->control.composition == 0)
		for(i = 1; i <= E->lmesh.nno; i++)
		  //      fprintf(E->filed[10],"%.5e %.4e %.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i]);
		  fprintf(E->filed[10], "%.5e %.4e %.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i]);
	      else if(E->control.composition)
		for(i = 1; i <= E->lmesh.nno; i++)
		  //      fprintf(E->filed[10],"%.5e %.4e %.4e %.4e\n",E->T[i],E->heatflux[i],E->heatflux_adv[i],E->C[i]);
		  fprintf(E->filed[10], "%.5e %.4e %.4e %.4e\n", E->T[i], E->V[3][i], E->heatflux_adv[i], E->C[i]);
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
	    fprintf(E->filed[11], "%6d %6d %.5e %.5e\n", E->lmesh.nsf, E->advection.timesteps, E->monitor.elapsed_time, E->slice.Nut);
	    for(i = 1; i <= E->lmesh.nsf; i++)
	      fprintf(E->filed[11], "%.5e %.5e %.5e %.5e\n", 
		      E->slice.tpg[i], E->slice.shflux[i], E->Fas410_b[i], E->Fas670_b[i]);
	    fclose(E->filed[11]);
	  }else{
	    sprintf(output_file, "%s/%d/th_t.%d.%d.gz", 
		    E->control.data_file2, file_number, E->parallel.me, file_number);
	    gzout = safe_gzopen(output_file, "w");
	    gzprintf(gzout, "%6d %6d %.5e %.5e\n", E->lmesh.nsf, E->advection.timesteps, E->monitor.elapsed_time, E->slice.Nut);
	    for(i = 1; i <= E->lmesh.nsf; i++)
	      gzprintf(gzout, "%.5e %.5e %.5e %.5e\n", 
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
	    fprintf(E->filed[11], "%6d %6d %.5e %.5e\n", E->lmesh.nsf, E->advection.timesteps, E->monitor.elapsed_time, E->slice.Nub);
	    for(i = 1; i <= E->lmesh.nsf; i++)
	      fprintf(E->filed[11], "%.5e %.5e %.5e %.5e\n", 
		      E->slice.tpgb[i], E->slice.bhflux[i], E->Fas410_b[i], E->Fas670_b[i]);
	    fclose(E->filed[11]);
	  }else{
	    sprintf(output_file, "%s/%d/th_b.%d.%d.gz", 
		    E->control.data_file2, file_number, E->parallel.me, file_number);
	    gzout = safe_gzopen(output_file, "w");
	    gzprintf(gzout, "%6d %6d %.5e %.5e\n", 
		     E->lmesh.nsf, E->advection.timesteps, E->monitor.elapsed_time, E->slice.Nub);
	    for(i = 1; i <= E->lmesh.nsf; i++)
	      gzprintf(gzout, "%.5e %.5e %.5e %.5e\n", E->slice.tpgb[i], E->slice.bhflux[i], E->Fas410_b[i], E->Fas670_b[i]);
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
	gzprintf(gzout, "%6d %6d %.5e\n", E->advection.markers, E->advection.timesteps, E->monitor.elapsed_time);
	for(i = 1; i <= E->advection.markers; i++)
	  gzprintf(gzout, "%g %g %g %d %d\n", E->XMC[1][i], E->XMC[2][i], E->XMC[3][i], E->CElement[i], E->C12[i]);
	for(i = 1; i <= E->lmesh.nel; i++)
	  fprintf(gzout, "%g\n", E->CE[i]);
	gzclose(gzout);
      }else{
	sprintf(output_file, "%s/%d/traces.%d", E->control.data_file2, file_number, E->parallel.me);
	E->filed[10] = safe_fopen(output_file, "w");
	fprintf(E->filed[10], "%6d %6d %.5e\n", E->advection.markers, E->advection.timesteps, E->monitor.elapsed_time);
	for(i = 1; i <= E->advection.markers; i++)
	  fprintf(E->filed[10], "%g %g %g %d %d\n", E->XMC[1][i], E->XMC[2][i], E->XMC[3][i], E->CElement[i], E->C12[i]);
	for(i = 1; i <= E->lmesh.nel; i++)
	  fprintf(E->filed[10], "%g\n", E->CE[i]);
	fclose(E->filed[10]);
      }
    }

if(E->parallel.me==0)fprintf(stderr,"vel output done\n");
  return;
}

/* safe gzopen function */
static gzFile *safe_gzopen(char *name,char *mode)
{
  gzFile *tmp;
  char m2[300];
  if((tmp=(gzFile *)gzopen(name,mode))==NULL){	
    fprintf(stderr,"error: cannot gzopen file %s, exiting\n",
	    name);
    parallel_process_termination();
  }
  return ((gzFile *)tmp);
}	

/* gzdir version, will not limit temperatures to [0;1] */
void process_restart_tc_gzdir(struct All_variables *E)
{

  int node, i, j, k, p;
  FILE *fp;
  float temp1, temp2, *temp;
  char input_s[200], input_file[255];
  gzFile gzin;

  temp = (float *)malloc((E->mesh.noz + 1) * sizeof(float));
  
  if(E->control.restart == 1 || E->control.restart == 3)
    {

      sprintf(input_file,"%s/%d/t.%d.%d.gz",
	      E->convection.old_T_file,  
	      E->monitor.solution_cycles,E->parallel.me, E->monitor.solution_cycles);
      gzin = safe_gzopen(input_file, "r");
      if(gzgets (gzin,input_s, 200) == Z_NULL){
	fprintf(stderr,"read error\n");
	parallel_process_termination();
      }
      sscanf(input_s, "%i %i %g", &i, &j, &E->monitor.elapsed_time);

      for(node = 1; node <= E->lmesh.nno; node++)
	{
	  gzgets (gzin,input_s, 200);
	  sscanf(input_s, "%g", &E->T[node]);
	  //if(E->SX[3][node] == 0)fprintf(stderr,"%g %g\n",E->SX[3][node],E->T[node]);
	  E->C[node] = 0;
	}
      gzclose(gzin);
      if(E->parallel.me == 0)
	fprintf(stderr,"restarting from %s timestep %i time %g\n",
		E->convection.old_T_file,E->monitor.solution_cycles,
		E->monitor.elapsed_time);
       
    }
  else 
    {
      fprintf(stderr,"error: restart mode -1 or 2 not implemented for gzdir/vtkout\n");
      parallel_process_termination();
    }
  
  // for composition
  
  if(E->control.composition && (E->control.restart == 1 || E->control.restart == 2))
    {
      convection_initial_markers(E,0);
    }
  
  else if(E->control.composition && (E->control.restart == 3))
    { 
      sprintf(input_file,"%s/%d/c.%d.%d.gz",
	      E->convection.old_T_file, E->monitor.solution_cycles,
	      E->parallel.me, E->monitor.solution_cycles);
      gzin = safe_gzopen(input_file, "r");
      gzgets (gzin,input_s, 200);
      sscanf(input_s, "%d %d ", &i, &j);
      
      for(node = 1; node <= E->lmesh.nno; node++)
	{  
	  gzgets (gzin,input_s, 200);
   	  sscanf(input_s, "%g", &E->C[node]);
	}
      gzclose(gzin);
      
      convection_initial_markers1(E);
    }
  
  E->advection.timesteps = E->monitor.solution_cycles;
  

  return;
}


