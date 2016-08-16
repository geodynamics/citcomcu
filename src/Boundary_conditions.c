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
	int node,top;
	const int bottom = 1;	/* bottom node number */

	for(lv = E->mesh.levmax; lv >= E->mesh.levmin; lv--)
	{
	  top = E->mesh.NOZ[lv]; /* top for this MG level */
		if(E->mesh.botvbc != 1) /* free slip on bottom */
		{
			horizontal_bc(E, E->VB, bottom, 1, 0.0, VBX, 0, lv);
			horizontal_bc(E, E->VB, bottom, 3, 0.0, VBZ, 1, lv);
			horizontal_bc(E, E->VB, bottom, 2, 0.0, VBY, 0, lv);
			horizontal_bc(E, E->VB, bottom, 1, E->control.VBXbotval, SBX, 1, lv);
			horizontal_bc(E, E->VB, bottom, 3, 0.0, SBZ, 0, lv);
			horizontal_bc(E, E->VB, bottom, 2, E->control.VBYbotval, SBY, 1, lv);
		}
		if(E->mesh.topvbc != 1) /* free slip on top */
		{
			horizontal_bc(E, E->VB, top, 1, 0.0, VBX, 0, lv);
			horizontal_bc(E, E->VB, top, 3, 0.0, VBZ, 1, lv);
			horizontal_bc(E, E->VB, top, 2, 0.0, VBY, 0, lv);
			horizontal_bc(E, E->VB, top, 1, E->control.VBXtopval, SBX, 1, lv);
			horizontal_bc(E, E->VB, top, 3, 0.0, SBZ, 0, lv);
			horizontal_bc(E, E->VB, top, 2, E->control.VBYtopval, SBY, 1, lv);
		}
	}

	if(E->mesh.periodic_x || E->mesh.periodic_y){

		velocity_apply_periodic_bcs(E);
	
	}else
		velocity_refl_vert_bc(E);	/* default */

	for(lv = E->mesh.levmax; lv >= E->mesh.levmin; lv--)
	{
	  top = E->mesh.NOZ[lv]; /* top node for this MG level */
		if(E->mesh.botvbc == 1) /* bottom velocity boundary condition */
		{
			horizontal_bc(E, E->VB, bottom, 1, E->control.VBXbotval, VBX, 1, lv);
			horizontal_bc(E, E->VB, bottom, 3, 0.0, VBZ, 1, lv);
			horizontal_bc(E, E->VB, bottom, 2, E->control.VBYbotval, VBY, 1, lv);
			horizontal_bc(E, E->VB, bottom, 1, 0.0, SBX, 0, lv);
			horizontal_bc(E, E->VB, bottom, 3, 0.0, SBZ, 0, lv);
			horizontal_bc(E, E->VB, bottom, 2, 0.0, SBY, 0, lv);
		}
		if(E->mesh.topvbc == 1) /* top velocity boundary condition */
		{
		        E->control.VBXtopval = E->control.plate_vel; /* this should be more explicit */
			horizontal_bc(E, E->VB, top, 1, E->control.VBXtopval, VBX, 1, lv);
			horizontal_bc(E, E->VB, top, 3, 0.0, VBZ, 1, lv);
			horizontal_bc(E, E->VB, top, 2, E->control.VBYtopval, VBY, 1, lv);
			horizontal_bc(E, E->VB, top, 1, 0.0, SBX, 0, lv);
			horizontal_bc(E, E->VB, top, 3, 0.0, SBZ, 0, lv);
			horizontal_bc(E, E->VB, top, 2, 0.0, SBY, 0, lv);
		}
	}
#ifdef USE_GGRD			/* velocities from grd */
	if(E->control.ggrd.vtop_control)
	  ggrd_read_vtop_from_file(E, (int)E->control.Rsphere);
#endif	
	if(E->mesh.periodic_pin_or_filter == 1){
	  if(E->mesh.periodic_x || E->mesh.periodic_y){
	    
	    
	    /* 
	       pin one node at top in lower left corner
	    */
	    for(lv = E->mesh.levmax; lv >= E->mesh.levmin; lv--){
	      node = E->mesh.NOZ[lv];
	      if(E->mesh.periodic_x){
		E->NODE[lv][node] = E->NODE[lv][node] & (~SBX); /* turn of stress */
		E->NODE[lv][node] = E->NODE[lv][node] | (VBX);/* turn on velocity BC */
		if(lv == E->mesh.levmax){	/* NB */
		  E->VB[1][node] = 0;	/* set to zero */
		}
	      }
	      if(E->mesh.periodic_y){
		E->NODE[lv][node] = E->NODE[lv][node] & (~SBY);
		E->NODE[lv][node] = E->NODE[lv][node] | (VBY);
		if(lv == E->mesh.levmax){	/* NB */
		  E->VB[2][node] = 0;
		}
	      }
	    }
	  }
	}

	if(E->mesh.slab_influx_side_bc) /* slab in-flux on side */
	  velocity_apply_slab_influx_side_bc(E);


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



void freeze_surface(struct All_variables *E)
{
  
  if(E->parallel.me == 0)
    fprintf(stderr,"WARNING: freezing surface boundary condition at time step %i (not working yet)\n",
	    E->monitor.solution_cycles);
  /* 

  WARNING: this is not working yet

  */
  /* no slip on top */
  E->mesh.topvbc = 1;E->control.VBXtopval=E->control.VBYtopval=E->control.plate_vel=0.0;
  velocity_boundary_conditions(E);  

  check_bc_consistency(E);
  construct_id(E);
  //construct_lm(E);
  //construct_sub_element(E);
  parallel_shuffle_ele_and_id(E);
  parallel_communication_routs(E);
  //construct_shape_functions(E);
  //mass_matrix(E);
  velocities_conform_bcs(E, E->U);
 
}
/* ========================================== */

void temperature_boundary_conditions(struct All_variables *E)
{
  int node;
  /* bottom */
  if(E->mesh.bottbc == 3){	/* variable bottom temperature condition */
    /* to left of domain */
    horizontal_bc_range(E, E->TB, 1, 3, E->control.TBCbotval, TBZ, 1, E->mesh.levmax,0,E->control.TBCbotval_side_xapply,0,E->mesh.layer[2]);
    horizontal_bc_range(E, E->TB, 1, 3, E->control.TBCbotval, FBZ, 0, E->mesh.levmax,0,E->control.TBCbotval_side_xapply,0,E->mesh.layer[2]);
    /* to right of domain, override node at boundary */
    horizontal_bc_range(E, E->TB, 1, 3, E->control.TBCbotval_side, TBZ, 1, E->mesh.levmax,E->control.TBCbotval_side_xapply,E->mesh.layer[1],0,E->mesh.layer[2]);
    horizontal_bc_range(E, E->TB, 1, 3, E->control.TBCbotval_side, FBZ, 0, E->mesh.levmax,E->control.TBCbotval_side_xapply,E->mesh.layer[1],0,E->mesh.layer[2]);
  }else if((E->mesh.bottbc == 1)||(E->mesh.bottbc == 2)){	
    /* regular temp boundary condition and ggrd temp override */
    horizontal_bc(E, E->TB, 1, 3, E->control.TBCbotval, TBZ, 1, E->mesh.levmax);
    horizontal_bc(E, E->TB, 1, 3, E->control.TBCbotval, FBZ, 0, E->mesh.levmax);
  }else if(E->mesh.bottbc == 0){ /* const flux */
    horizontal_bc(E, E->TB, 1, 3, E->control.TBCbotval, TBZ, 0, E->mesh.levmax);
    horizontal_bc(E, E->TB, 1, 3, E->control.TBCbotval, FBZ, 1, E->mesh.levmax);
  }else if(E->mesh.bottbc == -1){ /* const flux and fixed T*/
    /* flux on left */
    horizontal_bc_range(E, E->TB, 1, 3, E->control.TBCbotval, TBZ, 0, E->mesh.levmax,0,E->control.TBCbotval_side_xapply,0,E->mesh.layer[2]);
    horizontal_bc_range(E, E->TB, 1, 3, E->control.TBCbotval, FBZ, 1, E->mesh.levmax,0,E->control.TBCbotval_side_xapply,0,E->mesh.layer[2]);
    /* temp on right */
    horizontal_bc_range(E, E->TB, 1, 3, E->control.TBCbotval_side, TBZ, 1, E->mesh.levmax,E->control.TBCbotval_side_xapply,E->mesh.layer[1],0,E->mesh.layer[2]);
    horizontal_bc_range(E, E->TB, 1, 3, E->control.TBCbotval_side, FBZ, 0, E->mesh.levmax,E->control.TBCbotval_side_xapply,E->mesh.layer[1],0,E->mesh.layer[2]);
  }
  
  /* top */
  if(E->mesh.toptbc >= 1)	/* temp */
    {
      horizontal_bc(E, E->TB, E->mesh.noz, 3, E->control.TBCtopval, TBZ, 1, E->mesh.levmax);
      horizontal_bc(E, E->TB, E->mesh.noz, 3, E->control.TBCtopval, FBZ, 0, E->mesh.levmax);
    }
  else if(E->mesh.toptbc == 0) /* flux */
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
void velocity_refl_vert_bc(E)
     struct All_variables *E;
{
  int i,j,ii,jj;
  int node1,node2;
  int level,nox,noy,noz;

  /* except one side with XOZ and y=0, all others are not reflecting BC*/
  /* for two YOZ planes if 3-D, or two OZ side walls for 2-D */
  if (((!E->mesh.slab_influx_side_bc)&&(E->parallel.me_loc[1]==0)) || 
      (E->parallel.me_loc[1]==E->parallel.nprocx-1)){
    for(j=1;j<=E->lmesh.noy;j++)
      for(i=1;i<=E->lmesh.noz;i++)  {
	node1 = i  + (j-1)*E->lmesh.noz*E->lmesh.nox;
	node2 = node1 + (E->lmesh.nox-1)*E->lmesh.noz;
	  
	ii = i + E->lmesh.nzs - 1;
	if (E->parallel.me_loc[1]==0 )  {
	  E->VB[1][node1] = 0.0;
	  if((ii != 1) && (ii != E->mesh.noz))
	    E->VB[3][node1] = 0.0;  
	}
	if (E->parallel.me_loc[1]==E->parallel.nprocx-1)  {
	  E->VB[1][node2] = 0.0;
	  if((ii != 1) && (ii != E->mesh.noz))
	    E->VB[3][node2] = 0.0;
	}
      }      /* end loop for i and j */
  }
  /* for two XOZ planes if 3-D */
	
  if ((E->parallel.me_loc[2]==0) || (E->parallel.me_loc[2]==E->parallel.nprocy-1)){
    for(j=1;j<=E->lmesh.nox;j++)
      for(i=1;i<=E->lmesh.noz;i++)       {
	node1 = i +(j-1)*E->lmesh.noz;
	node2 = node1+(E->lmesh.noy-1)*E->lmesh.noz*E->lmesh.nox;
	ii = i + E->lmesh.nzs - 1;

	if (E->parallel.me_loc[2]==0)  {
	  E->VB[2][node1] = 0.0;
	  if((ii != 1) && (ii != E->mesh.noz))
	    E->VB[3][node1] = 0.0;  
	}
	if (E->parallel.me_loc[2]==E->parallel.nprocy-1)  {
	  E->VB[2][node2] = 0.0;
	  if((ii != 1) && (ii != E->mesh.noz))
	    E->VB[3][node2] = 0.0;  
	}
      }    /* end of loop i & j */
  }
  /* all vbc's apply at all levels  */
  for(level=E->mesh.levmax;level>=E->mesh.levmin;level--) {
    nox = E->lmesh.NOX[level] ;
    noz = E->lmesh.NOZ[level] ;
    noy = E->lmesh.NOY[level] ;
    if (((!E->mesh.slab_influx_side_bc)&&(E->parallel.me_loc[1]==0)) || 
	(E->parallel.me_loc[1]==E->parallel.nprocx-1)){
      for(j=1;j<=noy;j++)
	for(i=1;i<=noz;i++) {
	  node1 = i + (j-1)*noz*nox ;
	  node2 = node1 + (nox-1) * noz ;
	  ii = i + E->lmesh.NZS[level] - 1;
	  if (E->parallel.me_loc[1]==0 )  {
	    E->NODE[level][node1] = E->NODE[level][node1] & (~SBX);
	    E->NODE[level][node1] = E->NODE[level][node1] | (VBX);
	    if((ii!=1) && (ii!=E->mesh.NOZ[level])) {
	      E->NODE[level][node1] = E->NODE[level][node1] & (~VBY);
	      E->NODE[level][node1] = E->NODE[level][node1] | SBY;
	      E->NODE[level][node1] = E->NODE[level][node1] & (~ VBZ);
	      E->NODE[level][node1] = E->NODE[level][node1] | SBZ;    
	    }
	  }
	  if (E->parallel.me_loc[1]==E->parallel.nprocx-1)  {
	    E->NODE[level][node2] = E->NODE[level][node2] & (~SBX);
	    E->NODE[level][node2] = E->NODE[level][node2] | (VBX);
	    if((ii!=1) && (ii!=E->mesh.NOZ[level])) {
	      E->NODE[level][node2] = E->NODE[level][node2] & (~VBY);
	      E->NODE[level][node2] = E->NODE[level][node2] | SBY;
	      E->NODE[level][node2] = E->NODE[level][node2] & (~ VBZ);
	      E->NODE[level][node2] = E->NODE[level][node2] | SBZ;
	    }
	  }
	}   /* end for loop i & j */
    } 
    if ((E->parallel.me_loc[2]==0) || (E->parallel.me_loc[2]==E->parallel.nprocy-1)){
      for(j=1;j<=nox;j++) 
	for(i=1;i<=noz;i++) {
	  node1 = i + (j-1)*noz;
	  node2 = node1+(noy-1)*noz*nox;
	  ii = i + E->lmesh.NZS[level] - 1;
	  jj = j + E->lmesh.NXS[level] - 1;
	  if (E->parallel.me_loc[2]==0)  {
	    E->NODE[level][node1] = E->NODE[level][node1] | VBY;
	    E->NODE[level][node1] = E->NODE[level][node1] & (~SBY);
	    if((ii!= 1) && (ii != E->mesh.NOZ[level]))  {
	      E->NODE[level][node1] = E->NODE[level][node1] & (~VBZ);
	      E->NODE[level][node1] = E->NODE[level][node1] | SBZ;
	    } 
	    if((jj!=1) && (jj!=E->mesh.NOX[level]) && (ii!=1) && (ii!=E->mesh.NOZ[level])){
	      E->NODE[level][node1] = E->NODE[level][node1] & (~VBX);
	      E->NODE[level][node1] = E->NODE[level][node1] | SBX;
	    }
	  }
	  if (E->parallel.me_loc[2]==E->parallel.nprocy-1) {
	    E->NODE[level][node2] = E->NODE[level][node2] | VBY;
	    E->NODE[level][node2] = E->NODE[level][node2] & (~SBY);
	    if((ii!= 1) && (ii != E->mesh.NOZ[level]))  {
	      E->NODE[level][node2] = E->NODE[level][node2] & (~VBZ);
	      E->NODE[level][node2] = E->NODE[level][node2] | SBZ;
	    } 
	    if((jj!=1) && (jj!=E->mesh.NOX[level]) && (ii!=1) && (ii!=E->mesh.NOZ[level])){
	      E->NODE[level][node2] = E->NODE[level][node2] & (~VBX);
	      E->NODE[level][node2] = E->NODE[level][node2] | SBX;
	    }
	  }
	  
	}    /* end for loop i & j  */
    }
  }                   /* end for loop level */

  
  return;
}
void velocity_apply_slab_influx_side_bc(struct All_variables *E)
{
  
  int i,j,node,level,nox,noy,noz;
  float vminus,d1,d2;
  d1 = 1.0 - E->mesh.slab_influx_z2;
  d2 = E->mesh.slab_influx_z2 - E->mesh.slab_influx_z1;
  vminus = -2*(d1+d2/2.)/(1-d1-d2)*E->control.sub_vel;

  if(E->parallel.me == 0)
    fprintf(stderr,"WARNING: using experimental side influx boundary condition, v0: %g v-: %g z1: %g z2: %g\n",
	    E->control.sub_vel,vminus, E->mesh.slab_influx_z1, E->mesh.slab_influx_z2);
  
  if(E->mesh.periodic_x)
    myerror("error, need reflective boundary conditions on X for slab influx",E);
  if(E->control.Rsphere)
    myerror("error, need Cartesian geoemetry for slab influx",E);
  if (E->parallel.me_loc[1]==0){
    /* 
       assign velocity values 
    */
    for(j=1;j<=E->lmesh.noy;j++){
      for(i=1;i<=E->lmesh.noz;i++)  {
        node = i  + (j-1)*E->lmesh.noz*E->lmesh.nox;
	//if(E->X[2][node] <= E->mesh.slab_influx_y1){
	  /* apply VBc for slab */
	if(E->X[3][node] >= E->mesh.slab_influx_z2)
	  E->VB[1][node] = E->control.sub_vel;	/* plate velocity above z2 */
	else if(E->X[3][node] >= E->mesh.slab_influx_z1){
	  E->VB[1][node] = E->control.sub_vel * (E->X[3][node] - E->mesh.slab_influx_z1)/d2;	/* tapered down to z1 */
	}else{
	  E->VB[1][node] = vminus * (E->X[3][node]/E->mesh.slab_influx_z1);
	}
	
	/* y and z directions fixed */
	E->VB[2][node] = 0.0;
	E->VB[3][node] = 0.0;
      }
    }
    /* 
       
    set velocity boundary conditions at all levels on the complete left hand size
    
    */
    for(level=E->mesh.levmax;level>=E->mesh.levmin;level--) {
      nox = E->lmesh.NOX[level] ;
      noz = E->lmesh.NOZ[level] ;
      noy = E->lmesh.NOY[level] ;
      for(j=1;j<=noy;j++)
	for(i=1;i<=noz;i++) {
	  node = i + (j-1)*noz*nox ;
	  /* no slip all around */
	  E->NODE[level][node] = E->NODE[level][node] | (VBX);
	  E->NODE[level][node] = E->NODE[level][node] & (~SBX);

	  E->NODE[level][node] = E->NODE[level][node] | (VBY);
	  E->NODE[level][node] = E->NODE[level][node] & (~SBY);

	  E->NODE[level][node] = E->NODE[level][node] | (VBZ);
	  E->NODE[level][node] = E->NODE[level][node] & (~SBZ);

	}
      }
  } /* end only left-most processors branch */
}

void temperature_refl_vert_bc(struct All_variables *E)
{
	int i, j;
	int node1, node2;

	/* Temps and bc-values  at top level only */


	if((E->parallel.me_loc[1] == 0) || (E->parallel.me_loc[1] == E->parallel.nprocx - 1)) /* left or right */
		for(j = 1; j <= E->lmesh.noy; j++)
			for(i = 1; i <= E->lmesh.noz; i++)
			{
				node1 = i + (j - 1) * E->lmesh.noz * E->lmesh.nox;
				node2 = node1 + (E->lmesh.nox - 1) * E->lmesh.noz;
				if(E->parallel.me_loc[1] == 0) /* left zero flux */
				{
					E->node[node1] = E->node[node1] & (~TBX);
					E->node[node1] = E->node[node1] | FBX;
					E->TB[1][node1] = 0.0;
				}
				if(E->parallel.me_loc[1] == E->parallel.nprocx - 1) /* right zero flux */
				{
					E->node[node2] = E->node[node2] & (~TBX);
					E->node[node2] = E->node[node2] | FBX;
					E->TB[1][node2] = 0.0;
				}
			}					/* end for loop i & j */

	if((E->parallel.me_loc[2] == 0) || (E->parallel.me_loc[2] == E->parallel.nprocy - 1)) /* front or back */
		for(j = 1; j <= E->lmesh.nox; j++)
			for(i = 1; i <= E->lmesh.noz; i++)
			{
				node1 = i + (j - 1) * E->lmesh.noz;
				node2 = node1 + (E->lmesh.noy - 1) * E->lmesh.noz * E->lmesh.nox;
				if(E->parallel.me_loc[2] == 0) /* front zero flux */
				{
					E->node[node1] = E->node[node1] & (~TBY);
					E->node[node1] = E->node[node1] | FBY;
					E->TB[2][node1] = 0.0;
				}
				if(E->parallel.me_loc[2] == E->parallel.nprocy - 1) /* back zero flux */
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



	if(E->parallel.me_loc[3] == 0) /* I think this should be the bottom processor */
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

	if(((ROW == 1) && (E->parallel.me_loc[3] == 0)) || ((ROW == E->mesh.NOZ[level]) && (E->parallel.me_loc[3] == E->parallel.nprocz - 1))) /* bottom or top processor */
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

/* only apply the boundary condition between x1 <= x <= x2 y1 <= y <= y2 */
void horizontal_bc_range(struct All_variables *E, float *BC[], int ROW, int dirn, float value, unsigned int mask, char onoff, int level,
			 float x1, float x2, float y1, float y2)

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

	if(((ROW == 1) && (E->parallel.me_loc[3] == 0)) || ((ROW == E->mesh.NOZ[level]) && (E->parallel.me_loc[3] == E->parallel.nprocz - 1))) /* bottom or top processor */
	{

		/* turn bc marker to zero */
		if(onoff == 0)
			for(j = 1; j <= E->lmesh.NOY[level]; j++)
				for(i = 1; i <= E->lmesh.NOX[level]; i++)
				{
					node = rowl + (i - 1) * E->lmesh.NOZ[level] + (j - 1) * E->lmesh.NOZ[level] * E->lmesh.NOX[level];
					if(in_range_for_bc(x1,x2,y1,y2,E,node))
					  E->NODE[level][node] = E->NODE[level][node] & (~mask);
				}				/* end for loop i & j */

		/* turn bc marker to one */
		else
			for(j = 1; j <= E->lmesh.NOY[level]; j++)
				for(i = 1; i <= E->lmesh.NOX[level]; i++)
				{
					node = rowl + (i - 1) * E->lmesh.NOZ[level] + (j - 1) * E->lmesh.NOZ[level] * E->lmesh.NOX[level];
					if(in_range_for_bc(x1,x2,y1,y2,E,node)){
					  E->NODE[level][node] = E->NODE[level][node] | (mask);
					  if(level == E->mesh.levmax)	/* NB */
					    BC[dirn][node] = value;
					}
				}				/* end for loop i & j */

	}							/* end for if */

	return;
}
/* check if node within range */
int in_range_for_bc(float x1,float x2,float y1,float y2,struct All_variables *E,
		    int node)
{

  if(E->control.CART3D){
    if((E->X[1][node] >= x1) && (E->X[1][node] <= x2) &&
       (E->X[2][node] >= y1) && (E->X[2][node] <= y2))
      return 1;
    else
      return 0;
  }else{
    if((E->SX[1][node] >= x1) && (E->SX[1][node] <= x2) && /* theta */
       (E->SX[2][node] >= y1) && (E->SX[2][node] <= y2))   /* phi */
      return 1;
    else
      return 0;
  }

}


void velocity_apply_periodic_bcs(E)
    struct All_variables *E;
{

  int i,j,ii,jj;
  int node1,node2;
  int level,nox,noy,noz;
  
  fprintf(E->fp,"Periodic boundary conditions\n");

  //if (E->mesh.periodic_y && E->mesh.periodic_x)    {
  // return;
  //}

 if (E->mesh.periodic_y)     {

  /* except one side with XOZ and y=0, all others are not reflecting BC*/
  /* for two YOZ planes if 3-D, or two OZ side walls for 2-D */

  if (E->parallel.me_loc[1]==0 || E->parallel.me_loc[1]==E->parallel.nprocx-1)
    for(j=1;j<=E->lmesh.noy;j++)
      for(i=1;i<=E->lmesh.noz;i++)  {
        node1 = i  + (j-1)*E->lmesh.noz*E->lmesh.nox;
        node2 = node1 + (E->lmesh.nox-1)*E->lmesh.noz;

        ii = i + E->lmesh.nzs - 1;
        if (E->parallel.me_loc[1]==0 )  {
          E->VB[1][node1] = 0.0;
          if((ii != 1) && (ii != E->mesh.noz))
              E->VB[3][node1] = 0.0;
          }
        if (E->parallel.me_loc[1]==E->parallel.nprocx-1)  {
          E->VB[1][node2] = 0.0;
          if((ii != 1) && (ii != E->mesh.noz))
              E->VB[3][node2] = 0.0;
          }
        }      /* end loop for i and j */

   }

 if (E->mesh.periodic_x)     {

  /* for two XOZ planes if 3-D */

    if (E->parallel.me_loc[2]==0 || E->parallel.me_loc[2]==E->parallel.nprocy-1)
      for(j=1;j<=E->lmesh.nox;j++)
        for(i=1;i<=E->lmesh.noz;i++)       {
          node1 = i +(j-1)*E->lmesh.noz;
          node2 = node1+(E->lmesh.noy-1)*E->lmesh.noz*E->lmesh.nox;
          ii = i + E->lmesh.nzs - 1;

          if (E->parallel.me_loc[2]==0)  {
             E->VB[2][node1] = 0.0;
             if((ii != 1) && (ii != E->mesh.noz))
                E->VB[3][node1] = 0.0;
             }
          if (E->parallel.me_loc[2]==E->parallel.nprocy-1)  {
             E->VB[2][node2] = 0.0;
             if((ii != 1) && (ii != E->mesh.noz))
                E->VB[3][node2] = 0.0;
             }
          }    /* end of loop i & j */
   }

  /* all vbc's apply at all levels  */
 for(level=E->mesh.levmax;level>=E->mesh.levmin;level--) {
    nox = E->lmesh.NOX[level] ;
    noz = E->lmesh.NOZ[level] ;
    noy = E->lmesh.NOY[level] ;

  if (E->mesh.periodic_y)     {

    if (E->parallel.me_loc[1]==0 || E->parallel.me_loc[1]==E->parallel.nprocx-1)
      for(j=1;j<=noy;j++)
        for(i=1;i<=noz;i++) {
          node1 = i + (j-1)*noz*nox ;
          node2 = node1 + (nox-1) * noz ;
          ii = i + E->lmesh.NZS[level] - 1;
          if (E->parallel.me_loc[1]==0 )  {
              E->NODE[level][node1] = E->NODE[level][node1] & (~SBX);
              E->NODE[level][node1] = E->NODE[level][node1] | (VBX);
              if((ii!=1) && (ii!=E->mesh.NOZ[level])) {
                  E->NODE[level][node1] = E->NODE[level][node1] & (~VBY);
                  E->NODE[level][node1] = E->NODE[level][node1] | SBY;
                  E->NODE[level][node1] = E->NODE[level][node1] & (~ VBZ);
                  E->NODE[level][node1] = E->NODE[level][node1] | SBZ;
                  }
              }
          if (E->parallel.me_loc[1]==E->parallel.nprocx-1)  {
              E->NODE[level][node2] = E->NODE[level][node2] & (~SBX);
              E->NODE[level][node2] = E->NODE[level][node2] | (VBX);
              if((ii!=1) && (ii!=E->mesh.NOZ[level])) {
                  E->NODE[level][node2] = E->NODE[level][node2] & (~VBY);
                  E->NODE[level][node2] = E->NODE[level][node2] | SBY;
                  E->NODE[level][node2] = E->NODE[level][node2] & (~ VBZ);
                  E->NODE[level][node2] = E->NODE[level][node2] | SBZ;
                  }
              }
          }   /* end for loop i & j */
     }

  if (E->mesh.periodic_x)     {

      if (E->parallel.me_loc[2]==0 || E->parallel.me_loc[2]==E->parallel.nprocy-1)
        for(j=1;j<=nox;j++)
          for(i=1;i<=noz;i++) {
            node1 = i + (j-1)*noz;
            node2 = node1+(noy-1)*noz*nox;
            ii = i + E->lmesh.NZS[level] - 1;
            jj = j + E->lmesh.NXS[level] - 1;
            if (E->parallel.me_loc[2]==0)  {
               E->NODE[level][node1] = E->NODE[level][node1] | VBY;
               E->NODE[level][node1] = E->NODE[level][node1] & (~SBY);
               if((ii!= 1) && (ii != E->mesh.NOZ[level]))  {
                  E->NODE[level][node1] = E->NODE[level][node1] & (~VBZ);
                  E->NODE[level][node1] = E->NODE[level][node1] | SBZ;
                  }
               if((jj!=1) && (jj!=E->mesh.NOX[level]) && (ii!=1) && (ii!=E->mesh.NOZ[level])){
                  E->NODE[level][node1] = E->NODE[level][node1] & (~VBX);
                  E->NODE[level][node1] = E->NODE[level][node1] | SBX;
                  }
               }
            if (E->parallel.me_loc[2]==E->parallel.nprocy-1) {
               E->NODE[level][node2] = E->NODE[level][node2] | VBY;
               E->NODE[level][node2] = E->NODE[level][node2] & (~SBY);
               if((ii!= 1) && (ii != E->mesh.NOZ[level]))  {
                  E->NODE[level][node2] = E->NODE[level][node2] & (~VBZ);
                  E->NODE[level][node2] = E->NODE[level][node2] | SBZ;
                  }
               if((jj!=1) && (jj!=E->mesh.NOX[level]) && (ii!=1) && (ii!=E->mesh.NOZ[level])){
                  E->NODE[level][node2] = E->NODE[level][node2] & (~VBX);
                  E->NODE[level][node2] = E->NODE[level][node2] | SBX;
                  }
               }

            }    /* end for loop i & j  */
     }

  }                   /* end for loop level */

    return;

}

void temperature_apply_periodic_bcs(E)
    struct All_variables *E;
{
  int i,j;
  int node1,node2;

 /* Temps and bc-values  at top level only */
/* fixed temperature at x=0 */

  //if (E->mesh.periodic_x && E->mesh.periodic_y)  {
  //return;
  //}

 if (E->mesh.periodic_y)  {

  if (E->parallel.me_loc[1]==0 || E->parallel.me_loc[1]==E->parallel.nprocx-1)
    for(j=1;j<=E->lmesh.noy;j++)
      for(i=1;i<=E->lmesh.noz;i++) {
        node1 = i  + (j-1)*E->lmesh.noz*E->lmesh.nox;
        node2 = node1 + (E->lmesh.nox-1)*E->lmesh.noz;
        if (E->parallel.me_loc[1]==0 )                   {
          E->node[node1] = E->node[node1] & (~TBX);
          E->node[node1] = E->node[node1] | FBX;
          E->TB[1][node1] = 0.0;
          }
        if (E->parallel.me_loc[1]==E->parallel.nprocx-1)   {
          E->node[node2] = E->node[node2] & (~TBX);
          E->node[node2] = E->node[node2] | FBX;
          E->TB[1][node2] = 0.0;
          }
        }       /* end for loop i & j */
      }

 if (E->mesh.periodic_x)  {
    if (E->parallel.me_loc[2]==0 || E->parallel.me_loc[2]==E->parallel.nprocy-1)
      for(j=1;j<=E->lmesh.nox;j++)
        for(i=1;i<=E->lmesh.noz;i++) {
          node1 = i +(j-1)*E->lmesh.noz;
          node2 = node1 +  (E->lmesh.noy-1)*E->lmesh.noz*E->lmesh.nox;
          if (E->parallel.me_loc[2]==0 ) {
             E->node[node1] = E->node[node1] & (~TBY);
             E->node[node1] = E->node[node1] | FBY;
             E->TB[2][node1] = 0.0;
             }
          if (E->parallel.me_loc[2]==E->parallel.nprocy-1) {
             E->node[node2] = E->node[node2] & (~TBY);
             E->node[node2] = E->node[node2] | FBY;
             E->TB[2][node2] = 0.0;
             }
      }
    }

 fprintf(E->fp,"Periodic temperature boundary conditions\n");

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
  int node;
  
  const unsigned int typex = VBX;
  const unsigned int typez = VBZ;
  const unsigned int typey = VBY;
  
  const int nno = E->lmesh.nno;

  for(node = 1; node <= nno; node++){
    /* if(E->SX[3][node] > .999){ */
    /*   fprintf(stderr,"%g %g %g - %i - %i %i - %g %g\n", */
    /* 	      E->SX[1][node],E->SX[2][node],E->SX[3][node], */
    /* 	      node, */
    /* 	      (E->node[node] & typex), */
    /* 	      (E->node[node] & typey), */
    /* 	      E->VB[1][node],E->VB[2][node]); */
    /* } */

    if(E->node[node] & typex)
      U[E->id[node].doff[1]] = E->VB[1][node];
    if(E->node[node] & typey)
      U[E->id[node].doff[2]] = E->VB[2][node];
    if(E->node[node] & typez){
      U[E->id[node].doff[3]] = E->VB[3][node];
    }
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
