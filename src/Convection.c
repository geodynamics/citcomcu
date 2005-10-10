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

/* Assumes parameter list is opened and reads the things it needs. 
   Variables are initialized etc, default values are set */


#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include <stdlib.h> /* for "system" command */
#include <strings.h>

void set_convection_defaults(E)
     struct All_variables *E;
{
    void PG_timestep_with_melting();
    void PG_timestep();
    void PG_timestep_particle();
    void PG_process();
    void read_convection_settings();
    void convection_derived_values();
    void convection_allocate_memory();
    void convection_boundary_conditions();
    void node_locations();
    void convection_initial_fields();
    void twiddle_thumbs();
 
    input_int("composition",&(E->control.composition),"0");

    if (E->control.composition)
       E->next_buoyancy_field = PG_timestep_particle;
    else
       E->next_buoyancy_field = PG_timestep;

    E->special_process_new_velocity = PG_process; 
    E->special_process_new_buoyancy = twiddle_thumbs; 
    E->problem_settings = read_convection_settings;
    E->problem_derived_values = convection_derived_values;
    E->problem_allocate_vars = convection_allocate_memory;
    E->problem_boundary_conds = convection_boundary_conditions;
    E->problem_initial_fields = convection_initial_fields;
    E->problem_node_positions = node_locations;
    E->problem_update_node_positions = twiddle_thumbs;
    E->problem_update_bcs = twiddle_thumbs;

    sprintf(E->control.which_data_files,"Temp,Strf,Pres");
    sprintf(E->control.which_horiz_averages,"Temp,Visc,Vrms");
    sprintf(E->control.which_running_data,"Step,Time,");
    sprintf(E->control.which_observable_data,"Shfl");
 
return;
}

void read_convection_settings(E)
     struct All_variables *E;
    
{ 
    void advection_diffusion_parameters();
    char tmp_string[100],tmp1_string[100];
    
/* parameters */

    input_float("rayleigh",&(E->control.Atemp),"essential");

    input_float("rayleigh_comp",&(E->control.Acomp),"essential");
    
    input_boolean("halfspace",&(E->convection.half_space_cooling),"off");
    input_float("halfspage",&(E->convection.half_space_age),"nodefault");
    
    input_int("temperature_blobs",&(E->convection.temp_blobs),"0");
    input_float_vector("temperature_blobx",E->convection.temp_blobs,E->convection.temp_blob_x);
    input_float_vector("temperature_bloby",E->convection.temp_blobs,E->convection.temp_blob_y);
    input_float_vector("temperature_blobz",E->convection.temp_blobs,E->convection.temp_blob_z);
    input_float_vector("temperature_blobsize",E->convection.temp_blobs,E->convection.temp_blob_radius);
    input_float_vector("temperature_blobDT",E->convection.temp_blobs,E->convection.temp_blob_T);
    input_float_vector("temperature_blobbg",E->convection.temp_blobs,E->convection.temp_blob_bg);
    input_int_vector("temperature_blobsticky",E->convection.temp_blobs,E->convection.temp_blob_sticky);
    
    input_int("temperature_zones",&(E->convection.temp_zones),"0");
    input_float_vector("temperature_zonex1",E->convection.temp_zones,E->convection.temp_zonex1);
    input_float_vector("temperature_zonex2",E->convection.temp_zones,E->convection.temp_zonex2);
    input_float_vector("temperature_zonez1",E->convection.temp_zones,E->convection.temp_zonez1);
    input_float_vector("temperature_zonez2",E->convection.temp_zones,E->convection.temp_zonez2);
    input_float_vector("temperature_zoney1",E->convection.temp_zones,E->convection.temp_zoney1);
    input_float_vector("temperature_zoney2",E->convection.temp_zones,E->convection.temp_zoney2);
    input_float_vector("temperature_zoney2",E->convection.temp_zones,E->convection.temp_zoney2);
    input_float_vector("temperature_zoney2",E->convection.temp_zones,E->convection.temp_zoney2);
    input_float_vector("temperature_zonehw",E->convection.temp_zones,E->convection.temp_zonehw);
    input_float_vector("temperature_zonemag",E->convection.temp_zones,E->convection.temp_zonemag);
    input_int_vector("temperature_zonesticky",E->convection.temp_zones,E->convection.temp_zone_sticky);
    
    input_int("num_perturbations",&(E->convection.number_of_perturbations),"0,0,32");
    input_float_vector("perturbmag",E->convection.number_of_perturbations,E->convection.perturb_mag);
    input_float_vector("perturbk",E->convection.number_of_perturbations,E->convection.perturb_k);
    input_float_vector("perturbm",E->convection.number_of_perturbations,E->convection.perturb_mm);
    input_float_vector("perturbl",E->convection.number_of_perturbations,E->convection.perturb_ll);
    
    input_string("prevT",E->convection.old_T_file,"initialize");

  if (E->control.restart)  {
       input_int("restart_timesteps",&(E->monitor.solution_cycles),"0");

       input_string("oldfile",tmp1_string,"initialize");
       input_string("use_scratch",tmp_string,"local");
       if ( strcmp(tmp_string,"local") == 0 )
          strcpy( E->convection.old_T_file, tmp1_string ) ;
       else
          sprintf( E->convection.old_T_file, "/scratch_%s/%s/%s",
                  E->parallel.machinename, tmp_string, tmp1_string ) ;

       }


    
    advection_diffusion_parameters(E);
    
    return;
}

/* =================================================================
   Any setup which relates only to the convection stuff goes in here
   ================================================================= */

void convection_derived_values(E)  
     struct All_variables *E;
 
{ 

return;
}

void convection_allocate_memory(E)
     struct All_variables *E;

{ void advection_diffusion_allocate_memory();

  advection_diffusion_allocate_memory(E);

return;
}

/* ============================================ */
    
void convection_initial_fields(E)
     struct All_variables *E;

{ 
    void convection_initial_temperature();
  void parallel_process_termination();
  int i;

    if (E->control.composition)   {
         E->advection.markers = E->advection.markers_per_ele*E->mesh.nel;
         E->advection.markers = E->advection.markers*E->lmesh.volume/E->mesh.volume;
         E->advection.markers_uplimit = E->advection.markers*2;
 fprintf(stderr,"aaaa %d %g %g\n",E->advection.markers,E->lmesh.volume,E->mesh.volume);
       for (i=1;i<=E->mesh.nsd;i++)   {
         E->VO[i] = (float *) malloc((E->advection.markers_uplimit+1)*sizeof(float));
         E->Vpred[i] = (float *) malloc((E->advection.markers_uplimit+1)*sizeof(float));
         E->XMCpred[i] = (double *) malloc((E->advection.markers_uplimit+1)*sizeof(double));
         E->XMC[i] = (double *) malloc((E->advection.markers_uplimit+1)*sizeof(double));
         E->C12 = (int *) malloc((E->advection.markers_uplimit+1)*sizeof(int));
         E->traces_leave = (int *) malloc((E->advection.markers_uplimit+1)*sizeof(int));
         E->CElement = (int *) malloc((E->advection.markers_uplimit+1)*sizeof(int));
         }
       }

    report(E,"convection, initial temperature");
    convection_initial_temperature(E);

  return; }

/* =========================================== */

void convection_boundary_conditions(E)
     struct All_variables *E;

{
    void velocity_boundary_conditions();
    void temperature_boundary_conditions();
    void temperatures_conform_bcs();
 
    velocity_boundary_conditions(E);      /* universal */
    temperature_boundary_conditions(E);

    temperatures_conform_bcs(E);

    return;
}

/* ===============================
   Initialization of fields .....
   =============================== */

void convection_initial_temperature(E)
     struct All_variables *E;
{
    int ll,mm,i,j,k,p,node,ii;
    double temp,temp1,temp2,temp3,base,radius,radius2;
    double drand48();
    FILE *fp;
    void remove_horiz_ave();
    void vcopy();
    void temperatures_conform_bcs();
    void thermal_buoyancy();
    void parallel_process_termination();
    void convection_initial_markers();
    void process_restart_tc();
    double modified_plgndr_a(),x1,y1,z1,con,beta;
    
    int noz2,nfz,in1,in2,in3,instance,nox,noy,noz;
    char input_s[200],output_file[255];
    float weight,para1,plate_velocity,delta_temp,age;


    const int dims=E->mesh.nsd;

    noy=E->lmesh.noy;  
    noz=E->lmesh.noz;  
    nox=E->lmesh.nox;  

  p=0;

  if (E->control.restart==0)   {
    if (E->control.CART3D)  {
      for(i=1;i<=noy;i++)  
       for(j=1;j<=nox;j++) 
          for(k=1;k<=noz;k++)  {
             node=k+(j-1)*noz+(i-1)*noz*nox;
             x1=E->X[1][node];
             z1=E->X[3][node];
             y1=E->X[2][node];
             E->T[node] = 1-z1 + E->convection.perturb_mag[p] *
                           sin(M_PI*(1.0-z1))*
                      cos(E->convection.perturb_k[p] * M_PI*x1)*
                          ((E->mesh.nsd!=3)? 1.0 :
                        cos(E->convection.perturb_k[p] * M_PI*y1));

//             E->T[node] = 1-z1 + E->convection.perturb_mag[p] *
//                           sin(M_PI*(1.0-z1))*
//                      cos(E->convection.perturb_k[p] * M_PI*y1);
       

//      E->T[node] = 1.0-z1 + E->convection.perturb_mag[p]*drand48();

             E->C[node] = 0;

             E->node[node] = E->node[node] | (INTX | INTZ | INTY);
	     }
      }
    else if (E->control.Rsphere)  {
      mm = E->convection.perturb_mm[0];
      ll = E->convection.perturb_ll[0];
      con = E->convection.perturb_mag[0];
      con = (E->mesh.noz-1)/(E->sphere.ro-E->sphere.ri);
      beta = E->sphere.ri/(E->sphere.ri-E->sphere.ro);
      noz2 = (E->mesh.noz-1)/2+1;

      for(i=1;i<=noy;i++)  
       for(j=1;j<=nox;j++) 
          for(k=1;k<=noz;k++)  {
             ii = k + E->lmesh.nzs - 1;
             node=k+(j-1)*noz+(i-1)*noz*nox;
             x1=E->SX[1][node];
             z1=E->SX[3][node];
             y1=E->SX[2][node];
             E->T[node] = 0.0;
             E->T[node] = beta*(1.0-1.0/z1) + E->convection.perturb_mag[p]*drand48();

/*             if (ii==noz2)    {
                  E->T[node] = con*modified_plgndr_a(ll,mm,x1)*cos(mm*y1);
                  }
             if (ii==noz2)   
               E->T[node] =
              0.01*modified_plgndr_a(ll,mm,x1);
 */   

               E->T[node] = beta*(1.0-1.0/z1)+E->convection.perturb_mag[p]*modified_plgndr_a(ll,mm,x1)*
cos(mm*y1)*sin(M_PI*(z1-E->sphere.ri)/(E->sphere.ro-E->sphere.ri));

             E->T[node] = beta*(1.0-1.0/z1) + E->convection.perturb_mag[p] *
               sin(M_PI*(E->sphere.ro-z1)/(E->sphere.ro-E->sphere.ri))*
               cos(E->convection.perturb_k[p] * M_PI*(x1-E->XG1[1])/(E->XG2[1]-E->XG1[1]))*
                ((E->mesh.nsd!=3)? 1.0 :
               cos(E->convection.perturb_k[p] * M_PI*(y1-E->XG1[2])/(E->XG2[2]-E->XG1[2])));


             E->C[node] = 0.0;


             E->node[node] = E->node[node] | (INTX | INTZ | INTY);
        }
      }       // end for if else if of geometry


      if (E->control.composition)
          convection_initial_markers(E);

    }          // end for restart==0

  else if (E->control.restart)  {

    process_restart_tc(E);

    }

   temperatures_conform_bcs(E);

/*
   sprintf(output_file,"%s.XandT.%d",E->control.data_file,E->parallel.me);
   if ( (fp = fopen(output_file,"w")) != NULL) {
     if (E->control.Rsphere)
      for (j=1;j<=E->lmesh.nno;j++)
          fprintf(fp,"X[%06d] = %.4e Z[%06d] = %.4e Y[%06d] = %.4e T[%06d] = %.4e \n",j,E->SX[1][j],j,E->SX[2][j],j,E->SX[3][j],j,E->T[j]);
     else 
      for (j=1;j<=E->lmesh.nno;j++)
          fprintf(fp,"X[%06d] = %.4e Z[%06d] = %.4e Y[%06d] = %.4e T[%06d] = %.4e \n",j,E->X[1][j],j,E->X[2][j],j,E->X[3][j],j,E->T[j]);
      }        

   fclose(fp);

   exit(11); 
*/

    thermal_buoyancy(E);

    return; 
    }

void process_restart_tc(E)
     struct All_variables *E;
{

  int node,i,j,k,p;
  FILE *fp;
  float temp1,temp2,temp3,*temp;
  char input_s[200],output_file[255];
  void get_C_from_markers();
  void convection_initial_markers();
  void convection_initial_markers1();
  void parallel_process_termination();

  temp = (float *)malloc((E->mesh.noz+1)*sizeof(float));

  if (E->control.restart==1 || E->control.restart==3)  {
    sprintf(output_file,"%s.temp.%d.%d",E->convection.old_T_file,E->parallel.me,E->monitor.solution_cycles);
    fp = fopen(output_file,"r");
    fgets(input_s,200,fp);
    sscanf(input_s,"%d %d %g",&i,&j,&E->monitor.elapsed_time);
 
          for(node=1;node<=E->lmesh.nno;node++)  {
             fgets(input_s,200,fp);
             sscanf(input_s,"%g",&E->T[node]);
//             E->T[node] = min(1.0,E->T[node]);
             E->T[node] = max(0.0,E->T[node]);
             E->node[node] = E->node[node] | (INTX | INTZ | INTY);
	     }
    fclose(fp);  
    }
  else if (E->control.restart==-1)  {
    sprintf(output_file,"%s.temp.%d.%d",E->convection.old_T_file,E->parallel.me,E->monitor.solution_cycles);
    fp = fopen(output_file,"r");
    fgets(input_s,200,fp);
    sscanf(input_s,"%d %d %g",&i,&j,&E->monitor.elapsed_time);
 
          for(node=1;node<=E->lmesh.nno;node++)  {
             fgets(input_s,200,fp);
             sscanf(input_s,"%g %g %g %g",&E->T[node],&temp1,&temp2,&E->C[node]);
//             E->T[node] = min(1.0,E->T[node]);
             E->T[node] = max(0.0,E->T[node]);
             E->node[node] = E->node[node] | (INTX | INTZ | INTY);
	     }
    fclose(fp);  
    }
  else if (E->control.restart==2)  {
//    sprintf(output_file,"%s.ave.0.%d",E->convection.old_T_file,E->monitor.solution_cycles);
    sprintf(output_file,"%s.ave",E->control.data_file);
    fprintf(stderr,"%s %d\n",output_file,E->mesh.noz);
    fp = fopen(output_file,"r");
    fgets(input_s,200,fp);
    E->monitor.solution_cycles=0;
   for(node=1;node<=E->mesh.noz;node++)  {
      fgets(input_s,200,fp);
      sscanf(input_s,"%g %g",&temp1,&temp[node]);
      }
    fclose(fp);  
    for(k=1;k<=E->lmesh.noz;k++)  {
      i = k + E->lmesh.nzs - 1;
      E->Have.T[k] = temp[i];
      fprintf(E->fp,"ab %d %g\n",k,E->Have.T[k]);fflush(E->fp);
    }

   p = 0;
    for(i=1;i<=E->lmesh.noy;i++)  
       for(j=1;j<=E->lmesh.nox;j++) 
          for(k=1;k<=E->lmesh.noz;k++)  {
             node=k+(j-1)*E->lmesh.noz+(i-1)*E->lmesh.noz*E->lmesh.nox;
             E->T[node] = E->Have.T[k] +E->convection.perturb_mag[p]*drand48();
             E->node[node] = E->node[node] | (INTX | INTZ | INTY);
	     }
  }

// for composition

  if (E->control.composition && (E->control.restart==1 || E->control.restart==2)){
          convection_initial_markers(E);
  }

  else if (E->control.composition && (E->control.restart==3))  {
    sprintf(output_file,"%s.comp.%d.%d",E->convection.old_T_file,E->parallel.me,E->monitor.solution_cycles);
    fp = fopen(output_file,"r");
    fgets(input_s,200,fp);
    sscanf(input_s,"%d %d %g",&i,&j,&E->monitor.elapsed_time);
 
          for(node=1;node<=E->lmesh.nno;node++)  {
             fgets(input_s,200,fp);
             sscanf(input_s,"%g",&E->C[node]);
	  }
    fclose(fp);  

    convection_initial_markers1(E);
  }

  else if (E->control.composition && E->control.restart==-1)  {

    sprintf(output_file,"%s.traces.%d",E->convection.old_T_file,E->parallel.me);
    fp = fopen(output_file,"r");
    fgets(input_s,200,fp);
    sscanf(input_s,"%d %d %g",&E->advection.markers,&j,&temp1);
    for(i=1;i<=E->advection.markers;i++)  {
             fgets(input_s,200,fp);
             sscanf(input_s,"%lf %lf %lf %d %d",&E->XMC[1][i],&E->XMC[2][i],&E->XMC[3][i],&E->CElement[i],&E->C12[i]);
	     if (E->XMC[3][i]<E->XP[3][1])
		      E->XMC[3][i]=E->XP[3][1];
	     else if (E->XMC[3][i]>E->XP[3][E->lmesh.noz])
		      E->XMC[3][i]=E->XP[3][E->lmesh.noz];
	     if (E->XMC[2][i]<E->XP[2][1])
		      E->XMC[2][i]=E->XP[2][1];
	     else if (E->XMC[2][i]>E->XP[2][E->lmesh.noy])
		      E->XMC[2][i]=E->XP[2][E->lmesh.noy];
	     if (E->XMC[1][i]<E->XP[1][1])
		      E->XMC[1][i]=E->XP[1][1];
	     else if (E->XMC[1][i]>E->XP[1][E->lmesh.nox])
		      E->XMC[1][i]=E->XP[1][E->lmesh.nox];
    }
   for (i=1;i<=E->lmesh.nel;i++)  {
         fgets(input_s,200,fp);
         sscanf(input_s,"%g",&E->CE[i]);
         }
    fclose(fp);

//    get_C_from_markers(E,E->C);

  }

  E->advection.timesteps = E->monitor.solution_cycles;

  return;
}


void convection_initial_markers1(E)
   struct All_variables *E;
   {
   int *element,el,i,j,k,p,node,ii,jj;
   double  x,y,z,r,t,f,dX[4],dx,dr;
   char input_s[100],output_file[255];
   FILE *fp;
  void get_C_from_markers();
  void parallel_process_termination();
  int get_element();
  int in_the_domain();
  double drand48(),myatan();

    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];

    element = (int *)malloc((E->lmesh.nel+1)*sizeof(int));

  for (el=1;el<=E->lmesh.nel;el++)   {
     t = 0;
     for (j=1;j<=ends;j++) {
       t += E->C[E->ien[el].node[j]];
       }
     t = t/ends;
     E->CE[el] = t;
     element[el] = 0;
  }

  if (E->control.CART3D)  {
     node = 0;
     do {
        x = drand48()*(E->XG2[1]-E->XG1[1]);
        y = drand48()*(E->XG2[2]-E->XG1[2]);
        z = drand48()*(E->XG2[3]-E->XG1[3]);

	if ( (x>=E->XP[1][1]&&x<=E->XP[1][E->lmesh.nox])  &&
	     (y>=E->XP[2][1]&&y<=E->XP[2][E->lmesh.noy])  &&
	     (z>=E->XP[3][1]&&z<=E->XP[3][E->lmesh.noz]) )  {
    	  node ++;
          E->XMC[1][node] = x;
          E->XMC[2][node] = y;
          E->XMC[3][node] = z;

          el = get_element(E,E->XMC[1][node],E->XMC[2][node],E->XMC[3][node],dX);
          E->CElement[node] = el;
	  element[el] ++; 
          }
       }  while (node<E->advection.markers);
    }
  else if (E->control.Rsphere)  {
    node = 0;
    do {
       x = (drand48()-0.5)*2.0;
       y = drand48();
//       y = (drand48()-0.5)*2.0;
       z = (drand48()-0.5)*2.0;
       r = sqrt(x*x+y*y+z*z);
       t = acos(z/r);
       f = myatan(y,x);
       if ( (t>=E->XP[1][1]&&t<=E->XP[1][E->lmesh.nox])  &&
	    (f>=E->XP[2][1]&&f<=E->XP[2][E->lmesh.noy])  &&
	    (r>=E->XP[3][1]&&r<=E->XP[3][E->lmesh.noz]) )  {
	  node ++;
          E->XMC[1][node] = t;
          E->XMC[2][node] = f;
          E->XMC[3][node] = r;
          el = get_element(E,E->XMC[1][node],E->XMC[2][node],E->XMC[3][node],dX);
          E->CElement[node] = el;
	  element[el] ++; 
          }
        } while (node<E->advection.markers);
     }

  for (node=1;node<=E->advection.markers;node++)   {
          el = E->CElement[node];
          dx = 1.0/((double)element[el]);
          dr = 0.5*dx;
	  j = (1-E->CE[el])*element[el]; 
	  if (E->CE[el]<0.5)
             E->C12[node] = 0;
	  else 
             E->C12[node] = 1;

/*	  if (E->CE[el]<dx)
             E->C12[node] = 0;
	  else if (j>=1) {
             E->C12[node] = 0;
	     E->CE[el] +=dx;
	  }
	  else {
             E->C12[node] = 1;
	  }
*/

      }

  get_C_from_markers(E,E->C);

  return;
  }

void convection_initial_markers(E)
   struct All_variables *E;
   {
   int el,i,j,k,p,node,ii,jj;
   double  x,y,z,r,t,f,dX[4],dx,dr;
   char input_s[100],output_file[255];
   FILE *fp;
  void get_C_from_markers();
  void parallel_process_termination();
  int get_element();
  int in_the_domain();
  double drand48(),myatan();

  if (E->control.CART3D)  {
     node = 0;
     do {
        x = drand48()*(E->XG2[1]-E->XG1[1]);
        y = drand48()*(E->XG2[2]-E->XG1[2]);
        z = drand48()*(E->XG2[3]-E->XG1[3]);

	if ( (x>=E->XP[1][1]&&x<=E->XP[1][E->lmesh.nox])  &&
	     (y>=E->XP[2][1]&&y<=E->XP[2][E->lmesh.noy])  &&
	     (z>=E->XP[3][1]&&z<=E->XP[3][E->lmesh.noz]) )  {
    	  node ++;
          E->XMC[1][node] = x;
          E->XMC[2][node] = y;
          E->XMC[3][node] = z;

          el = get_element(E,E->XMC[1][node],E->XMC[2][node],E->XMC[3][node],dX);
          E->CElement[node] = el;
          if (E->XMC[3][node]>E->viscosity.zcomp)
             E->C12[node] = 0;
          else
             E->C12[node] = 1;
          }
       }  while (node<E->advection.markers);
    }
  else if (E->control.Rsphere)  {
    node = 0;
    do {
       x = (drand48()-0.5)*2.0;
       y = drand48();
//       y = (drand48()-0.5)*2.0;
       z = (drand48()-0.5)*2.0;
       r = sqrt(x*x+y*y+z*z);
       t = acos(z/r);
       f = myatan(y,x);
       if ( (t>=E->XP[1][1]&&t<=E->XP[1][E->lmesh.nox])  &&
	    (f>=E->XP[2][1]&&f<=E->XP[2][E->lmesh.noy])  &&
	    (r>=E->XP[3][1]&&r<=E->XP[3][E->lmesh.noz]) )  {
	  node ++;
          E->XMC[1][node] = t;
          E->XMC[2][node] = f;
          E->XMC[3][node] = r;
          el = get_element(E,E->XMC[1][node],E->XMC[2][node],E->XMC[3][node],dX);
          E->CElement[node] = el;
          if (r>E->viscosity.zcomp)
             E->C12[node] = 0;
          else
             E->C12[node] = 1;
          }
        } while (node<E->advection.markers);
     }


  get_C_from_markers(E,E->C);

  return;
  }

void setup_plume_problem(E)
     struct All_variables *E;
{
    int i;
    FILE *fp;
    
    int l,noz;
    char output_file[255];
    double temp,temp1,time_scale,velo_scale;

    noz=E->lmesh.noz;  


    E->control.plate_vel = E->control.plate_vel*0.01/3.1536e7;  
           /* now  with a unit of m/sec */
    velo_scale = E->data.therm_diff/(E->monitor.length_scale);
            /* in m/sec */ 
    E->control.plate_vel = E->control.plate_vel/velo_scale;  
            /* dimensionless */

    time_scale = E->monitor.length_scale*E->monitor.length_scale/E->data.therm_diff;
            /* in sec */ 

    E->control.plate_age = E->control.plate_age*1.0e6*3.1536e7;
            /* in sec */ 
    E->control.plate_age = E->control.plate_age/time_scale;
            /* dimensionless */


    for(i=1;i<=noz;i++)    {
       temp = 0.5*E->X[3][i]/sqrt(E->control.plate_age);
       E->segment.Tz[i] = erf(temp);
       }

    temp = 0.0;
    for(i=1;i<=noz;i++)    {
       l = E->mat[i];
       if (E->viscosity.RHEOL==1)
         E->Have.Vi[i] = E->viscosity.N0[l-1]*exp(E->viscosity.E[l-1]/(E->viscosity.T[l-1]+E->segment.Tz[i]));
       else if (E->viscosity.RHEOL==2)
         E->Have.Vi[i] = E->viscosity.N0[l-1]*exp((E->viscosity.E[l-1]+(1-E->X[3][i])*E->viscosity.Z[l-1])/(E->viscosity.T[l-1]+E->segment.Tz[i]));
       else if (E->viscosity.RHEOL==3)
         E->Have.Vi[i] = E->viscosity.N0[l-1]*exp((E->viscosity.E[l-1]+(1-E->X[3][i])*E->viscosity.Z[l-1])/E->segment.Tz[i]-(E->viscosity.E[l-1]+E->viscosity.Z[l-1]));
      if (i>1) {
        temp += (1/E->Have.Vi[i] + 1/E->Have.Vi[i-1])*(E->X[3][i]-E->X[3][i-1])*0.5;
        }
      }
    temp = E->control.plate_vel/temp;

    temp1 = 0.0;
    for(i=1;i<=noz;i++)    {
       if (i>1) {
         temp1 += (1/E->Have.Vi[i] + 1/E->Have.Vi[i-1])*(E->X[3][i]-E->X[3][i-1])*0.5;
        }
       E->segment.Vx[i] = -temp*temp1 + E->control.plate_vel;
      fprintf(E->fp,"visc %d %g %g %g\n",i,E->Have.Vi[i],E->segment.Tz[i],E->segment.Vx[i]);
       }

 /* non-dimensionalize other parameters */
    E->segment.plume_radius /= E->monitor.length_scale;
    E->segment.plume_DT /= E->data.ref_temperature;
    E->segment.plume_coord[1] /= E->monitor.length_scale;
    if (E->mesh.nsd==3)
      E->segment.plume_coord[3] /= E->monitor.length_scale;

  return;
  }



void PG_process(E,ii)
     struct All_variables *E;
     int ii;
{ 
    void remove_horiz_ave();
    void strain_rate_2_inv();
    void generic_data_storage();
    float return_bulk_value();
    
    float *P,*P2;
    float visc[9];

    int i,j,k,p,a1,nint,n,el;
    int node,hnode;        
    int this_node,that_node;

    struct Shape_function GN;
    struct Shape_function_dx GNx;
    struct Shape_function_dA dOmega;
    
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];

    float x1,x2,int1,int2;

    P=(float *)malloc((1+E->lmesh.nno)*sizeof(double));
    P2=(float *)malloc((1+E->lmesh.nel)*vpts*sizeof(double));
 
    /* This is an ideal point to calculate the energy consistency integral that
       Slava suggested to me.. NOTE, the buoyancy and velocity are most closely in step
       at this point, before advection. */

    strain_rate_2_inv(E,P2,0);  /* strain rate invariant squared */
    v_to_nodes(E,P2,P,E->mesh.levmax);
    
    int1=0.0;
    int2=0.0;

    for(el=1;el<=E->lmesh.nel;el++) {
	
	for(i=1;i<=ends;i++) {
	    n = E->ien[el].node[i];
	    for(j=1;j<=vpts;j++) {
		int1 += E->EVi[(el-1)*vpts+j] * P[n] * E->N.vpt[GNVINDEX(i,j)] * E->gDA[el].vpt[j];
		int2 += E->V[3][n]*E->buoyancy[n] * E->N.vpt[GNVINDEX(i,j)] * E->gDA[el].vpt[j];
	    }
	}
    }

    
    if(int1 != 0.0)
	fprintf(E->fp,"Energy balance integral accurate to  %g%% in %e\n",100.0*fabs((int1-int2)/int1),int1);
  
    
    free((void *)P);
    free((void *)P2);
    return;
}



