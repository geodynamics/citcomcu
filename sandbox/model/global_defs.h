#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

struct MESH_DATA
{
    int nsd;                            /* number of spatial dimensions */
    int nno;
    int nel;

    int nox;                            /* nodes in x-direction */
    int noy;                            /* nodes in y-direction */
    int noz;                            /* nodes in z-direction */

    int elx;
    int ely;
    int elz;

    int nsf;                            /* number of surface nodes */

    float layer[4];                     /* dimensionless dimensions */
};

struct All_variables
{
    /* MPI stuff */
    int  mpi_size;
    int  mpi_rank;
    char mpi_name[MPI_MAX_PROCESSOR_NAME];
    int  mpi_namelen;
    /* parallel stuff */
    int nproc;
    int nprocx;
    int nprocy;
    int nprocz;
    int me_loc[4];
    /* mesh data stuff */
    struct MESH_DATA mesh;
    struct MESH_DATA lmesh;
    /* monitor stuff */
    float cpu_start;
    float cpu_time_elapsed;
    float cpu_time_overhead;
    float elapsed_time;
    int solution_cycles;
    /* control stuff */
    int keep_going;
    int record_every;
    int record_all_until;
    int verbose;
    char data_file[100];
    /* advection ... */
    int min_timesteps;
    int max_timesteps;
    int total_timesteps;
    float timestep;
    /* actual data */
    double *X[4];
    double *V[4];
    double *T;
    /* debugging */
    FILE *fp;
    FILE *fp_info;
    FILE *fp_time;
};

#ifndef __CPROTO__
#include "prototypes.h"
#endif
