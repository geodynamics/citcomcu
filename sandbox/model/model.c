
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <hdf5.h>
#include "global_defs.h"

#define malloc0(a)      safe_malloc((a), __FILE__, __LINE__)

#ifdef DEBUG
#   define  report0(a,b)    report((a),(b))
#else
#   define  report0(a,b)
#endif

int Emergency_stop;


int main(int argc, char *argv[])
{
    struct All_variables *E;
    double start_time;

    E = (struct All_variables *)malloc0(sizeof(struct All_variables));

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &(E->mpi_size));
    MPI_Comm_rank(MPI_COMM_WORLD, &(E->mpi_rank));
    MPI_Get_processor_name(E->mpi_name, &(E->mpi_namelen));

    if(argc < 2)
    {
        fprintf(stderr, "Usage: %s PARAMETERFILE\n", argv[0]);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    read_instructions(E, argv[1]);

    /* prepare for main loop */
    Emergency_stop = 0;
    E->cpu_start = MPI_Wtime();

    /* main loop */
    process(E);
    output(E);
    report_cycle(E);
    while(E->keep_going && (Emergency_stop == 0))
    {
        E->solution_cycles++;
        process(E);
        output(E);
        report_cycle(E);
    }
    report_summary(E);

    MPI_Finalize();
    finalize(E);
    free(E);

    return EXIT_SUCCESS;
}

void finalize(struct All_variables *E)
{
    finalize_hdf5(E);
    fclose(E->fp);
    fclose(E->fp_time);
    if(E->verbose)
        fclose(E->fp_info);
}

void read_instructions(struct All_variables *E, char *filename)
{
    double time = MPI_Wtime();

    setup_parser(filename, (E->mpi_rank == 0));

    global_default_values(E);
    read_initial_settings(E);

    open_log(E);
    open_time(E);
    if(E->verbose)
        open_info(E);

    global_derived_values(E);
    
    parallel_processor_setup(E);
    parallel_domain_decomp0(E);

    allocate_vars(E);

    shutdown_parser();

    /* report on overhead */
    if(E->mpi_rank == 0)
    {
        if(E->verbose)
        {
            fprintf(stderr, "Input parameters taken from file %s\n",
                    filename);
            fprintf(stderr, "Initialization complete after %g seconds\n",
                    MPI_Wtime() - time);
        }
        fprintf(E->fp, "Initialization complete after %g seconds\n",
                MPI_Wtime() - time);
        fflush(E->fp);
    }
    report0(E, "read_instructions: ok");

    /* calculate overhead last */
    E->cpu_time_overhead = MPI_Wtime() - time;
}

void global_default_values(struct All_variables *E)
{
    /* initialize some important variables */
    E->mesh.nsd = 3;
    E->lmesh.nsd = 3;

    /* initialize some monitor & control variables */
    E->keep_going = 1;
    E->total_timesteps = 0;
    E->min_timesteps = 1;
    E->max_timesteps = 100;
    E->timestep = 1.0;
    E->elapsed_time = 0.0;

    /* others */
    E->nproc  = E->mpi_size;
    E->nprocx = 1;
    E->nprocy = 1;
    E->nprocz = 1;
    E->me_loc[1] = 0;
    E->me_loc[2] = 0;
    E->me_loc[3] = 0;

    report0(E, "global_default_values: ok");
}

void read_initial_settings(struct All_variables *E)
{
    input_string("datafile", E->data_file, "model");
    input_boolean("verbose", &(E->verbose), "off");
    input_int("storage_spacing", &(E->record_every), "10");
    
    input_int("nprocx", &(E->nprocx), "1");
    input_int("nprocy", &(E->nprocy), "1");
    input_int("nprocz", &(E->nprocz), "1");

    input_float("dimenx", &(E->mesh.layer[1]), "nodefault");
    input_float("dimeny", &(E->mesh.layer[2]), "nodefault");
    input_float("dimenz", &(E->mesh.layer[3]), "nodefault");
    
    input_int("nodex", &(E->mesh.nox), "nodefault,1,nomax");
    input_int("nodey", &(E->mesh.noy), "nodefault,1,nomax");
    input_int("nodez", &(E->mesh.noz), "nodefault,1,nomax");

    /* advection ... */
    input_int("minstep", &(E->min_timesteps), "1");
    input_int("maxstep", &(E->max_timesteps), "1000");
    input_float("fixed_timestep", &(E->timestep), "1.0");

    report0(E, "read_initial_settings: ok");
}


void global_derived_values(struct All_variables *E)
{
    E->mesh.elx = E->mesh.nox - 1;
    E->mesh.ely = E->mesh.noy - 1;
    E->mesh.elz = E->mesh.noz - 1;

    E->mesh.nno = E->mesh.nox * E->mesh.noy * E->mesh.noz;
    E->mesh.nsf = E->mesh.nox * E->mesh.noy;
    E->mesh.nel = E->mesh.elx * E->mesh.ely * E->mesh.elz;

    report0(E, "global_derived_values: ok");
}

void parallel_processor_setup(struct All_variables *E)
{
    int i,j,k;
    int me = E->mpi_rank;
    int nprocx = E->nprocx;
    int nprocy = E->nprocy;
    int nprocz = E->nprocz;

    /* determine the location of processors in each cap:
     * 
     * note that given (i,j,k), the processor rank would be
     *
     *   rank = k + i*nprocz + j*nprocz*nprocx
     * 
     * solving for i,j,k yields the following formulas:
     *
     */
    k = me % nprocz;
    i = ((me - k)/nprocz) % nprocx;
    j = ((((me - k)/nprocz) - i)/nprocx) % nprocy;

    E->me_loc[1] = i;
    E->me_loc[2] = j;
    E->me_loc[3] = k;

    if(E->verbose)
    {
        fprintf(E->fp_info, "me=%d loc1=%d loc2=%d loc3=%d\n", me,
                E->me_loc[1], E->me_loc[2], E->me_loc[3]);
    }

    report0(E, "parallel_processor_setup: ok");
}

/* get element information for each processor */
void parallel_domain_decomp0(struct All_variables *E)
{
    int me = E->mpi_rank;

    E->lmesh.elx = E->mesh.elx / E->nprocx;
    E->lmesh.elz = E->mesh.elz / E->nprocz;
    E->lmesh.ely = E->mesh.ely / E->nprocy;

    E->lmesh.nox = E->lmesh.elx + 1;
    E->lmesh.noz = E->lmesh.elz + 1;
    E->lmesh.noy = E->lmesh.ely + 1;

    E->lmesh.nno = E->lmesh.noz * E->lmesh.nox * E->lmesh.noy;
    E->lmesh.nel = E->lmesh.elz * E->lmesh.elx * E->lmesh.ely;
    E->lmesh.nsf = E->lmesh.nox * E->lmesh.noy;

    report0(E, "parallel_domain_decomp0: ok");
}

void allocate_vars(struct All_variables *E)
{
    int nsd = E->lmesh.nsd;
    int nno = E->lmesh.nno;
    //int nsf = E->lmesh.nsf;
    //int nox = E->lmesh.nox;
    //int noz = E->lmesh.noz;
    //int noy = E->lmesh.noy;
    //int elx = E->lmesh.elx;
    //int ely = E->lmesh.ely;

    int d;
    int n;

    /* allocate */
    for(d=1; d <= nsd; d++)
        E->X[d] = (double *)malloc0((nno+1)*sizeof(double));

    for(d=1; d <= nsd; d++)
        E->V[d] = (double *)malloc0((nno+1)*sizeof(double));

    E->T = (double *)malloc0((nno+1)*sizeof(double));

    /* initialize - (XXX move this elsewhere) */
    for(n = 1; n <= nno; n++)
        for(d = 1; d <= nsd; d++)
            E->X[d][n] = 0.0;

    for(n = 1; n <= nno; n++)
        for(d = 1; d <= nsd; d++)
            E->V[d][n] = 0.0;

    for(n = 1; n <= nno; n++)
        E->T[n] = 0.0;

    report0(E, "allocate_vars: ok");
}


/*
 * Data Processing
 */

void process(struct All_variables *E)
{
    int me = E->mpi_rank;
    int nno = E->lmesh.nno;
    int cycles = E->solution_cycles;

    int n;
    double x = cycles/1000.0;


    /* Data gets processed here.
     * Nothing to be done, yet.
     */
    for(n = 1; n <= nno; n++)
        E->T[n] = me + x;

    for(n = 1; n <= nno; n++)
    {
        E->V[1][n] = me + x;
        E->V[2][n] = me + x/10;
        E->V[3][n] = me + x/100;
    }


    /*
     * Now, some important bookkeeping.
     */
    E->total_timesteps++;
    E->elapsed_time += E->timestep;

    if(E->solution_cycles < E->max_timesteps)
        E->keep_going = 1;
    else
        E->keep_going = 0;
}


/*
 * Reporting Routines
 */

void report_cycle(struct All_variables *E)
{
    static double time;

    if(E->solution_cycles == 0)
        time = MPI_Wtime();

    if(E->mpi_rank == 0)
    {
        /* write to log */
        fprintf(E->fp, "CPU total = %g & CPU = %g for step %d time = %.4e "
                "dt = %.4e\n", (time - E->cpu_start), (MPI_Wtime() - time),
                E->solution_cycles, E->elapsed_time, E->timestep);
        /* write to time output file */
        fprintf(E->fp_time, "%d %.4e %.4e %.4e %.4e\n",
                E->solution_cycles, E->elapsed_time, E->timestep,
                (MPI_Wtime() - E->cpu_start), (MPI_Wtime() - time));
    }

    /* prepare for next call */
    time = MPI_Wtime();
}


void report_summary(struct All_variables *E)
{
    double time = MPI_Wtime();
    if(E->mpi_rank == 0)
    {
        if(E->verbose)
        {
            fprintf(stderr, "cycles=%d\n", E->solution_cycles);
        }
        fprintf(E->fp, "Initialization overhead = %f\n",
                E->cpu_time_overhead);
        fprintf(E->fp, "Average cpu time taken per step = %f\n",
                (time - E->cpu_start)/((float)(E->solution_cycles)));
    }
}



/*
 * Logging Routines
 */
void open_log(struct All_variables *E)
{
    char logfile[255];
    sprintf(logfile, "%s.log", E->data_file);
    E->fp = fopen2(logfile);
}
void open_info(struct All_variables *E)
{
    char infofile[255];
    sprintf(infofile, "%s.info.%d", E->data_file, E->mpi_rank);
    E->fp_info = fopen2(infofile);
}
void open_time(struct All_variables *E)
{
    char timefile[255];
    sprintf(timefile, "%s.time", E->data_file);
    E->fp_time = fopen2(timefile);
}
void report(struct All_variables *E, char *string)
{
    if(E->verbose && E->mpi_rank == 0)
    {
        fprintf(stderr, "%s\n", string);
        fflush(stderr);
    }
}
void record(struct All_variables *E, char *string)
{
    if(E->verbose)
    {
        fprintf(E->fp, "%s\n", string);
        fflush(E->fp);
    }
}



/*
 * Output Routines
 *
 */

void output(struct All_variables *E)
{
    int cycles = E->solution_cycles;
    
    if(cycles % E->record_every != 0)
        return;

    if(cycles == 0)
        output_coord(E);

    output_velocity(E);
    output_temperature(E);
}

void output_coord(struct All_variables *E)
{
    FILE *fp;
    char output_file[255];

    int cycles = E->solution_cycles;
    int me = E->mpi_rank;

    int nno = E->lmesh.nno;
    int n;

    sprintf(output_file, "%s.coord.%d", E->data_file, me);
    fp = fopen2(output_file);
    
    fprintf(fp, "%6d\n", nno);
    for(n = 1; n <= nno; n++)
        fprintf(fp, "%.6e %.6e %.6e\n", E->X[1][n], E->X[2][n], E->X[3][n]);
    fclose(fp);
}

void output_velocity(struct All_variables *E)
{
    FILE *fp;
    char output_file[255];

    int me = E->mpi_rank;
    int cycles = E->solution_cycles;
    
    int nno = E->lmesh.nno;
    int n;
    
    sprintf(output_file, "%s.velo.%d.%d", E->data_file, me, cycles);
    fp = fopen2(output_file);

    fprintf(fp, "%d %d %.6e\n", cycles, nno, E->elapsed_time);
    for(n = 1; n <= nno; n++)
        fprintf(fp, "%.6e %.6e %.6e\n", E->V[1][n], E->V[2][n], E->V[3][n]);
    fclose(fp);
}

void output_temperature(struct All_variables *E)
{
    FILE *fp;
    char output_file[255];

    int me = E->mpi_rank;
    int cycles = E->solution_cycles;

    int nno = E->lmesh.nno;
    int n;

    sprintf(output_file, "%s.temp.%d.%d", E->data_file, me, cycles);
    fp = fopen2(output_file);

    fprintf(fp, "%d %d %.6e\n", cycles, nno, E->elapsed_time);
    for(n = 1; n <= nno; n++)
        fprintf(fp, "%.6e\n", E->T[n]);
    fclose(fp);
}


/*
 * Binary Output Routines (HDF5)
 */

void init_hdf5(struct All_variables *E)
{

}

void output_hdf5(struct All_variables *E)
{
}

void output_velocity_hdf5(struct All_variables *E)
{

}

void output_temperature_hdf5(struct All_variables *E)
{

}

void finalize_hdf5(struct All_variables *E)
{

}


/*
 * Some small routines.
 */
void parallel_process_termination(int code)
{
    MPI_Finalize();
    exit(code);
}

FILE *fopen2(char *filename)
{
    FILE *fp;
    if(*filename)
    {
        fp = fopen(filename, "w");
        if(!fp)
        {
            fprintf(stderr, "Cannot open file '%s'\n", filename);
            parallel_process_termination(EXIT_FAILURE);
        }
    }
    else
    {
        fp = stderr;
    }
    return fp;
}


/* non-runaway malloc */
void *safe_malloc(int bytes, char *file, int line)
{
    void *ptr;
    ptr = malloc((size_t)bytes);
    if(ptr == (void *)NULL)
    {
        fprintf(stderr, "safe_malloc: cannot allocate another %d bytes "
                "(line %d of file %s)\n", bytes, line, file);
        exit(0);
    }
#ifdef DEBUG
    fprintf(stderr, 
            "safe_malloc: allocated %d bytes (line %d of file %s)\n",
            bytes, line, file);
#endif
    return ptr;
}
