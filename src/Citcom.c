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

/*CITCOM: A finite element convection program written at Caltech 1992 */
/*Aims to include an iterative matrix solver based on Multigrid techniques */
/*To do this requires the use of a mixed method and a conjugate-gradient */
/*approach to determining the */

#include <mpi.h>

#include <math.h>
#include <malloc.h>
#include <unistd.h>
#include <sys/types.h>

#include "element_definitions.h"
#include "global_defs.h"

extern int Emergency_stop;

int main(int argc, char **argv)
{
    struct All_variables E;

    double initial_time;    /* start of variable initializations */
    double start_time;      /* start of calculations */
    double time;

    parallel_process_initialization(&E, argc, argv);
    gethostname(E.parallel.machinename, 160);

    E.monitor.solution_cycles = 0;

    if(E.parallel.me == 0)
    {
        start_time = time = CPU_time0();
    }

    read_instructions(&E, argc, argv);

    E.control.keep_going = 1;

    if(E.parallel.me == 0)
    {
        fprintf(E.fp, "Input parameters taken from file '%s'\n", argv[1]);
        fprintf(E.fp, "Initialization complete after %g seconds\n\n", 
                CPU_time0() - time);
        fflush(E.fp);
        initial_time = CPU_time0() - time;
        fprintf(E.fp, "Initialization overhead = %f\n", initial_time);
        initial_time = CPU_time0();
    }

    general_stokes_solver(&E);
    process_temp_field(&E, E.monitor.solution_cycles);
    process_new_velocity(&E, E.monitor.solution_cycles);

    if(E.control.stokes)
    {
        E.control.keep_going = 0;
        E.monitor.solution_cycles++;
    }


    while(E.control.keep_going && (Emergency_stop == 0))
    {
        process_heating(&E);

        E.monitor.solution_cycles++;
        if(E.monitor.solution_cycles > E.control.print_convergence)
            E.control.print_convergence = 1;

        report(&E, "Update buoyancy for further `timesteps'");
        (E.next_buoyancy_field) (&E);

        report(&E, "Process results of buoyancy update");
        process_temp_field(&E, E.monitor.solution_cycles);

        general_stokes_solver(&E);

        if(E.control.composition)
            (E.next_buoyancy_field)(&E);   /* correct with R-G */

        report(&E, "Process results of velocity solver");
        process_new_velocity(&E, E.monitor.solution_cycles);


        if(E.monitor.T_interior > 1.5)
        {
            fprintf(E.fp, "quit due to maxT = %.4e sub_iteration%d\n",
                    E.monitor.T_interior, E.advection.last_sub_iterations);
            parallel_process_termination(8);
        }

        if(E.parallel.me == 0)
        {
            fprintf(E.fp, 
                    "CPU total = %g & CPU = %g for step %d "
                    "time = %.4e dt = %.4e  maxT = %.4e "
                    "sub_iteration%d markers=%d\n", 
                    CPU_time0() - start_time, CPU_time0() - time,
                    E.monitor.solution_cycles, E.monitor.elapsed_time,
                    E.advection.timestep, E.monitor.T_interior,
                    E.advection.last_sub_iterations, E.advection.markers_g);
            time = CPU_time0();
        }

    }

    if(E.parallel.me == 0)
    {
        time = CPU_time0() - initial_time;
        fprintf(E.fp, "Average cpu time taken for velocity step = %f\n",
                time / ((float)(E.monitor.solution_cycles - 1)));
        fprintf(stderr, "Average cpu time taken for velocity step = %f\n",
                time / ((float)(E.monitor.solution_cycles - 1)));
    }

    fclose(E.fp);

    parallel_process_termination(0);

    return 0;
}
