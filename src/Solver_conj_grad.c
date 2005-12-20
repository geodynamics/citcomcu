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

void set_cg_defaults(struct All_variables *E)
{
	E->build_forcing_term = assemble_forces_iterative;
	E->solve_stokes_problem = solve_constrained_flow_iterative;
	E->solver_allocate_vars = cg_allocate_vars;
	return;
}

void cg_allocate_vars(struct All_variables *E)
{
	/* Nothing required ONLY by conj-grad stuff  */
	/* printf("here here\n"); */
	return;
}

void assemble_forces_iterative(struct All_variables *E)
{
	assemble_forces(E, 0);
	strip_bcs_from_residual(E, E->F, E->mesh.levmax);
	return;
}
