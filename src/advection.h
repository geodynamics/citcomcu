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
struct ADVECTION
{
	int ADVECTION;

	float gamma;
	float timestep;
	float fine_tune_dt;
	float dt_reduced;
	float fixed_timestep;
	float max_dimensionless_time;

	int markers_g, markers, markers1, markers_uplimit, markers_per_ele;
	float marker_maxdist, marker_mindist;
	int markerIX, markerIY, markerIZ;


	int min_timesteps;
	int max_timesteps;
	int max_total_timesteps;
	int timesteps;
	int total_timesteps;
	int temp_iterations;
	int max_substeps;
	int sub_iterations;
	int last_sub_iterations;

	float vel_substep_aggression;
	float temp_updatedness;
	float visc_updatedness;

	float lid_defining_velocity;
	float sub_layer_sample_level;

} advection;
