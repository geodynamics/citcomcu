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

struct CONVECTION
{								/* information controlling convection problems */
	char old_T_file[100];

	float temp_blob_x[40];
	float temp_blob_y[40];
	float temp_blob_z[40];
	float temp_blob_radius[40];	/* +/- */
	float temp_blob_T[40];
	float temp_blob_bg[40];		/* Reference level if sticky */
	int temp_blob_sticky[40];
	int temp_blobs;

	float temp_zonex1[40];
	float temp_zonex2[40];
	float temp_zonez1[40];
	float temp_zonez2[40];
	float temp_zoney1[40];
	float temp_zoney2[40];
	float temp_zonehw[40];
	float temp_zonemag[40];
	int temp_zone_sticky[40];
	int temp_zones;

  int random_t_init;

	float half_space_age;
	int half_space_cooling;

	int number_of_perturbations;
	float perturb_mag[33];
	float perturb_k[33];
	float perturb_ll[33];
	float perturb_mm[33];

	struct SOURCES
	{
		int number;
		float t_offset;
		float Q[10];
		float lambda[10];
	} heat_sources;

	float elasticity1;

} convection;
