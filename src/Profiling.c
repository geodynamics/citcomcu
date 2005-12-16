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

/* Profiling functions .... return elapsed CPU time etc. 
   These functions seem the most likely to get broken by
   different architectures/operating systems &c */

#if defined(__osf__) || defined(__aix__) || defined(__sunos__) || defined(__sgi)
#include <sys/time.h>
#include <sys/resource.h>
#define RUSAGE_STYLE_TIME
#elif defined(_UNICOS) || defined(__hpux)
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#define TIMES_STYLE_TIME
#endif

/* ===============================================
   Function to return currently elapsed CPU usage
   =============================================== */

float CPU_time(void)
{
#if defined(RUSAGE_STYLE_TIME)
	struct rusage rusage;
	double time;

	getrusage(RUSAGE_SELF, &rusage);
	time = rusage.ru_utime.tv_sec + 1.0e-6 * rusage.ru_utime.tv_usec;
#elif defined(TIMES_STYLE_TIME)

	struct tms time_now;
	time_t utime;
	long sometime;

	float time;
	static float initial_time;
	static int visit = 0;

	if(visit == 0)
	{
		sometime = times(&time_now);
		initial_time = (float)time_now.tms_utime / (float)CLK_TCK;
		visit++;
	}

	sometime = times(&time_now);
	time = (float)time_now.tms_utime / (float)CLK_TCK - initial_time;
	return (time);

#else /* stupid, break nothing "timer" */
	static float time;
	time += 0.0001;
#endif

	return (float)time;
}
