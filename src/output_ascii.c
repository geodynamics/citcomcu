#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "output_ascii.h"

#define MAX_FILENAME_LENGTH 256


static bool output_ascii_ave( struct All_variables * E );
static bool output_ascii_coords( struct All_variables * E );
static bool output_ascii_temp( struct All_variables * E );
static bool output_ascii_th_b( struct All_variables * E );

static bool output_ascii_th_t( struct All_variables * E );
static bool output_ascii_traces( struct All_variables * E );
static bool output_ascii_velocity( struct All_variables * E );
static bool output_print_coords( FILE * fd, int number_of_coords, float * x0, float * x1, float * x2 );


/**
 * interface for the output of results in ascii format
 * @param E pointer to struct All_variables
 * @return success or failure
 */
bool output_ascii( struct All_variables * E )
{
	bool written = true;
	static bool first_time = true;

	if( first_time )
	{
		first_time = false;
		written = output_ascii_coords( E );
		assert( written );
	}
	
	if( E -> parallel.me < E -> parallel.nprocz )
	{
		written = output_ascii_ave( E );
		assert( written );
	}

	if( ( E -> monitor.solution_cycles % E -> control.record_every ) == 0 )
	{
		written = output_ascii_temp( E );
		assert( written );

		written = output_ascii_velocity( E );
		assert( written );


		if( E -> parallel.me_loc[3] == E -> parallel.nprocz - 1 )
		{
			written = output_ascii_th_t( E );
			assert( written );
		}

		if( E -> parallel.me_loc[3] == 0 )
		{
			written = output_ascii_th_b( E );
			assert( written );
		}
	}

	if( E -> control.composition && E -> monitor.solution_cycles % ( 10 * E -> control.record_every ) == 0 )
	{
		written = output_ascii_traces( E );
		assert( E );
	}

	return written;
}


static bool output_ascii_ave( struct All_variables * E )
{
	/* TODO: ave == averages? of what? */
	bool written = true;

	int i;
	char output_filename[ MAX_FILENAME_LENGTH ];
	FILE *fd = NULL;

	memset( output_filename, '\0', MAX_FILENAME_LENGTH );
	sprintf( output_filename, "%s/%s%05d_%08d.ave",
				E -> control.output_path,
				E -> control.experiment_name,
				E -> parallel.me,
				E -> monitor.solution_cycles );
	assert( strlen( output_filename ) < MAX_FILENAME_LENGTH );

	if( ( fd = fopen( output_filename, "w" ) ) != NULL )
	{
		fprintf( fd, "%6d %6d %.5e %.5e %.5e %.4e %.4e %.5e %.5e\n",
				E -> lmesh.noz,
				E -> advection.timesteps, E -> monitor.elapsed_time,
				E -> slice.Nut, E -> slice.Nub,
				E -> data.T_adi0, E -> data.T_adi1,
				E -> monitor.Sigma_interior, E -> monitor.Sigma_max );

		for( i = 1; i <= E -> lmesh.noz; i++ )
		{
			fprintf( fd, "%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n",
					E -> Have.T[i],
					E -> Have.vrms[i],
					E -> Have.Vi[i],
					E -> Have.Rho[i],
					E -> Have.F[i],
					E -> Have.f[i],
					E -> Have.C[i],
					E -> Have.Tadi[i] );
		}
	}
	else
	{
		written = false;
		fprintf( stderr, "failed to open file %s\n", output_filename );
		fflush( stderr );
	}

	fclose( fd );
	return written;
}

static bool output_ascii_coords( struct All_variables * E )
{
	bool written = true;

	char output_filename[ MAX_FILENAME_LENGTH ];
	FILE *fd = NULL;

	memset( output_filename, 0, MAX_FILENAME_LENGTH );
	sprintf( output_filename, "%s/%s%05d.coordinates",
				E -> control.output_path,
				E -> control.experiment_name,
				E -> parallel.me );
	assert( strlen( output_filename ) < MAX_FILENAME_LENGTH );

	if( ( fd = fopen( output_filename, "w" ) ) != NULL )
	{
		fprintf( fd, "%6d\n", E -> lmesh.nno );
		if( E -> control.CART3D )
		{
			written = output_print_coords( fd, E -> lmesh.nno, E -> X[1], E -> X[2], E -> X[3] );
		}
		else
		{
			written = output_print_coords( fd, E -> lmesh.nno, E -> SX[1], E -> SX[2], E -> SX[3] );
		}
	}
	else
	{
		written = false;
		fprintf( stderr, "failed to open file %s\n", output_filename );
		fflush( stderr );
	}

	fclose( fd );
	return written;
}

static bool output_print_coords( FILE * fd, int number_of_coords, float * x_0, float * x_1, float * x_2 )
{
	int i = 1;
	int result = 0;

	for( i = 1; i <= number_of_coords; ++i )
	{
		result = fprintf( fd, "%.5e %.5e %.5e\n", x_0[i], x_1[i], x_2[i] );
	}

	return ( result > 0 ? true : false );
}

static bool output_ascii_temp( struct All_variables * E )
{
	bool written = true;

	int i;
	char output_filename[ MAX_FILENAME_LENGTH ];
	FILE *fd = NULL;

	memset( output_filename, 0, MAX_FILENAME_LENGTH );
	sprintf( output_filename, "%s/%s%05d_%08d.temperature",
				E -> control.output_path,
				E -> control.experiment_name,
				E -> parallel.me,
				E -> monitor.solution_cycles );
	assert( strlen( output_filename ) < MAX_FILENAME_LENGTH );

	if( ( fd = fopen( output_filename, "w" ) ) != NULL )
	{
		fprintf( fd, "%6d %6d %.5e\n", E -> lmesh.nno, E -> advection.timesteps, E -> monitor.elapsed_time );

		if( E -> monitor.solution_cycles % ( 20 * E -> control.record_every ) == 0 )
		{
			if( E -> control.composition == 0 )
			{
				for( i = 1; i <= E -> lmesh.nno; i++ )
				{
					/* note in original, the second arg to this fprintf was E -> V[3][i], not sure why */
					fprintf( fd, "%.5e %.4e %.4e\n", E -> T[i], E -> heatflux[i], E -> heatflux_adv[i] );
				}
			}
			else if( E -> control.composition )
			{
				for( i = 1; i <= E -> lmesh.nno; i++ )
				{
					/* note in original, the second arg to this fprintf was E -> V[3][i], not sure why */
					fprintf( fd, "%.5e %.4e %.4e %.4e\n", E -> T[i], E -> heatflux[i], E -> heatflux_adv[i], E -> C[i] );
				}
			}
		}
		else
		{
			if( E -> control.composition == 0 )
			{
				for( i = 1; i <= E -> lmesh.nno; i++ )
				{
					fprintf( fd, "%.5e %.4e %.4e\n", E -> T[i], E -> V[3][i], E -> heatflux_adv[i] );
				}
			}
			else if( E -> control.composition )
			{
				for( i = 1; i <= E -> lmesh.nno; i++ )
				{
					fprintf( fd, "%.5e %.4e %.4e %.4e\n", E -> T[i], E -> V[3][i], E -> heatflux_adv[i], E -> C[i] );
				}
			}
		}
	}
	else
	{
		written = false;
		fprintf( stderr, "failed to open file %s\n", output_filename );
		fflush( stderr );
	}

	fclose( fd );
	return written;
}

static bool output_ascii_th_b( struct All_variables * E )
{
	/* TODO: find out what th_b is */
	bool written = true;

	int i;
	char output_filename[ MAX_FILENAME_LENGTH ];
	FILE *fd = NULL;

	memset( output_filename, 0, MAX_FILENAME_LENGTH );
	sprintf( output_filename, "%s/%s%05d_%08d.th_b",
				E -> control.output_path,
				E -> control.experiment_name,
				E -> parallel.me,
				E -> monitor.solution_cycles );
	assert( strlen( output_filename ) < MAX_FILENAME_LENGTH );

	if( ( fd = fopen( output_filename, "w" ) ) != NULL )
	{
		fprintf( fd, "%6d %6d %.5e %.5e\n", E -> lmesh.nsf, E -> advection.timesteps, E -> monitor.elapsed_time, E -> slice.Nub );

		for( i = 1; i <= E -> lmesh.nsf; i++ )
		{
			fprintf( fd, "%.5e %.5e %.5e %.5e\n", E -> slice.tpgb[i], E -> slice.bhflux[i], E -> Fas410_b[i], E -> Fas670_b[i] );
		}
	}
	else
	{
		written = false;
		fprintf( stderr, "failed to open file %s\n", output_filename );
		fflush( stderr );
	}

	fclose( fd );
	return written;
}

static bool output_ascii_th_t( struct All_variables * E )
{
	/* TODO: find out what th_t is */
	bool written = true;

	int i;
	char output_filename[ MAX_FILENAME_LENGTH ];
	FILE *fd = NULL;

	memset( output_filename, 0, MAX_FILENAME_LENGTH );
	sprintf( output_filename, "%s/%s%05d_%08d.th_t",
				E -> control.output_path,
				E -> control.experiment_name,
				E -> parallel.me,
				E -> monitor.solution_cycles );
	assert( strlen( output_filename ) < MAX_FILENAME_LENGTH );

	if( ( fd = fopen( output_filename, "w" ) ) != NULL )
	{
		fprintf( fd, "%6d %6d %.5e %.5e\n", E -> lmesh.nsf, E -> advection.timesteps, E -> monitor.elapsed_time, E -> slice.Nut );

		for( i = 1; i <= E -> lmesh.nsf; i++ )
		{
			fprintf( fd, "%.5e %.5e %.5e %.5e\n", E -> slice.tpg[i], E -> slice.shflux[i], E -> Fas410_b[i], E -> Fas670_b[i] );
		}
	}
	else
	{
		written = false;
		fprintf( stderr, "failed to open file %s\n", output_filename );
		fflush( stderr );
	}

	fclose( fd );
	return written;
}


static bool output_ascii_traces( struct All_variables * E )
{
	bool written = true;

	int i;
	char output_filename[ MAX_FILENAME_LENGTH ];
	FILE *fd = NULL;

	memset( output_filename, 0, MAX_FILENAME_LENGTH );
	sprintf( output_filename, "%s/%s%05d_%08d.traces",
				E -> control.output_path,
				E -> control.experiment_name,
				E -> parallel.me,
				E -> monitor.solution_cycles );
	assert( strlen( output_filename ) < MAX_FILENAME_LENGTH );

	if( ( fd = fopen( output_filename, "w" ) ) != NULL )
	{
		fprintf( fd, "%6d %6d %.5e\n", E -> advection.markers, E -> advection.timesteps, E -> monitor.elapsed_time );

		for( i = 1; i <= E -> advection.markers; i++ )
		{
			fprintf( fd, "%g %g %g %d %d\n", E -> XMC[1][i], E -> XMC[2][i], E -> XMC[3][i], E -> CElement[i], E -> C12[i] );
		}

		for( i = 1; i <= E -> lmesh.nel; i++ )
		{
			fprintf( fd, "%g\n", E -> CE[i] );
		}

	}
	else
	{
		written = false;
		fprintf( stderr, "failed to open file %s\n", output_filename );
		fflush( stderr );
	}

	fclose( fd );
	return written;
}

static bool output_ascii_velocity( struct All_variables * E )
{
	bool written = true;

	int i;
	char output_filename[ MAX_FILENAME_LENGTH ];
	FILE *fd = NULL;

	memset( output_filename, 0, MAX_FILENAME_LENGTH );
	sprintf( output_filename, "%s/%s%05d_%08d.velocity",
				E -> control.output_path,
				E -> control.experiment_name,
				E -> parallel.me,
				E -> monitor.solution_cycles );
	assert( strlen( output_filename ) < MAX_FILENAME_LENGTH );

	if( ( fd = fopen( output_filename, "w" ) ) != NULL )
	{

		fprintf( fd, "%6d %6d %.5e\n", E -> lmesh.nno, E -> advection.timesteps, E -> monitor.elapsed_time );

		for( i = 1; i <= E -> lmesh.nno; i++ )
		{
			fprintf( fd, "%.6e %.6e %.6e\n", E -> V[1][i], E -> V[2][i], E -> V[3][i] );
		}
	}
	else
	{
		written = false;
		fprintf( stderr, "failed to open file %s\n", output_filename );
		fflush( stderr );
	}
	fclose( fd );

	return written;
}
