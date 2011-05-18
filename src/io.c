#include <assert.h>
#include <errno.h>
#include <math.h>
#include <mpi.h>

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>

#include <unistd.h>

#include "io.h"
#include "output_ascii.h"
#include "output_vtk.h"


#define MAXPATHLENGTH 256
#define MAXFILENAMELENGTH 256 /* TODO: sync it with global_defs.h */


static bool io_setup_filename_root( struct All_variables * );

/**
 * create a directory
 * @param const char *, directory name
 * @return true or false
 */
bool io_directory_create( const char * dir_name )
{
	int i = 0;
	int length = 0;
	int status = 0;

	char * tmp = NULL;

	assert( dir_name );
	length = strlen( dir_name );
	tmp = (char *) malloc( (length+1) * sizeof( char ) );
	assert( tmp );
	memset( tmp, '\0', (length+1) );


	for( i = 0; i != length; ++i )
	{
		tmp[i] = dir_name[i];

		if( dir_name[i] == '/' )
			if( !io_directory_exists( tmp ) )
				status = mkdir( tmp, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

		if( status != 0 )
		{
			io_error( );
			return false;
		}
	}

	free( tmp );
	/* 
	 * if dirName does not contain trailing / above
	 * loop does not create final directory
	 */
	if( !io_directory_exists( dir_name ) )
	{
		status = mkdir( dir_name, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
	}

	if( status == 0 )
	{
		return true;
	}
	
	io_error( );

	return false;
}


/**
 * check if a directory exists
 * @param const char *, directory name
 * @return true or false
 */
bool io_directory_exists( const char * dir_name )
{
	struct stat dir_info;

	assert( dir_name );
	if( stat( dir_name, &dir_info ) < 0 )
	{
		return false;
	}
	else
	{
		return ( dir_info.st_mode & S_IFDIR ); /* S_IFDIR requires -D_GNU_SOURCE */
	}
}

/**
 *
 */
void io_error( )
{
}

/**
 * inferface for the output of results
 * @param struct All_variables * E
 * @return success or failure
 */ 
bool io_results( struct All_variables * E )
{
	static bool first_time = true;
	bool result = false;

	if( first_time )
	{
		int proceed = 1;
		if( E -> parallel.me == 0 )
		{
			if( !io_setup_path( E ) )
			{
				fprintf( stdout, "failed to setup output directory %s\n", E -> control.output_path );
				fflush( stdout );
				
				proceed = 0;
			}
		}
		MPI_Bcast( &proceed, 1, MPI_INT, 0, MPI_COMM_WORLD );

		if( !proceed )
		{
			MPI_Finalize( );
			exit( EXIT_FAILURE );
		}

		io_setup_filename_root( E );
		first_time = false;
	}

	if( E -> control.output_ascii )
	{
		result = output_ascii( E );
		assert( result == true );
	}

	if( E -> control.output_vtk )
	{
		result = output_vtk( E );
		assert( result == true );
	}

	return result;
}

/**
 * setup the output path for results
 * @param struct All_variables * E
 * @retrun success or failure
 */
bool io_setup_path( struct All_variables * E )
{
	char output_path[ MAXPATHLENGTH ];
	int length = 0;

	if( E -> control.output_path[0] != '/' )
	{
		sprintf( output_path, "%s", "./" );
	}

	length = strlen( E -> control.output_path );
	if( length > MAXPATHLENGTH - 3 )
	{
		fprintf( stderr, "path length (%d), too long\n", length );
		exit( EXIT_FAILURE );
	}

	strncat( output_path, E -> control.output_path, length );

	if( !io_directory_exists( output_path ) )
	{
		if( !io_directory_create( output_path ) )
		{
			return false;
		}
	}

	sprintf( E -> control.output_path, "%s", output_path );

	return true;
}

static bool io_setup_filename_root( struct All_variables * E )
{
	char time_string[13];

	struct timeval tv;
	time_t curtime;

	/* 15 = "_" + "201011312359\0" + "_" */
	assert( MAXFILENAMELENGTH > 15 );  /* sanity check, in case it is unintentionally modified */
	assert( strlen( E -> control.experiment_name ) < MAXFILENAMELENGTH - 15 );

	strcat( E -> control.experiment_name, "_" );

	gettimeofday( &tv, NULL );
	curtime = tv.tv_sec;
	strftime( time_string, 13, "%Y%m%d%H%M", localtime( &curtime ) );
	strcat( E -> control.experiment_name, time_string );

	strcat( E -> control.experiment_name, "_" );

	return true;
}
