#include <assert.h>
#include <mpi.h>
#include <stdbool.h>
#include <string.h>

#include "output_vtk.h"


#define LONG_VTK_LINE 256


static void output_vtk_calculate_extents( struct All_variables * E, int ** extent );
static void output_vtk_print_vts_coords( struct All_variables * E, FILE * data_file, float * w1, float * w2, float * w3 );

static inline void output_vtk_sort_vts_coords_w1( struct All_variables * E, float * restrict x1, float * restrict xcoords );
static inline void output_vtk_sort_vts_coords_w2( struct All_variables * E, float * restrict x2, float * restrict ycoords );
static inline void output_vtk_sort_vts_coords_w3( struct All_variables * E, float * restrict x3, float * restrict zcoords );

static bool output_vtk_write_pvts( struct All_variables * E );
static void output_vtk_write_vts_coords( struct All_variables * E, FILE * data_file, float * w1, float * w2, float * w3 );
static void output_vtk_write_vts_data( struct All_variables * E, FILE * data_file );

static bool output_vtk_write_vts_output( struct All_variables *E );
static void output_vtk_write_vts_epilog( FILE *data_file );
static void output_vtk_write_vts_prolog( struct All_variables *E, FILE *data_file );

__attribute__((unused)) static void output_vtk_write_vts_scalar_d( struct All_variables * E, FILE * data_file, char * name, double * w );
static void output_vtk_write_vts_scalar_f( struct All_variables * E, FILE * data_file, char * name, float * w );
__attribute__((unused)) static void output_vtk_write_vts_scalar_i( struct All_variables * E, FILE * data_file, char * name, int * w );

__attribute__((unused)) static void output_vtk_write_vts_vector_d( struct All_variables * E, FILE * data_file, char * name, double * w1, double * w2, double * w3 );
static void output_vtk_write_vts_vector_f( struct All_variables * E, FILE * data_file, char * name, float * w1, float * w2, float * w3 );
__attribute__((unused)) static void output_vtk_write_vts_vector_i( struct All_variables * E, FILE * data_file, char * name, int * w1, int * w2, int * w3 );


bool output_vtk( struct All_variables * E )
{
	bool written = false;
	written = output_vtk_write_vts_output( E );  assert( written );
	written = output_vtk_write_pvts( E );  assert( written );

	return written;
}


/* static functions below */


static void output_vtk_calculate_extents( struct All_variables * E, int ** extent )
{
	int current_proc;
	int elx, ely, elz;
	int w1, w2, w3;

	elx = E -> mesh.elx / E -> parallel.nprocx;
	ely = E -> mesh.ely / E -> parallel.nprocy;
	elz = E -> mesh.elz / E -> parallel.nprocz;

	current_proc = 0;
	for( w3 = 0; w3 < E -> parallel.nprocz; ++w3 )
	{
		for( w2 = 0; w2 < E -> parallel.nprocy; ++w2 )
		{
			for( w1 = 0; w1 < E -> parallel.nprocx; ++w1 )
			{
				extent[ current_proc ][0] = elx * ( w1 );
				extent[ current_proc ][1] = elx * ( w1 + 1 );
				extent[ current_proc ][2] = ely * ( w2 );
				extent[ current_proc ][3] = ely * ( w2 + 1 );
				extent[ current_proc ][4] = elz * ( w3 );
				extent[ current_proc ][5] = elz * ( w3 + 1 );

				++current_proc;
			}
		}
	}
}


static void output_vtk_print_vts_coords( struct All_variables *E, FILE *data_file, float * w1, float * w2, float * w3 )
{
	int i, j, k;

	fprintf( data_file, "%s\n", "<Points>" );
	fprintf( data_file, "%s\n", "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" );

	for( k = 0; k != E -> lmesh.noz; ++k )
	{
		for( j = 0; j != E -> lmesh.noy; ++j )
		{
			for( i = 0; i != E -> lmesh.nox; ++i )
			{
				fprintf( data_file, "%.5e %.5e %.5e\n", w1[i], w2[j], w3[k] );
			}
		}
	}

	fprintf( data_file, "%s\n", "</DataArray>" );
	fprintf( data_file, "%s\n", "</Points>" );
}


static inline void output_vtk_sort_vts_coords_w1( struct All_variables * E, float * restrict w1, float * restrict xcoords )
{
	int i, j;

	/*  x-coords */
	for( i = 1, j = 0; i <= E -> lmesh.noz * E -> lmesh.nox; ++i )
	{
		if( i % E -> lmesh.noz == 0 )
		{
			xcoords[j] = w1[i];
			++j;
		}
	}
}


static inline void output_vtk_sort_vts_coords_w2( struct All_variables * E, float * restrict w2, float * restrict ycoords )
{
	int i, j;

	/*  y-coords */
	for( i = 1, j = 0; i <= E -> lmesh.nno; ++i )
	{
		if( i % ( E -> lmesh.noz * E -> lmesh.nox ) == 0 )
		{
			ycoords[j] = w2[i];
			++j;
		}
	}
}


static inline void output_vtk_sort_vts_coords_w3( struct All_variables * E, float * restrict w3, float * restrict zcoords )
{
	int i, j;

	/*  z-coords */
	for( i = 1, j = 0; i <= E -> lmesh.noz; ++i, ++j )
	{
		zcoords[j] = w3[i];
	}
}


static bool output_vtk_write_pvts( struct All_variables * E )
{
	int i = 0;
	int ret_i = 0;
	bool written = false;

	char output_file[ 512 ];
	FILE * pvts_file = NULL;
	int ** extent = NULL;

	if( E -> parallel.me == 0 )
	{
		extent = ( int ** ) malloc( ( E -> parallel.nproc ) * sizeof( int * ) );  assert( extent != NULL );
		for( i = 0; i < E -> parallel.nproc; ++i )
		{
			extent[i] = NULL;
			extent[i] = ( int * ) malloc( 6 * sizeof( int ) );  assert( extent[i] != NULL );
		}
		output_vtk_calculate_extents( E, extent );

		sprintf( output_file, "%s/%s%08d.pvts",
					E -> control.output_path,
					E -> control.experiment_name,
					E -> monitor.solution_cycles );

		if( ( pvts_file = fopen( output_file, "w" ) ) != NULL )
		{
			fprintf( pvts_file, "<?xml version=\"1.0\"?>\n" );
			fprintf( pvts_file, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n" );
			fprintf( pvts_file, "<PStructuredGrid GhostLevel=\"0\" WholeExtent=\"0 %d 0 %d 0 %d\">\n",
										E -> mesh.nox - 1,
										E -> mesh.noy - 1,
										E -> mesh.noz - 1 );

			fprintf( pvts_file, "  <PCellData>\n" );
			fprintf( pvts_file, "  </PCellData>\n" );

			fprintf( pvts_file, "  <PPoints>\n" );
			fprintf( pvts_file, "     <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n" );
			fprintf( pvts_file, "  </PPoints>\n" );
			
			fprintf( pvts_file, "  <PPointData Scalars=\"Temperature\" Vectors=\"Velocity\">\n" );
			fprintf( pvts_file, "     <PDataArray type=\"Float32\" Name=\"Temperature\"   format=\"ascii\" NumberOfComponents=\"1\"/>\n" );
			fprintf( pvts_file, "     <PDataArray type=\"Float32\" Name=\"Heatflux\"      format=\"ascii\" NumberOfComponents=\"1\"/>\n" );
			fprintf( pvts_file, "     <PDataArray type=\"Float32\" Name=\"Heatflux_adv\"  format=\"ascii\" NumberOfComponents=\"1\"/>\n" );
			fprintf( pvts_file, "     <PDataArray type=\"Float32\" Name=\"Velocity\"      format=\"ascii\" NumberOfComponents=\"3\"/>\n" );
			fprintf( pvts_file, "  </PPointData>\n" );

			for( i = 0; i < E -> parallel.nproc; ++i )
			{
				fprintf( pvts_file, "    <Piece Extent=\" %d %d %d %d %d %d\" Source=\"%s%05d_%08d.vts\"/>\n",
										extent[i][0], extent[i][1],
					 					extent[i][2], extent[i][3],
					 					extent[i][4], extent[i][5],
										E -> control.experiment_name,
										i,
										E -> monitor.solution_cycles );
			}
			fprintf( pvts_file, "  </PStructuredGrid>\n" );
			fprintf( pvts_file, "</VTKFile>\n" );

			fclose( pvts_file );

			for( i = 0; i < E -> parallel.nproc; ++i )
			{
				free( extent[i] );
				extent[i] = NULL;
			}
			free( extent );
			
			written = true;
		}

		ret_i = written ? 1 : 0;
	}

	MPI_Bcast( &ret_i, 1, MPI_INT, 0, MPI_COMM_WORLD );
	written = ( ret_i == 1 ? true : false );

	return written;
}


static void output_vtk_write_vts_coords( struct All_variables * E, FILE * data_file, float * w1, float * w2, float * w3 )
{
	/* v1, v2, v3 will contain the sorted (for vtk/vts) cartesian or spherical coordinates */
	float * v1 = NULL;
   	float * v2 = NULL;
	float * v3 = NULL;

	v1 = malloc( sizeof( float ) * ( E -> lmesh.nox ) );  assert( v1 != NULL );
	v2 = malloc( sizeof( float ) * ( E -> lmesh.noy ) );  assert( v2 != NULL );
	v3 = malloc( sizeof( float ) * ( E -> lmesh.noz ) );  assert( v3 != NULL );

	output_vtk_sort_vts_coords_w1( E, w1, v1 );
	output_vtk_sort_vts_coords_w2( E, w2, v2 );
	output_vtk_sort_vts_coords_w3( E, w3, v3 );

	output_vtk_print_vts_coords( E, data_file, v1, v2, v3 );

	free( v3 );
	free( v2 );
	free( v1 );
}


static void output_vtk_write_vts_data( struct All_variables * E, FILE * data_file )
{
	fprintf( data_file, "%s\n", "<PointData Scalars = \"Temperature\" Vectors = \"Velocity\">" );
	output_vtk_write_vts_scalar_f( E, data_file, "Temperature", E -> T );
	output_vtk_write_vts_scalar_f( E, data_file, "Heatflux", E -> heatflux );
	output_vtk_write_vts_scalar_f( E, data_file, "Heatflux_adv", E -> heatflux_adv );
	output_vtk_write_vts_vector_f( E, data_file, "Velocity", E -> V[1], E -> V[2], E -> V[3] );
	fprintf( data_file, "%s\n", "</PointData>" );
}


static bool output_vtk_write_vts_output( struct All_variables *E )
{
	bool written = true;
	char output_file[ 512 ];
	FILE *data_file = NULL;

	sprintf( output_file, "%s/%s%05d_%08d.vts",
			E -> control.output_path,
			E -> control.experiment_name,
			E -> parallel.me,
			E -> monitor.solution_cycles );

	if( ( data_file = fopen( output_file, "w" ) ) != NULL )
	{
		output_vtk_write_vts_prolog( E, data_file );
		if( E -> control.CART3D )
		{
			output_vtk_write_vts_coords( E, data_file, E -> X[1], E -> X[2], E -> X[3] );
		}
		else
		{
			output_vtk_write_vts_coords( E, data_file, E -> SX[1], E -> SX[2], E -> SX[3] );
		}

		output_vtk_write_vts_data( E, data_file );

		output_vtk_write_vts_epilog( data_file );

		fclose( data_file );
	}
	else
	{
		written = false;
	}

	return written;
}


static void output_vtk_write_vts_epilog( FILE *data_file )
{
	fprintf( data_file, "%s\n", "</Piece>" );
	fprintf( data_file, "%s\n", "</StructuredGrid> " );
	fprintf( data_file, "%s\n", "</VTKFile>" );
}


static void output_vtk_write_vts_prolog( struct All_variables *E, FILE *data_file )
{
	fprintf( data_file, "%s\n", "<?xml version=\"1.0\"?>" );
	fprintf( data_file, "%s\n", "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" );

	fprintf( data_file, "%s%d %d %d %d %d %d%s\n",
				"<StructuredGrid WholeExtent=\"",
				( E -> mesh.nox - 1 ) * 0,
				( E -> mesh.nox - 1 ) * 1,
				( E -> mesh.noy - 1 ) * 0,
				( E -> mesh.noy - 1 ) * 1,
				( E -> mesh.noz - 1 ) * 0,
				( E -> mesh.noz - 1 ) * 1,
				"\">" );

	fprintf( data_file, "%s%d %d %d %d %d %d%s\n",
				"<Piece Extent=\"",
				( E -> lmesh.nox - 1 ) * ( E -> parallel.me_loc[1] ),
				( E -> lmesh.nox - 1 ) * ( E -> parallel.me_loc[1] + 1 ),
				( E -> lmesh.noy - 1 ) * ( E -> parallel.me_loc[2] ),
				( E -> lmesh.noy - 1 ) * ( E -> parallel.me_loc[2] + 1 ),
				( E -> lmesh.noz - 1 ) * ( E -> parallel.me_loc[3] ),
				( E -> lmesh.noz - 1 ) * ( E -> parallel.me_loc[3] + 1 ),
				"\">" );

	fprintf( data_file, "%s\n", "<CellData></CellData>" );
}



__attribute__((unused)) static void output_vtk_write_vts_scalar_d( struct All_variables * E, FILE * data_file, char * name, double * w )
{
	int i, j, k, l;
	char line[ LONG_VTK_LINE ];

	sprintf( line, "%s", "<DataArray type=\"Float32\" Name=\"" );
	strcat( line, name );
	strcat( line, "\" format=\"ascii\" NumberOfComponents=\"1\">" );
	fprintf( data_file, "%s\n", line );

	/* data is stored in zxy direction, paraview wants it in xyz */
	for( k = 0; k < E -> lmesh.noz; ++k )
	{
		for( j = 0; j < E -> lmesh.noy; ++j )
		{
			for( i = 0; i < E -> lmesh.nox; ++i )
			{
				l = k + j * E -> lmesh.nox * E -> lmesh.noz + i * E -> lmesh.noz;
				l = ( l % ( E -> lmesh.nno ) ) + 1;
				fprintf( data_file, "%lf\n", w[l] );
			}
		}
	}

	fprintf( data_file, "%s\n", "</DataArray>" );
}


static void output_vtk_write_vts_scalar_f( struct All_variables * E, FILE * data_file, char * name, float * w )
{
	int i, j, k, l;
	char line[ LONG_VTK_LINE ];

	sprintf( line, "%s", "<DataArray type=\"Float32\" Name=\"" );
	strcat( line, name );
	strcat( line, "\" format=\"ascii\" NumberOfComponents=\"1\">" );
	fprintf( data_file, "%s\n", line );

	/* data is stored in zxy direction, paraview wants it in xyz */
	for( k = 0; k < E -> lmesh.noz; ++k )
	{
		for( j = 0; j < E -> lmesh.noy; ++j )
		{
			for( i = 0; i < E -> lmesh.nox; ++i )
			{
				l = k + j * E -> lmesh.nox * E -> lmesh.noz + i * E -> lmesh.noz;
				l = ( l % ( E -> lmesh.nno ) ) + 1;
				fprintf( data_file, "%f\n", w[l] );
			}
		}
	}

	fprintf( data_file, "%s\n", "</DataArray>" );
}


__attribute__((unused)) static void output_vtk_write_vts_scalar_i( struct All_variables * E, FILE * data_file, char * name, int * w )
{
	int i, j, k, l;
	char line[ LONG_VTK_LINE ];

	sprintf( line, "%s", "<DataArray type=\"Float32\" Name=\"" );
	strcat( line, name );
	strcat( line, "\" format=\"ascii\" NumberOfComponents=\"1\">" );
	fprintf( data_file, "%s\n", line );

	/* data is stored in zxy direction, paraview wants it in xyz */
	for( k = 0; k < E -> lmesh.noz; ++k )
	{
		for( j = 0; j < E -> lmesh.noy; ++j )
		{
			for( i = 0; i < E -> lmesh.nox; ++i )
			{
				l = k + j * E -> lmesh.nox * E -> lmesh.noz + i * E -> lmesh.noz;
				l = ( l % ( E -> lmesh.nno ) ) + 1;
				fprintf( data_file, "%d\n", w[l] );
			}
		}
	}

	fprintf( data_file, "%s\n", "</DataArray>" );
}


__attribute__((unused)) static void output_vtk_write_vts_vector_d( struct All_variables * E, FILE * data_file, char * name, double * w1, double * w2, double * w3 )
{
	int i, j, k, l;
	char line[ LONG_VTK_LINE ];

	sprintf( line, "%s", "<DataArray type=\"Float32\" Name=\"" );
	strcat( line, name );
	strcat( line, "\" format=\"ascii\" NumberOfComponents=\"3\">" );
	fprintf( data_file, "%s\n", line );

	/* data is stored in zxy direction, paraview wants it in xyz */
	for( k = 0; k < E -> lmesh.noz; ++k )
	{
		for( j = 0; j < E -> lmesh.noy; ++j )
		{
			for( i = 0; i < E -> lmesh.nox; ++i )
			{
				l = k + j * E -> lmesh.nox * E -> lmesh.noz + i * E -> lmesh.noz;
				l = ( l % ( E -> lmesh.nno ) ) + 1;
				fprintf( data_file, "%.6e %.6e %.6e\n", w1[l], w2[l], w3[l] );
			}
		}
	}

	fprintf( data_file, "%s\n", "</DataArray>" );
}


__attribute__((unused)) static void output_vtk_write_vts_vector_f( struct All_variables * E, FILE * data_file, char * name, float * w1, float * w2, float * w3 )
{
	int i, j, k, l;
	char line[ LONG_VTK_LINE ];

	sprintf( line, "%s", "<DataArray type=\"Float32\" Name=\"" );
	strcat( line, name );
	strcat( line, "\" format=\"ascii\" NumberOfComponents=\"3\">" );
	fprintf( data_file, "%s\n", line );

	/*  data is stored in zxy direction, paraview wants it in xyz */
	for( k = 0; k < E -> lmesh.noz; ++k )
	{
		for( j = 0; j < E -> lmesh.noy; ++j )
		{
			for( i = 0; i < E -> lmesh.nox; ++i )
			{
				l = k + j * E -> lmesh.nox * E -> lmesh.noz + i * E -> lmesh.noz;
				l = ( l % ( E -> lmesh.nno ) ) + 1;
				fprintf( data_file, "%.6e %.6e %.6e\n", w1[l], w2[l], w3[l] );
			}
		}
	}

	fprintf( data_file, "%s\n", "</DataArray>" );
}


__attribute__((unused)) static void output_vtk_write_vts_vector_i( struct All_variables * E, FILE * data_file, char * name, int * w1, int * w2, int * w3 )
{
	int i, j, k, l;
	char line[ LONG_VTK_LINE ];

	sprintf( line, "%s", "<DataArray type=\"Float32\" Name=\"" );
	strcat( line, name );
	strcat( line, "\" format=\"ascii\" NumberOfComponents=\"3\">" );
	fprintf( data_file, "%s\n", line );

	/*  data is stored in zxy direction, paraview wants it in xyz */
	for( k = 0; k < E -> lmesh.noz; ++k )
	{
		for( j = 0; j < E -> lmesh.noy; ++j )
		{
			for( i = 0; i < E -> lmesh.nox; ++i )
			{
				l = k + j * E -> lmesh.nox * E -> lmesh.noz + i * E -> lmesh.noz;
				l = ( l % ( E -> lmesh.nno ) ) + 1;
				fprintf( data_file, "%d %d %d\n", w1[l], w2[l], w3[l] );
			}
		}
	}

	fprintf( data_file, "%s\n", "</DataArray>" );
}
