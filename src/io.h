#ifndef IO_HEADER_H
#define IO_HEADER_H

#include <stdbool.h>

#include "global_defs.h"

/**
 * create a directory
 * @param const char *, directory name
 * @return true or false
 */
bool io_directory_create( const char * );

/**
 * check if a directory exists
 * @param const char *, directory name
 * @return true or false
 */
bool io_directory_exists( const char * );

/**
 *
 */
void io_error( );

/**
 * inferface for the output of results
 * @param struct All_variables * E
 * @return success or failure
 */ 
bool io_results( struct All_variables * );

/**
 * setup the output path for results
 * @param struct All_variables * E
 * @retrun success or failure
 */
bool io_setup_path( struct All_variables * );


#endif
