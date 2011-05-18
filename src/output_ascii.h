#ifndef OUTPUT_ASCII_H
#define OUTPUT_ASCII_H

#include <stdbool.h>

#include "global_defs.h"

/**
 * interface for the output of results in ascii format
 * @param E pointer to struct All_variables
 * @return success or failure
 */
bool output_ascii( struct All_variables * E );

#endif
