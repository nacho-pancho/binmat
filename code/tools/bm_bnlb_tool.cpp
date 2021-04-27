/*
  * Copyright (c) 2019 Ignacio Francisco Ramírez Paulino and Ignacio Hounie
  * This program is free software: you can redistribute it and/or modify it
  * under the terms of the GNU Affero General Public License as published by
  * the Free Software Foundation, either version 3 of the License, or (at
  * your option) any later version.
  * This program is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero
  * General Public License for more details.
  * You should have received a copy of the GNU Affero General Public License
  *  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/
/**
 * \file bnlb_tool.c
 * \brief Command line interface
 *
 *  Uses GNU argp to parse arguments
 *  see http://www.gnu.org/software/libc/manual/html_node/Argp.html
 */
#include <stdlib.h>
#include <time.h>                     // basic benchmarking
#include <argp.h>                     // argument parsing
#include "bm_binmat.h"
#include "bm_binimage.h"

//----------------------------------------------------------------------------

typedef struct {
    unsigned radius; // patch radius
    float perr; // error probability [0,1)
    const char* input_file;
    const char* output_file;
} bnlb_config;

static void _parse_args ( int argc, char** argv, bnlb_config* cfg);

//----------------------------------------------------------------------------

/**
 * main function
 */
int main ( int argc, char **argv ) {

    bnlb_config cfg;
    //paco_info ( "Parsing arguments...\n" );
    _parse_args ( argc, argv, &cfg );
}

//----------------------------------------------------------------------------

static void _parse_args ( int argc, char** argv, bnlb_config* cfg) {

}
