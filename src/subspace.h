//
//  subspace.h
//  libmsym
//
//  Created by Marcus Johansson on 28/05/15.
//  Copyright (c) 2014 Marcus Johansson.
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#ifndef __MSYM__SUBSPACE_h
#define __MSYM__SUBSPACE_h

#include <stdio.h>
#include "msym.h"
#include "character_table.h"

void printTransform(int r, int c, double M[r][c]);
void freeSubrepresentationSpaces(int srsl, msym_subrepresentation_space_t *srs);

#endif /* defined(__MSYM__SUBSPACE_h) */

