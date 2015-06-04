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

void printSubspace(CharacterTable *ct, msym_subspace_t *ss);
void freeSubspace(msym_subspace_t *ss);
int filterSubspace(msym_subspace_t *ss);
void printTransform(int r, int c, double M[r][c]);

#endif /* defined(__MSYM__SUBSPACE_h) */

