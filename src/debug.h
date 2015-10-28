//
//  debug.h
//  libmsym
//
//  Created by Marcus Johansson on 07/11/14.
//  Copyright (c) 2014 Marcus Johansson.
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#ifndef __MSYM__DEBUG_h
#define __MSYM__DEBUG_h

#include "msym.h"

#ifdef LIBMSYM_DEBUG
#include <stdlib.h>
#include "permutation.h"

#define clean_debug_printf(...) fprintf (stderr, __VA_ARGS__)
#define debug_printTransform(R,C,M) do { printTransform((R),(C),(M)); } while(0)
#define debug_printSubspace(C,L,S) do { printSubspace((C),(L),(S)); } while(0)
#define debug_printPermutation(P) do { printPermutation((P)); } while(0)
#define debug_printCharacterTable(C) do { printCharacterTable((C)); } while(0)


void printTransform(int r, int c, double M[r][c]);
void printPermutation(msym_permutation_t *perm);
void printSubspace(msym_character_table_t *ct, int l, msym_subrepresentation_space_t srs[l]);
void printCharacterTable(msym_character_table_t *ct);

#else
#define clean_debug_printf(fmt, ...) do{} while(0)
#define debug_printTransform(R,C,M) do {} while(0)
#define debug_printSubspace(C,L,S) do {} while(0)
#define debug_printPermutation(P) do {} while(0)
#define debug_printCharacterTable(C) do {} while(0)


#endif

#endif /* defined(__MSYM__DEBUG_h) */
