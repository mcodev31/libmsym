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
#include <stdlib.h>

#ifdef __LIBMSYM_DEBUG__

#include "permutation.h"

#define dbg_printf(...) fprintf (stderr, __VA_ARGS__)
#define dbg_print_permutation(P) do { printPermutation((P)); } while(0)

void printPermutation(msym_permutation_t *perm);

#ifndef __LIBMSYM_NO_VLA__

#define dbg_print_matrix(R,C,M) do { printMatrix((R),(C),(M)); } while(0)
#define dbg_print_space(C,L,S) do { printSubspace((C),(L),(S)); } while(0)
#define dbg_print_character_table(C) do { printCharacterTable((C)); } while(0)
void printSubspace(msym_character_table_t *ct, int l, msym_subrepresentation_space_t srs[l]);
void printCharacterTable(msym_character_table_t *ct);
void printMatrix(int r, int c, double M[r][c]);

#else
#define dbg_print_matrix(R,C,M) do {} while(0)
#define dbg_print_space(C,L,S) do {} while(0)
#define dbg_print_character_table(C) do {} while(0)
#endif /* ifndef __LIBMSYM_NO_VLA__ */

#else

#define dbg_printf(fmt, ...) do{} while(0)
#define dbg_print_matrix(R,C,M) do {} while(0)
#define dbg_print_space(C,L,S) do {} while(0)
#define dbg_print_permutation(P) do {} while(0)
#define dbg_print_character_table(C) do {} while(0)

#endif /* ifdef __LIBMSYM_DEBUG__ */

#endif /* defined(__MSYM__DEBUG_h) */
