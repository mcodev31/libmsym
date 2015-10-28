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
#include "point_group.h"

void freeSubrepresentationSpaces(int srsl, msym_subrepresentation_space_t *srs);
msym_error_t generateSubrepresentationSpaces(msym_point_group_t *pg, int sgl, const msym_subgroup_t sg[sgl], int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, int basisl, msym_basis_function_t basis[basisl], msym_element_t *elements, msym_equivalence_set_t **esmap, msym_thresholds_t *thresholds, int *osrsl, msym_subrepresentation_space_t **osrs, msym_basis_function_t ***osrsbf, int **ospan);
msym_error_t symmetrySpeciesComponents(msym_point_group_t *pg, int srsl, msym_subrepresentation_space_t *srs, int basisl, msym_basis_function_t *basis, double *wf, double *s);

#endif /* defined(__MSYM__SUBSPACE_h) */

