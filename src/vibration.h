//
//  vibration.h
//  libmsym
//
//  Created by Marcus Johansson on 26/05/15.
//  Copyright (c) 2014 Marcus Johansson.
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#ifndef __MSYM__VIBRATION_h
#define __MSYM__VIBRATION_h

#include <stdio.h>
#include "msym.h"
#include "symop.h"
#include "character_table.h"
#include "permutation.h"
#include "point_group.h"

msym_error_t generateDisplacementSubspaces(msym_point_group_t *pg, int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, msym_thresholds_t *thresholds, int *subspacel, msym_subspace_t **subspace, int **pspan);

#endif /* defined(__MSYM__VIBRATION_h) */