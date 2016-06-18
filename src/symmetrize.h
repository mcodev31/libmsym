//
//  symmetrize.h
//  Symmetry
//
//  Created by Marcus Johansson on 04/02/15.
//  Copyright (c) 2015 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#ifndef __MSYM__SYMMETRIZE_h
#define __MSYM__SYMMETRIZE_h

#include <stdio.h>

#include "msym.h"
#include "point_group.h"
#include "permutation.h"

msym_error_t symmetrizeElements(msym_point_group_t *pg, int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, msym_thresholds_t *thresholds, double *err);
msym_error_t symmetrizeTranslation(msym_point_group_t *pg, msym_equivalence_set_t *es, msym_permutation_t *perm, int pi, double translation[3]);
#ifndef __LIBMSYM_NO_VLA__
msym_error_t symmetrizeWavefunctions(msym_point_group_t *pg, int srsl, msym_subrepresentation_space_t *srs, int *span, int basisl, msym_basis_function_t *basis, double (*wf)[basisl], double (*symwf)[basisl], int *species, msym_partner_function_t *pfo);
#else
msym_error_t symmetrizeWavefunctions(msym_point_group_t *pg, int srsl, msym_subrepresentation_space_t *srs, int *span, int basisl, msym_basis_function_t *basis, void *wf, void *symwf, int *species, msym_partner_function_t *pfo);
#endif /* fndef __LIBMSYM_NO_VLA__ */
#endif /* defined(__MSYM__SYMMETRIZE_h) */
