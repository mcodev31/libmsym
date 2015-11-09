//
//  context.h
//  libmsym
//
//  Created by Marcus Johansson on 30/01/15.
//  Copyright (c) 2015 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#ifndef __MSYM__CONTEXT_h
#define __MSYM__CONTEXT_h

#include <stdio.h>
#include "msym.h"
#include "point_group.h"
#include "basis_function.h"

#define DEFAULT_ZERO_THRESHOLD 1.0e-3
#define DEFAULT_GEOMETRY_THRESHOLD 1.0e-3
#define DEFAULT_ANGLE_THRESHOLD 1.0e-3
#define DEFAULT_EQUIVALENCE_THRESHOLD 5.0e-4
#define DEFAULT_EIGFACT_THRESHOLD 1.0e-3
#define DEFAULT_PERMUTATION_THRESHOLD 5.0e-3
#define DEFAULT_ORTHOGONALIZATION_THRESHOLD 1.0e-2


msym_error_t ctxGetThresholds(msym_context ctx, msym_thresholds_t **thresholds);
msym_error_t ctxSetElements(msym_context ctx, int length, msym_element_t *elements);
msym_error_t ctxGetElements(msym_context, int *l, msym_element_t **elements);
msym_error_t ctxGetExternalElements(msym_context ctx, int *l, msym_element_t **elements);
msym_error_t ctxUpdateExternalElementCoordinates(msym_context ctx);
msym_error_t ctxGetElementPtrs(msym_context ctx, int *l, msym_element_t ***pelements);
msym_error_t ctxGetInternalElement(msym_context ctx, msym_element_t *ext, msym_element_t **element);
msym_error_t ctxGetInternalSubgroup(msym_context ctx, msym_subgroup_t *ext, msym_subgroup_t **sg);
msym_error_t ctxSetPointGroup(msym_context ctx, msym_point_group_t *pg);
msym_error_t ctxGetPointGroup(msym_context ctx, msym_point_group_t **pg);
msym_error_t ctxReduceLinearPointGroup(msym_context ctx);
msym_error_t ctxSetSubgroups(msym_context ctx, int sgl, msym_subgroup_t *sg);
msym_error_t ctxGetSubgroups(msym_context ctx, int *sgl, msym_subgroup_t **sg);
msym_error_t ctxSetEquivalenceSets(msym_context ctx, int esl, msym_equivalence_set_t *es);
msym_error_t ctxGetEquivalenceSets(msym_context ctx, int *esl, msym_equivalence_set_t **es);
msym_error_t ctxGetExternalEquivalenceSets(msym_context ctx, int *esl, msym_equivalence_set_t **es);
msym_error_t ctxGetExternalElementEquivalenceSetMap(msym_context ctx, msym_equivalence_set_t ***eesmap);
msym_error_t ctxSetEquivalenceSetPermutations(msym_context ctx, int r, int c, msym_permutation_t **perm);
msym_error_t ctxGetEquivalenceSetPermutations(msym_context ctx, int *r, int *c, msym_permutation_t ***perm);
msym_error_t ctxGetBasisFunctions(msym_context ctx, int *l, msym_basis_function_t **basis);
msym_error_t ctxGetSubrepresentationSpaces(msym_context ctx, int *srsl, msym_subrepresentation_space_t **srs, int **span);
msym_error_t ctxSetSubrepresentationSpaces(msym_context ctx, int srsl, msym_subrepresentation_space_t *srs, msym_basis_function_t **srsbf, int *span);
msym_error_t ctxUpdateGeometry(msym_context ctx);
msym_error_t ctxGetGeometry(msym_context ctx, msym_geometry_t *g, double eigval[3], double eigvec[3][3]);
msym_error_t ctxDestroyElements(msym_context ctx);
msym_error_t ctxDestroyEquivalcenceSets(msym_context ctx);
msym_error_t ctxDestroyEquivalcenceSetPermutations(msym_context ctx);
msym_error_t ctxDestroyPointGroup(msym_context ctx);
msym_error_t ctxDestroySubgroups(msym_context ctx);
msym_error_t ctxDestroyBasisFunctions(msym_context ctx);
msym_error_t ctxDestroyOrbitalSubspaces(msym_context ctx);
msym_error_t ctxDestroyDisplacementSubspaces(msym_context ctx);
msym_error_t ctxDestroySubrepresentationSpaces(msym_context ctx);



#endif /* defined(__MSYM__CONTEXT_h) */
