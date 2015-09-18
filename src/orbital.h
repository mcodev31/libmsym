//
//  orbital.h
//  libmsym
//
//  Created by Marcus Johansson on 07/11/14.
//  Copyright (c) 2014 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#ifndef __MSYM__ORBITAL_h
#define __MSYM__ORBITAL_h

#include <stdio.h>
#include "msym.h"
#include "symop.h"
#include "character_table.h"
#include "permutation.h"
#include "point_group.h"

static double spolynomial [][1] = {
    [0] = {1}
};

static double ppolynomial [][3] = {
    [0] = {0,1,0}, //Py 
    [1] = {0,0,1}, //Pz
    [2] = {1,0,0}  //Px
};


#define M_1_2SQRT3 (1/(2*1.732050807568877293527446341505872366942805253810380628))

static double dpolynomial [][9] = {
    //      x^2         xy   xz   yx    y^2         yz   zx   zy  z^2
    [0] = { 0.0,        0.5, 0.0, 0.5,  0.0,        0.0, 0.0, 0.0, 0.0           },  //d2- (xy)
    [1] = { 0.0,        0.0, 0.0, 0.0,  0.0,        0.5, 0.0, 0.5, 0.0           },  //d1- (yz)
    [2] = {-M_1_2SQRT3, 0.0, 0.0, 0.0, -M_1_2SQRT3, 0.0, 0.0, 0.0, 2.0*M_1_2SQRT3},  //d0 (z^2)
    [3] = { 0.0,        0.0, 0.5, 0.0,  0.0,        0.0, 0.5, 0.0, 0.0           },  //d1+ (xz)
    [4] = { 0.5,        0.0, 0.0, 0.0, -0.5,        0.0, 0.0, 0.0, 0.0           }   //d2+ (x^2-y^2)
};


msym_error_t basisFunctionFromName(char *, msym_basis_function_t *bf);
msym_error_t basisFunctionFromQuantumNumbers(int n, int l, int m, msym_basis_function_t *bf);

msym_error_t generateBasisFunctionTransforms(int sopsl, msym_symmetry_operation_t sops[sopsl], int l, double transform[sopsl][2*l+1][2*l+1]);
msym_error_t generateSALCSubspaces(msym_point_group_t *pg, int sgl, msym_subgroup_t sg[sgl], int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, int basisl, msym_basis_function_t basis[basisl], msym_element_t *elements, msym_equivalence_set_t **esmap, msym_thresholds_t *thresholds, int *ossl, msym_subspace_t **oss, int **ospan);


#endif /* defined(__MSYM__ORBITAL_h) */
