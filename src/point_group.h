//
//  point_group.h
//  Symmetry
//
//  Created by Marcus Johansson on 28/11/14.
//  Copyright (c) 2014 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#ifndef __MSYM__POINT_GROUP_h
#define __MSYM__POINT_GROUP_h

#include "symop.h"
#include "geometry.h"
#include "character_table.h"
#include "permutation.h"

enum pg_type {
    POINT_GROUP_Ci,
    POINT_GROUP_Cs,
    POINT_GROUP_Cn,
    POINT_GROUP_Cnh,
    POINT_GROUP_Cnv,
    POINT_GROUP_Dn,
    POINT_GROUP_Dnh,
    POINT_GROUP_Dnd,
    POINT_GROUP_S2n,
    POINT_GROUP_T,
    POINT_GROUP_Td,
    POINT_GROUP_Th,
    POINT_GROUP_O,
    POINT_GROUP_Oh,
    POINT_GROUP_I,
    POINT_GROUP_Ih,
    POINT_GROUP_K,
    POINT_GROUP_Kh
};

//We cant handle names larger than e.g. D999v (should be more than sufficient, we can see the order in the group anyways)
typedef struct {
    enum pg_type type;
    int n;
    int order;
    msym_symmetry_operation_t *primary;
    msym_symmetry_operation_t *sops;
    msym_permutation_t *perm;
    int sopsl;
    double transform[3][3];
    CharacterTable *ct;
    char name[6];
} msym_point_group_t;

msym_error_t findPointGroup(int sopsl, msym_symmetry_operation_t *sops, enum geometry g, msym_thresholds_t *thresholds, msym_point_group_t **pg);
msym_error_t findCharacterTable(msym_point_group_t *pg);
msym_error_t generatePointGroup(char *name, msym_thresholds_t *thresholds, msym_point_group_t **opg);

#endif /* defined(__MSYM__POINT_GROUP_h) */
