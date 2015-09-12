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

#include "msym.h"
#include "symop.h"
#include "geometry.h"
#include "character_table.h"
#include "permutation.h"


msym_error_t findPointGroup(int sopsl, msym_symmetry_operation_t *sops, msym_thresholds_t *thresholds, msym_point_group_t **pg);
msym_error_t findSubgroup(msym_subgroup_t *subgroup, msym_thresholds_t *thresholds);
msym_error_t findCharacterTable(msym_point_group_t *pg);
msym_error_t generatePointGroupFromName(const char *name, msym_thresholds_t *thresholds, msym_point_group_t **opg);
msym_error_t pointGroupFromSubgroup(msym_subgroup_t *sg, msym_thresholds_t *thresholds, msym_point_group_t **opg);
int numberOfSubgroups(msym_point_group_t *pg);

#endif /* defined(__MSYM__POINT_GROUP_h) */
