//
//  geometry.h
//  libmsym
//
//  Created by Marcus Johansson on 28/11/14.
//  Copyright (c) 2014 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#ifndef __MSYM__GEOMETRY_h
#define __MSYM__GEOMETRY_h

#include <stdio.h>
#include "msym.h"

enum geometry {
    GEOMETRY_UNKNOWN = -1,
    SPHERICAL,
    LINEAR,
    PLANAR_REGULAR,
    PLANAR_IRREGULAR,
    POLYHEDRAL_PROLATE,
    POLYHEDRAL_OBLATE,
    ASSYMETRIC
};


typedef enum geometry geometry_t;

enum geometry find_geometry(double e[3]);
msym_error_t findGeometry(int length, msym_element_t *elements[length], double cm[3], msym_thresholds_t *thresholds, geometry_t *g, double v[3][3]);
msym_error_t findCenterOfMass(int length, msym_element_t *elements[length], double v[3]);
void printGeometry(geometry_t g);

#endif /* defined(__MSYM__GEOMETRY_h) */
