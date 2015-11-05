//
//  msym_error.c
//  libmsym
//
//  Created by Marcus Johansson on 30/01/15.
//  Copyright (c) 2015 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include <stdio.h>
#include <stdarg.h>
#include "msym_error.h"


#define MSYM_ERROR_DETAILS_MAX_LENGTH 1024

const char * invalid = "Invalid error code";

char err_details[MSYM_ERROR_DETAILS_MAX_LENGTH];
char err_details_ext[MSYM_ERROR_DETAILS_MAX_LENGTH];

const struct _errordesc {
    msym_error_t code;
    char *message;
} error_desc[] = {
    { MSYM_SUCCESS, "Success" },
    { MSYM_INVALID_INPUT, "Invalid input" },
    { MSYM_INVALID_CONTEXT, "Invalid context" },
    { MSYM_INVALID_THRESHOLD, "Invalid threshold" },
    { MSYM_INVALID_ELEMENTS, "Invalid elements" },
    { MSYM_INVALID_BASIS_FUNCTIONS, "Invalid basis functions" },
    { MSYM_INVALID_POINT_GROUP, "Invalid point group" },
    { MSYM_INVALID_PERMUTATION, "Invalid permutation" },
    { MSYM_INVALID_EQUIVALENCE_SET, "Invalid equivalence set" },
    { MSYM_INVALID_GEOMETRY, "Invalid geometry" },
    { MSYM_INVALID_CHARACTER_TABLE, "Invalid character table" },
    { MSYM_INVALID_SUBSPACE, "Invalid subspace" },
    { MSYM_INVALID_SUBGROUPS, "Invalid subgroups" },
    { MSYM_INVALID_AXES, "Invalid axes" },
    { MSYM_SYMMETRY_ERROR, "Error determining symmetry operations" },
    { MSYM_PERMUTATION_ERROR, "Error determining permutation" },
    { MSYM_POINT_GROUP_ERROR, "Error determining point group" },
    { MSYM_SYMMETRIZATION_ERROR, "Error symmetrizing molecule/orbtials" },
    { MSYM_SUBSPACE_ERROR, "Error generating subspaces" },
    { MSYM_MEMORY_ERROR, "Error allocating memory" }
};

void msymSetErrorDetails(const char *format, ...){
    va_list args;
    va_start(args, format);
    vsnprintf(err_details, sizeof(err_details), format, args);
    va_end(args);
}

const char MSYM_EXPORT *msymGetErrorDetails(){
    snprintf(err_details_ext, sizeof(err_details_ext), "%s",err_details); // Not really neccessary
    msymSetErrorDetails("");
    return err_details_ext;
}

const char MSYM_EXPORT *msymErrorString(msym_error_t error){
    const char *ret = invalid;
    int length = sizeof(error_desc) / sizeof(error_desc[0]);
    for(int i = 0; i < length;i++){
        if(error == error_desc[i].code){
            ret = error_desc[i].message;
            break;
        }
    }
    return ret;
}
