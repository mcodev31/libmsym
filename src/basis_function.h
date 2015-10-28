//
//  basis_function.h
//  libmsym
//
//  Created by Marcus Johansson on 07/11/14.
//  Copyright (c) 2014 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#ifndef __MSYM__BASIS_FUNCTION_h
#define __MSYM__BASIS_FUNCTION_h

#include <stdio.h>
#include "msym.h"

msym_error_t basisFunctionFromName(char *, msym_basis_function_t *bf);
msym_error_t basisFunctionFromQuantumNumbers(int n, int l, int m, msym_basis_function_t *bf);


#endif /* defined(__MSYM__BASIS_FUNCTION_h) */
