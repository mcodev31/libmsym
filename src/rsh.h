//
//  rsh.h
//  libmsym
//
//  Created by Marcus Johansson on 21/10/15.
//  Copyright (c) 2015 Marcus Johansson.
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include "msym.h"

typedef struct _rsh_representations {
    int d;
    void *t;
} rsh_representations_t;

msym_error_t generateRSHRepresentations(int sopsl, msym_symmetry_operation_t sops[sopsl], int lmax, rsh_representations_t *lrs);
