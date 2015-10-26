//
//  example.h
//  libmsym
//
//  Created by Marcus Johansson on 24/04/15.
//  Copyright (c) 2015 Marcus Johansson.
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

// Use <libmsym/msym.h> if installed
#include "msym.h"

int example(const char* in_file, msym_thresholds_t *thresholds);
int read_xyz(const char *name, msym_element_t **ratoms);
