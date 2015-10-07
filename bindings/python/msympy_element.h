//
//  msympy_element.h
//  libmsympy
//
//  Created by Marcus Johansson on 07/10/15.
//  Copyright (c) 2015 Marcus Johansson.
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include <Python.h>
#include "libmsym/msym.h"

void initmsympy_element(PyObject *m);
int initmsympy_element_ready();
int msympy_elements_from_sequence(PyObject *data, int *out_size, msym_element_t **out_elements);
PyObject *msympy_elements_to_PyList(int len, msym_element_t *elements);
