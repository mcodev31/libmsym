//
//  msympymodule.c
//  libmsympy
//
//  Created by Marcus Johansson on 07/10/15.
//  Copyright (c) 2015 Marcus Johansson.
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include <Python.h>
#include "libmsym/msym.h"

#include "msympy_context.h"
#include "msympy_element.h"

#define DEBUG_FUNCTION() printf("%s:%d\n",__FUNCTION__,__LINE__)

static PyObject *msympy_getDefaultThresholds(PyObject *self, PyObject *args);

static PyMethodDef MSymPy_methods[] = {
    {"getDefaultThresholds", msympy_getDefaultThresholds, METH_VARARGS, "Get default thresholds"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static PyObject *msympy_getDefaultThresholds(PyObject *self, PyObject *args){
    return NULL;
}

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC initmsympy(void) {
    printf("hello\n");
    DEBUG_FUNCTION();
    PyObject* m;
    
    if(initmsympy_context_ready() < 0 ||
       initmsympy_element_ready() < 0)
        return;
    
    m = Py_InitModule3("msympy", MSymPy_methods,
                       "libmsym python module");
    
    if (m == NULL)
        return;
    
    initmsympy_context(m);
    initmsympy_element(m);
}