//
//  msympy_context.c
//  libmsympy
//
//  Created by Marcus Johansson on 07/10/15.
//  Copyright (c) 2015 Marcus Johansson.
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include <Python.h>
#include "structmember.h"
#include "msympy_context.h"
#include "msympy_element.h"
#include "libmsym/msym.h"

#define DEBUG_FUNCTION() printf("%s:%d\n",__FUNCTION__,__LINE__)

typedef struct {
    PyObject_HEAD
    msym_context ctx;
    PyObject *elements;
} MSymPyContext;

static PyObject *MSymPyContextError;

// Initialization
static void MSymPyContext_dealloc(MSymPyContext* self);
static PyObject *MSymPyContext_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static int MSymPyContext_init(MSymPyContext *self, PyObject *args, PyObject *kwds);

// Getters & Setters
static PyObject *MSymPyContext_getElements(MSymPyContext *self, void *closure);
static int MSymPyContext_setElements(MSymPyContext *self, PyObject *data, void *closure);

// Methods
static PyObject *MSymPyContext_findSymmetry(MSymPyContext* self);

static PyMethodDef MSymPyContext_methods[] = {
    {"findSymmetry", (PyCFunction)MSymPyContext_findSymmetry, METH_NOARGS, "determine symmetry of current elements"},
    {NULL}  /* Sentinel */
};

static PyGetSetDef MSymPyContext_getseters[] = {
    {"elements", (getter)MSymPyContext_getElements, (setter)MSymPyContext_setElements, "context elements", NULL},
    {NULL}  /* Sentinel */
};

static PyTypeObject MSymPyContextType = {
    PyObject_HEAD_INIT(NULL)
    0,                                  /*ob_size*/
    "msympy.Context",                   /*tp_name*/
    sizeof(MSymPyContext),              /*tp_basicsize*/
    0,                                  /*tp_itemsize*/
    (destructor)MSymPyContext_dealloc,  /*tp_dealloc*/
    0,                                  /*tp_print*/
    0,                                  /*tp_getattr*/
    0,                                  /*tp_setattr*/
    0,                                  /*tp_compare*/
    0,                                  /*tp_repr*/
    0,                                  /*tp_as_number*/
    0,                                  /*tp_as_sequence*/
    0,                                  /*tp_as_mapping*/
    0,                                  /*tp_hash */
    0,                                  /*tp_call*/
    0,                                  /*tp_str*/
    0,                                  /*tp_getattro*/
    0,                                  /*tp_setattro*/
    0,                                  /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "MSymPyContext objects",            /* tp_doc */
    0,                                  /* tp_traverse */
    0,                                  /* tp_clear */
    0,                                  /* tp_richcompare */
    0,                                  /* tp_weaklistoffset */
    0,                                  /* tp_iter */
    0,                                  /* tp_iternext */
    MSymPyContext_methods,              /* tp_methods */
    0,
    MSymPyContext_getseters,            /* tp_getset */
    0,                                  /* tp_base */
    0,                                  /* tp_dict */
    0,                                  /* tp_descr_get */
    0,                                  /* tp_descr_set */
    0,                                  /* tp_dictoffset */
    (initproc)MSymPyContext_init,       /* tp_init */
    0,                                  /* tp_alloc */
    MSymPyContext_new,                  /* tp_new */
};

void initmsympy_context(PyObject *m){
    DEBUG_FUNCTION();
    Py_INCREF(&MSymPyContextType);
    PyModule_AddObject(m, "Context", (PyObject *)&MSymPyContextType);
    MSymPyContextError = PyErr_NewException("msympy.error", NULL, NULL);
    Py_INCREF(MSymPyContextError);
    PyModule_AddObject(m, "error", MSymPyContextError);
}

int initmsympy_context_ready(){
    DEBUG_FUNCTION();
    return PyType_Ready(&MSymPyContextType);
}

static void MSymPyContext_dealloc(MSymPyContext* self){
    DEBUG_FUNCTION();
    msymReleaseContext(self->ctx);
    Py_XDECREF(self->elements);
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *MSymPyContext_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
    DEBUG_FUNCTION();
    MSymPyContext *self;
    
    self = (MSymPyContext *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->ctx = msymCreateContext();
        if (self->ctx == NULL){
            Py_DECREF(self);
            self = NULL;
        } else {
            self->elements = Py_BuildValue("()");
        }
    }
    
    return (PyObject *)self;
}

static int MSymPyContext_init(MSymPyContext *self, PyObject *args, PyObject *kwds){
    DEBUG_FUNCTION();
    return 0;
}



static PyObject *MSymPyContext_getElements(MSymPyContext *self, void *closure){
    //Make new elements when we call symmetrize and generate, otherwise just keep what we have
    // generate needs set element id function in msym
    //Py_INCREF(self->elements);
    return NULL;
}

static int MSymPyContext_setElements(MSymPyContext *self, PyObject *data, void *closure)
{
    msym_element_t *elements = NULL;
    int elementsl = 0;
    //const msym_element_t *melements = NULL;
    //int melementsl = 0;
    //PyObject **eobjs = NULL;
    
    if(msympy_elements_from_sequence(data, &elementsl, &elements) < 0){
        return -1;
    }
    
    /* might not need to decrease the reference counter since there is one in the elements sequence
    if(MSYM_SUCCESS == msymGetElements(self->ctx,&melementsl,&melements)){
        eobjs = (PyObject**) PyMem_Malloc(melementsl * sizeof(PyObject *));
        for(int i = 0;i < melementsl;i++){
            eobjs[i] = (PyObject *) melements[i].id;
        }
    }*/
    
    if(MSYM_SUCCESS != msymSetElements(self->ctx,elementsl,elements)){
        PyMem_Free(elements);
        //PyMem_Free(eobjs);
        PyErr_SetString(PyExc_ValueError, msymGetErrorDetails());
        return -1;
    }
    
    /* might not need to increase the reference counter since there is one in the elements sequence
    for(int i = 0;i < elementsl;i++) {
        Py_INCREF((PyObject *) elements[i].id);
    }
    
    for(int i = 0;i < melementsl;i++){
        Py_DECREF(eobjs[i]);
    }*/
    
    PyMem_Free(elements);
    
    Py_INCREF(data);
    Py_XDECREF(self->elements);
    self->elements = data;
    
    return 0;
}

static PyObject *MSymPyContext_findSymmetry(MSymPyContext* self){
    DEBUG_FUNCTION();
    msym_error_t ret = MSYM_SUCCESS;
    char buf[16];
    if(MSYM_SUCCESS != (ret = msymFindSymmetry(self->ctx))){
        PyErr_SetString(MSymPyContextError, msymGetErrorDetails());
        return NULL;
    }
    if(MSYM_SUCCESS != (ret = msymGetPointGroupName(self->ctx, sizeof(buf),buf))){
        PyErr_SetString(MSymPyContextError, msymGetErrorDetails());
        return NULL;
    }
    PyObject *name = PyString_FromString(buf);
    return name;
}
