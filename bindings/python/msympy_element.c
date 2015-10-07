//
//  msympy_element.c
//  libmsympy
//
//  Created by Marcus Johansson on 07/10/15.
//  Copyright (c) 2015 Marcus Johansson.
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include <Python.h>
#include "structmember.h"
#include "msympy_element.h"
#include "libmsym/msym.h"

#define DEBUG_FUNCTION() printf("%s:%d\n",__FUNCTION__,__LINE__)

typedef struct {
    PyObject_HEAD
    double m;
    double v[3];
    int n;
    PyObject *name;
} MSymPyElement;

static void MSymPyElement_dealloc(MSymPyElement* self);
static PyObject *MSymPyElement_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static int MSymPyElement_init(MSymPyElement *self, PyObject *args, PyObject *kwds);

// Getters & Setters
static PyObject *MSymPyElement_getCoordinates(MSymPyElement *self, void *closure);
static int MSymPyElement_setCoordinates(MSymPyElement *self, PyObject *data, void *closure);
static PyObject *MSymPyElement_getName(MSymPyElement *self, void *closure);
static int MSymPyElement_setName(MSymPyElement *self, PyObject *data, void *closure);

static PyMemberDef MSymPyElement_members[] = {
    {"charge", T_INT, offsetof(MSymPyElement, n), 0, "element nuclear charge"},
    {"mass", T_DOUBLE, offsetof(MSymPyElement, m), 0, "element mass"},
    {NULL}  /* Sentinel */
};

static PyGetSetDef MSymPyElement_getseters[] = {
    {"coordinates", (getter)MSymPyElement_getCoordinates, (setter)MSymPyElement_setCoordinates, "element coordinates", NULL},
    {"name", (getter)MSymPyElement_getName, (setter)MSymPyElement_setName, "element name", NULL},
    {NULL}  /* Sentinel */
};

static PyTypeObject MSymPyElementType = {
    PyObject_HEAD_INIT(NULL)
    0,                                  /*ob_size*/
    "msympy.Element",                   /*tp_name*/
    sizeof(MSymPyElement),              /*tp_basicsize*/
    0,                                  /*tp_itemsize*/
    (destructor)MSymPyElement_dealloc,  /*tp_dealloc*/
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
    "MSymPyElement objects",            /* tp_doc */
    0,                                  /* tp_traverse */
    0,                                  /* tp_clear */
    0,                                  /* tp_richcompare */
    0,                                  /* tp_weaklistoffset */
    0,                                  /* tp_iter */
    0,                                  /* tp_iternext */
    0,                                  /* tp_methods */
    MSymPyElement_members,              /* tp_members */
    MSymPyElement_getseters,            /* tp_getset */
    0,                                  /* tp_base */
    0,                                  /* tp_dict */
    0,                                  /* tp_descr_get */
    0,                                  /* tp_descr_set */
    0,                                  /* tp_dictoffset */
    (initproc)MSymPyElement_init,       /* tp_init */
    0,                                  /* tp_alloc */
    MSymPyElement_new,                  /* tp_new */
};

/*********************
 * Utility functions *
 *********************/

void initmsympy_element(PyObject *m){
    DEBUG_FUNCTION();
    Py_INCREF(&MSymPyElementType);
    PyModule_AddObject(m, "Element", (PyObject *)&MSymPyElementType);
}

int initmsympy_element_ready(){
    DEBUG_FUNCTION();
    return PyType_Ready(&MSymPyElementType);
}

int msympy_elements_from_sequence(PyObject *data, int *out_size, msym_element_t **out_elements){
    
    msym_element_t *elements;
    int size = 0;
    PyObject* seq = NULL;
    
    seq = PySequence_Fast(data, "expected a sequence");
    
    if(NULL == seq){
        return -1;
    }
    
    size = PySequence_Size(seq);
    if (size < 0){
        PyErr_SetString(PyExc_TypeError, "expected a sequence of size > 0");
        return -1;
    }
    
    elements = (msym_element_t*) PyMem_Malloc(size * sizeof(msym_element_t));
    
    if (NULL == elements) {
        Py_DECREF(seq);
        PyErr_NoMemory();
        return -1;
    }
    
    for (int i = 0; i < size; i++){
        PyObject *obj = PySequence_Fast_GET_ITEM(seq, i);
        if(NULL == obj || !PyObject_TypeCheck(obj, &MSymPyElementType)){
            PyErr_SetString(PyExc_TypeError, "expected a sequence of msympy.Element");
            PyMem_Free(elements);
            Py_DECREF(seq);
            return -1;
        }
        MSymPyElement *pyelement = (MSymPyElement *) obj;
        
        elements[i].id = obj;
        elements[i].n = pyelement->n;
        elements[i].m = pyelement->m;
        elements[i].v[0] = pyelement->v[0];
        elements[i].v[1] = pyelement->v[1];
        elements[i].v[2] = pyelement->v[2];
        snprintf(elements[i].name,sizeof(elements[i].name), "%s", PyString_AsString(pyelement->name));
    }
    
    Py_DECREF(seq);
    
    if (PyErr_Occurred()) {
        return -1;
    }
    
    *out_elements = elements;
    *out_size = size;
    
    return 0;
}

PyObject *msympy_elements_to_PyList(int len, msym_element_t *elements){
    PyObject* list = PyList_New((Py_ssize_t) len);
    if(NULL == list){
        return NULL;
    }

    for(int i = 0;i < len;i++){
        MSymPyElement *pyelement = PyObject_New(MSymPyElement, &MSymPyElementType);
        
        if(NULL == pyelement || PyList_SetItem(list, (Py_ssize_t) i, (PyObject*) pyelement) < 0){
            Py_XDECREF(list);
            list = NULL;
            break;
        }
        
        pyelement->n = elements[i].n;
        pyelement->m = elements[i].m;
        pyelement->v[0] = elements[i].v[0];
        pyelement->v[1] = elements[i].v[1];
        pyelement->v[2] = elements[i].v[2];
        pyelement->name = PyString_FromString(elements[i].name);
    }
    
    return list;
}




/********************
 * Python functions *
 ********************/


static void MSymPyElement_dealloc(MSymPyElement* self){
    DEBUG_FUNCTION();
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *MSymPyElement_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
    DEBUG_FUNCTION();
    MSymPyElement *self = (MSymPyElement *)type->tp_alloc(type, 0);
    
    if (self != NULL) {
        self->name = PyString_FromString("");
        if(self->name == NULL){
            Py_DECREF(self);
            self = NULL;
        } else {
            self->m = 0.0;
            self->n = 0;
            self->v[0] = self->v[1] = self->v[2] = 0.0;
        }
    }
    
    return (PyObject *)self;
}

static int MSymPyElement_init(MSymPyElement *self, PyObject *args, PyObject *kwds){
    DEBUG_FUNCTION();
    
    PyObject *name = NULL;
    
    static char *kwlist[] = {"name", "mass", "charge", "coordinates", NULL};
    
    if (! PyArg_ParseTupleAndKeywords(args, kwds, "|Sdi(ddd)", kwlist,
                                      &name,
                                      &self->m,
                                      &self->n,
                                      &self->v[0], &self->v[1], &self->v[2]))
        return -1;
    
    if (name) {
        Py_INCREF(name);
        Py_XDECREF(self->name);
        self->name = name;
    }
    
    return 0;
}

static PyObject *MSymPyElement_getName(MSymPyElement *self, void *closure) {
    Py_INCREF(self->name);
    return self->name;
}

static int MSymPyElement_setName(MSymPyElement *self, PyObject *name, void *closure) {
    if (name == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete name");
        return -1;
    }
    
    if (!PyString_Check(name)) {
        PyErr_SetString(PyExc_TypeError, "Name must be a string");
        return -1;
    }
    
    Py_DECREF(self->name);
    Py_INCREF(name);
    self->name = name;
    
    return 0;
}

static PyObject *MSymPyElement_getCoordinates(MSymPyElement *self, void *closure) {
    return Py_BuildValue("(ddd)",self->v[0],self->v[1],self->v[2]);
}

static int MSymPyElement_setCoordinates(MSymPyElement *self, PyObject *value, void *closure) {
    PyObject* seq = 0;
    int size = 0;
    double v[3];
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete coordinates");
        return -1;
    }
    
    seq = PySequence_Fast(value, "expected a sequence of double");
    if (!seq)
        return -1;
    
    size = PySequence_Size(seq);
    if (size != 3){
        PyErr_SetString(PyExc_ValueError, "expected a sequence of 3 double");
        return -1;
    }
    
    for (int i = 0; i < 3; i++) v[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
    
    Py_DECREF(seq);
    
    if (PyErr_Occurred()) {
        return -1;
    }
    
    for (int i = 0; i < 3; i++) self->v[i] = v[i];
    

    
    /*if (!PyArg_ParseTuple(value, "(ddd)", &self->v[0],&self->v[1],&self->v[2])){
        PyErr_SetString(PyExc_TypeError, "Coordinates must be sequence of 3 double");
        return -1;
    }*/
    
    return 0;
}