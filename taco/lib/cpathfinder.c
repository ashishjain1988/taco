//
//  cpathfinder.c
//  TACO
//
//  Created by Balaji Pandian on 4/2/16.
//  Copyright © 2016 Balaji Pandian. All rights reserved.
//
//  TODO:
//  

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

#include <Python.h>

// Constant minimum path score
#define MIN_SCORE (0.0000000001)

typedef struct Path_and_expr_t {
    int test;
} Path_and_expr_t;

Path_and_expr_t* find_path(long* order, int order_len, double* exprs, int exprs_len, PyObject** succs, int succs_len, int source, int sink) {
    double* min_exprs = (double*)malloc(sizeof(double) * exprs_len);
    double* sum_exprs = (double*)malloc(sizeof(double) * exprs_len);
    int* lengths = (int*)malloc(sizeof(int) * exprs_len);
    int* prevs = (int*)malloc(sizeof(int) * exprs_len);

    for (int i = 0; i < exprs_len; i++) {
        min_exprs[i] = MIN_SCORE;
        sum_exprs[i] = MIN_SCORE;
        lengths[i] = 1;
        prevs[i] = sink;
    }
    min_exprs[source] = exprs[source];
    sum_exprs[source] = exprs[source];

    double min_expr;
    double new_min_expr;
    double sum_expr;
    double new_sum_expr;
    int length;

    // Traverse nodes in topological sort order
    for (int i = 0; i < order_len; i++) {
        long cur_order = order[i];
        min_expr = min_exprs[cur_order];
        sum_expr = sum_exprs[cur_order];
        length = lengths[cur_order];

        int set_len = PySet_Size(succs[cur_order]);
        for (int j = 0; j < set_len; j++) {
            long 
            if (min_expr < exprs[])
        }
    }

    return NULL;
}

static PyObject* py_cpathfinder(PyObject* self, PyObject* args) {
    float G_exprs_G_SOURCE_ID;
    PyObject* G_top_sort_list_obj;
    PyObject* G_exprs_list_obj;
    PyObject* G_succs_list_obj;
    float path_frac;
    int max_paths;
    int source;
    int sink;

    if (!PyArg_ParseTuple(args, "fO!O!O!fiii", &G_exprs_G_SOURCE_ID, &PyList_Type, &G_top_sort_list_obj, &PyList_Type, &G_exprs_list_obj, &PyList_Type, &G_succs_list_obj, &path_frac, &max_paths, &source, &sink)) {
        return NULL;
    }

    if (!G_top_sort_list_obj || !G_exprs_list_obj || !G_succs_list_obj) {
        return NULL;
    }

    printf("\n\nEverything is working\n\n");

    // Don't run if all nodes are zero
    if (G_exprs_G_SOURCE_ID < MIN_SCORE) {
        return Py_BuildValue("[]");
    }

    int G_top_sort_list_obj_len = (int)PyList_Size(G_top_sort_list_obj);
    int G_exprs_list_obj_len = (int)PyList_Size(G_exprs_list_obj);
    int G_succs_list_obj_len = (int)PyList_Size(G_succs_list_obj);

    printf("\n");
    printf("G_top_sort_list_obj_len: %d\n", G_top_sort_list_obj_len);
    printf("G_exprs_list_obj_len: %d\n", G_exprs_list_obj_len);
    printf("G_succs_list_obj_len: %d\n", G_succs_list_obj_len);


    long* order = (long*)calloc(G_top_sort_list_obj_len, sizeof(long));
    double* exprs = (double*)calloc(G_exprs_list_obj_len, sizeof(double));
    PyObject** succs = (PyObject**)calloc(G_succs_list_obj_len, sizeof(PyObject*));

    for (int i = 0; i < G_top_sort_list_obj_len; i++) {
        order[i] = PyInt_AsLong(PyList_GetItem(G_top_sort_list_obj, i));
    }
    for (int i = 0; i < G_exprs_list_obj_len; i++) {
        exprs[i] = PyFloat_AsDouble(PyList_GetItem(G_exprs_list_obj, i));
    }
    for (int i = 0; i < G_succs_list_obj_len; i++) {
        succs[i] = PyList_GetItem(G_succs_list_obj, i);
    }

    find_path(order, G_top_sort_list_obj_len, exprs, G_exprs_list_obj_len, succs, G_succs_list_obj_len, source, sink);

    free(order);
    free(exprs);
    free(succs);
    // Return value -- returning "NULL" is an error
    return Py_BuildValue("i", 0);
}

// Mapping between python and c function names
static PyMethodDef CPathfinderMethods[] = {
    {"cpathfinder", py_cpathfinder, METH_VARARGS},
    {NULL, NULL}
};

// Module initialisation routine
PyMODINIT_FUNC initcpathfinder(void)
{
    (void) Py_InitModule("cpathfinder", CPathfinderMethods);
}
