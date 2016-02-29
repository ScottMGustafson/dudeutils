#include <Python.h>
#include <numpy/arrayobject.h>
#include "spectrum.h"

/* Docstrings */
static char module_docstring[] =
    "";
static char spectrum_docstring[] =
    "";

/* Available functions */
static PyObject *spectrum_spectrum(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = {
    {"spectrum", spectrum_spectrum, METH_VARARGS, spectrum_docstring},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef _spectrum ={
    PyModuleDef_HEAD_INIT,
    "spectrum",
    module_docstring,
    -1,  
    module_methods
};

/* Initialize the module */
PyMODINIT_FUNC PyInit__spectrum(void)
{
    PyObject *m = PyModule_Create(&_spectrum);
    //PyObject *m = Py_InitModule3("_spectrum", module_methods, module_docstring);
    if (m == NULL)
        return;

    /* Load `numpy` functionality. */
    import_array();
}

static PyObject *spectrum_spectrum(PyObject *self, PyObject *args)
{
    PyObject *wave_obj, *flux_obj, *err_obj;
    PyObject *x_obj, *y_obj;
    PyObject *N_obj, *b_obj, *z_obj;
    PyObject *rest_obj, *gamma_obj, *f_obj;
    PyObject *starts_obj, *ends_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOOOOOOOOOOOO",&wave_obj, &flux_obj, &err_obj,
                                                &x_obj, &y_obj, 
                                                &N_obj, &b_obj,&z_obj,
                                                &rest_obj,&gamma_obj,&f_obj,
                                                &starts_obj, &ends_obj))
        return NULL;

    /* Interpret the input objects as numpy arrays. */
    PyObject *wave_arr = PyArray_FROM_OTF(wave_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *flux_arr = PyArray_FROM_OTF(flux_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *err_arr = PyArray_FROM_OTF(err_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    //PyObject *cont_arr = PyArray_FROM_OTF(cont_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *x_arr    = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *y_arr    = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *N_arr    = PyArray_FROM_OTF(N_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *b_arr    = PyArray_FROM_OTF(b_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *z_arr    = PyArray_FROM_OTF(z_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *rest_arr = PyArray_FROM_OTF(rest_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *gamma_arr= PyArray_FROM_OTF(gamma_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *f_arr    = PyArray_FROM_OTF(f_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *starts_arr    = PyArray_FROM_OTF(starts_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *ends_arr    = PyArray_FROM_OTF(ends_obj, NPY_DOUBLE, NPY_IN_ARRAY);


    /* If that didn't work, throw an exception. */
    if (wave_arr == NULL || flux_arr == NULL || x_arr == NULL|| 
            y_arr == NULL|| N_arr == NULL|| b_arr == NULL|| z_arr == NULL|| 
            rest_arr == NULL|| gamma_arr == NULL|| f_arr == NULL
            || starts_arr == NULL|| ends_arr == NULL) {
        Py_XDECREF(wave_arr);
        Py_XDECREF(flux_arr);
        Py_XDECREF(err_arr);
        //Py_XDECREF(cont_arr);
        Py_XDECREF(x_arr);
        Py_XDECREF(y_arr);
        Py_XDECREF(N_arr);
        Py_XDECREF(b_arr);
        Py_XDECREF(z_arr);
        Py_XDECREF(rest_arr);
        Py_XDECREF(gamma_arr);
        Py_XDECREF(f_arr);
        Py_XDECREF(starts_arr);
        Py_XDECREF(ends_arr);
        return NULL;
    }
 /* How many data points are there? */
    int len_cont_points = (int)PyArray_DIM(x_arr, 0);
    int len_arr = (int)PyArray_DIM(wave_arr, 0);
    int len_abs = (int)PyArray_DIM(N_arr, 0);
    int len_pairs = (int)PyArray_DIM(starts_arr, 0);


    /* Get pointers to the data as C-types. */
    const double *wave    = (double*)PyArray_DATA(wave_arr);
    const double *flux    = (double*)PyArray_DATA(flux_arr);
    const double *err    = (double*)PyArray_DATA(err_arr);
    //double *continuum = (double*)PyArray_DATA(cont_arr);
    double *x       = (double*)PyArray_DATA(x_arr);
    double *y       = (double*)PyArray_DATA(y_arr);
    double *N       = (double*)PyArray_DATA(N_arr);
    double *b       = (double*)PyArray_DATA(b_arr);
    double *z       = (double*)PyArray_DATA(z_arr);
    double *rest    = (double*)PyArray_DATA(rest_arr);
    double *gamma   = (double*)PyArray_DATA(gamma_arr);
    double *f       = (double*)PyArray_DATA(f_arr);
    double *starts   = (double*)PyArray_DATA(starts_arr);
    double *ends      = (double*)PyArray_DATA(ends_arr);

    int dims[1],i;//abs_dims[1],i;
    dims[0]=len_arr;
    //abs_dims[0]=len_abs;

    PyArrayObject* contin_out=(PyArrayObject*) PyArray_FromDims(1,dims, NPY_DOUBLE);
    PyArrayObject* abs_out=(PyArrayObject*) PyArray_FromDims(1,dims, NPY_DOUBLE);
    //PyArrayObject* N_out=(PyArrayObject*) PyArray_FromDims(1,len_dims, NPY_DOUBLE);
    //PyArrayObject* b_out=(PyArrayObject*) PyArray_FromDims(1,len_dims, NPY_DOUBLE);
    //PyArrayObject* z_out=(PyArrayObject*) PyArray_FromDims(1,len_dims, NPY_DOUBLE);

/*TODO:  output N,b,z.  implement as outputing a struct of N,b,z*/

    double* cont_data=(double*)contin_out->data;
    double* abs_data=(double*)abs_out->data;
    //double* N_data=(double*)N_out->data;
    //double* b_data=(double*)b_out->data;
    //double* z_data=(double*)z_out->data;

    /*double* value = (double*)spectrum(wave,flux, err,continuum,x, y,N,b,z,rest,gamma,f, starts, ends,len_cont_points,len_arr,len_abs,len_pairs)
    for(i=0;i<len_arr;++i){
        data[i]=value[i];
    }*/

    double* cont = (double*)return_continuum(x,y, wave, len_cont_points, len_arr);
    for(i=0;i<len_arr;++i){
        cont_data[i]=cont[i];
    }

    double* ab = (double*)return_absorption(cont, wave, N,b, z, rest, gamma, f,len_abs, len_arr);
    for(i=0;i<len_arr;++i){
        abs_data[i]=ab[i];
    }


    double chi2 = (double)get_chi2(ab, flux, err, wave, starts, ends, len_pairs, len_arr);

    if (chi2 < 0.0) {
        PyErr_SetString(PyExc_RuntimeError,
                    "Chi-squared returned an impossible value.");
        return NULL;
    }
    /* Build the output tuple */


    /* Clean up. */
    Py_DECREF(wave_arr);
    Py_DECREF(flux_arr);
    Py_DECREF(err_arr);
    Py_DECREF(x_arr);
    Py_DECREF(y_arr);
    Py_DECREF(N_arr);
    Py_DECREF(b_arr);
    Py_DECREF(z_arr);
    Py_DECREF(rest_arr);
    Py_DECREF(gamma_arr);
    Py_DECREF(f_arr);
    Py_DECREF(starts_arr);
    Py_DECREF(ends_arr);

    free(ab);
    ab=NULL;

    free(cont);
    cont=NULL;

    return Py_BuildValue("OOd", contin_out, abs_out, chi2);
}





