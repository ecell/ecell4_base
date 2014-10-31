#ifndef PEER_COMPAT_H
#define PEER_COMPAT_H

#if PY_VERSION_HEX < 0x02050000

typedef inquiry lenfunc;
typedef intargfunc ssizeargfunc;
typedef intintargfunc ssizessizeargfunc;
typedef intobjargproc ssizeobjargproc;
typedef intintobjargproc ssizessizeobjargproc;
typedef int Py_ssize_t;

static Py_ssize_t PyNumber_AsSsize_t(PyObject *item, PyObject* err)
{
    Py_ssize_t retval = PyInt_AsLong(item);
    if (retval == -1 || PyErr_Occurred()) {
        if (err) {
            PyErr_Format(err,"cannot convert '%.200s' to Py_ssize_t",
                    item->ob_type->tp_name);
        }
    }
    return retval;
}

static Py_ssize_t _PyObject_LengthHint(PyObject *obj)
{
    return -1;
}

#endif /* PY_VERSION_HEX < 0x02050000 */

#if PY_VERSION_HEX < 0x02060000
#define compat_PyObject_LengthHint(a) _PyObject_LengthHint(a)
#else
#define compat_PyObject_LengthHint(a) _PyObject_LengthHint(a, -1)
#endif

#endif /* PEER_COMPAT_H */
