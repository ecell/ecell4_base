#ifndef __ECELL4_PYHANDLER_HPP
#define __ECELL4_PYHANDLER_HPP

#include "exceptions.hpp"
#include "types.hpp"

namespace ecell4
{

struct PyObjectHandler
{
    typedef void (*single_argument_func_type)(void*);
    typedef void *pyobject_type;

    PyObjectHandler(single_argument_func_type inc_ref, single_argument_func_type dec_ref)
        : inc_ref_func(inc_ref), dec_ref_func(dec_ref)
    {
        ;
    }

    ~PyObjectHandler();

    void inc_ref(pyobject_type obj)
    {
        if (this->inc_ref_func == 0)
        {
            throw IllegalState("Functions to Operate python reference counts are not registered");
        }
        this->inc_ref_func((void*)obj);
    }

    void dec_ref(pyobject_type obj)
    {
        if (this->dec_ref_func == 0)
        {
            throw IllegalState("Functions to Operate python reference counts are not registered");
        }
        this->dec_ref_func((void*)obj);
    }

    bool is_available(void)
    {
        if (this->inc_ref_func == 0)
        {
            return false;
        }
        else if (this->dec_ref_func == 0)
        {
            return false;
        }

        // add criteria for checking if you have a new handler function
        return true;
    }

    single_argument_func_type inc_ref_func, dec_ref_func;
};

} // namespace ecell4

#endif /*  __ECELL4_PYHANDLER_HPP */
