#ifndef __ECELL4_ODE_RATELOW_HPP
#define __ECELL4_ODE_RATELOW_HPP

#include "exceptions.hpp"

namespace ecell4
{


class PythonAPIBind
{
public:
    typedef void (*PyRefFuncType)(void*);

public:
    PythonAPIBind();
    PythonAPIBind(PyRefFuncType inc_ref, PyRefFuncType dec_ref);
    void inc_ref(void* object)
    {
        if (this->inc_ref_ == NULL)
        {
            throw IllegalState("Python API is not binded");
        }
        this->inc_ref_(object);
    }
    void dec_ref(void* object)
    {
        if (this->dec_ref_ == NULL)
        {
            throw IllegalState("Python API is not binded");
        }
        this->dec_ref_(object);
    }

    bool is_initialized()
    {
        if ( (this->inc_ref_ != NULL) && (this->dec_ref_ != NULL) )
        {
            return true;
        } else {
            return false;
        }
    }

protected:
    PyRefFuncType inc_ref_, dec_ref_;
};



} // ecell4
#endif  //__ECELL4_ODE_RATELOW_HPP




