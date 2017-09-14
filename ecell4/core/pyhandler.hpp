#ifndef __ECELL4_PYHANDLER_HPP
#define __ECELL4_PYHANDLER_HPP


#include "exceptions.hpp"
#include "types.hpp"

namespace ecell4 {


struct PyObjectHandler {
    typedef void (*PyObjectHandler_1arg)(void*);
    typedef void *myPython_Object; 

    PyObjectHandler(PyObjectHandler_1arg inc_ref, PyObjectHandler_1arg dec_ref) :
        inc_ref_(inc_ref), dec_ref_(dec_ref)
    {;}
    ~PyObjectHandler();
    void inc_ref(myPython_Object obj)
    {
        if (this->inc_ref_ == 0)
        {
            throw IllegalState("Functions to Operate python reference counts are not registered");
        }
        this->inc_ref_( (void*) obj);
    }
    void dec_ref(myPython_Object obj)
    {
        if (this->dec_ref_ == 0)
        {
            throw IllegalState("Functions to Operate python reference counts are not registered");
        }
        this->dec_ref_( (void*) obj);
    }
    bool is_available(void)
    {
        if (this->inc_ref_ == 0) {
            return false;
        }
        if (this->dec_ref_ == 0) {
            return false;
        }
        // add criteria for checking if you have a new handler function
        return true;
    }
    
    PyObjectHandler_1arg inc_ref_, dec_ref_;

};


} // namespace ecell4
#endif /*  __ECELL4_PYHANDLER_HPP */
