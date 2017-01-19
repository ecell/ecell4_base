#ifndef __ECELL4_CALLBACK_HPP
#define __ECELL4_CALLBACK_HPP

#include <vector>
#include <ecell4/core/types.hpp>
#include <ecell4/core/exceptions.hpp>
#include "Space.hpp"
#include "CompartmentSpace.hpp"
#include "Species.hpp"
#include <boost/shared_ptr.hpp>

namespace ecell4
{


//template <typename T>
class CallbackWrapper
{
public:
    typedef void *pyfunc_type;
    typedef int (*stepladder_type)(pyfunc_type, std::vector<double>);

    CallbackWrapper(stepladder_type stepladder, pyfunc_type pyfunc);
    ~CallbackWrapper();
    bool is_available()
    {
        return (this->pyfunc_ != NULL && this->stepladder_ != NULL);
    }
    void call();
protected:
    pyfunc_type pyfunc_;
    stepladder_type stepladder_;
};

struct PythonNumberHooker
{
    //typedef std::vector<std::vector<Real> > data_container_type;
    //typedef std::vector<Species> species_container_type;

    // Callback Related TypeDefs
    typedef void *pyfunc_type;
    typedef bool (*stepladder_type)(pyfunc_type,boost::shared_ptr<CompartmentSpaceVectorImpl>);

    PythonNumberHooker(
            stepladder_type indirect, pyfunc_type pyfunc)
        :indirect_(indirect), pyfunc_(pyfunc) 
    {
        ;
    }
    ~PythonNumberHooker()
    {
        ;
    }
    bool is_available()
    {
        return (this->pyfunc_ != NULL && this->indirect_ != NULL);
    }

    // XXX What does this return value mean ???
    bool call(boost::shared_ptr<CompartmentSpaceVectorImpl> s)
    {
        if (this->is_available()) {
            return this->indirect_(this->pyfunc_, s);
        } else {
            throw IllegalState("Callback failed.");
        }
    }
    pyfunc_type pyfunc_;
    stepladder_type indirect_;
};

template <typename T1>
struct PythonHook_1arg
{
    typedef void *pyfunc_type;
    typedef bool (*stepladder_type)(pyfunc_type, T1);

    PythonHook_1arg(
            stepladder_type stepladder, pyfunc_type pyfunc)
        :stepladder_(stepladder), pyfunc_(pyfunc)
    {;}
    ~PythonHook_1arg()
    {
        ;
    }
    bool is_available()
    {
        return (this->pyfunc_ != NULL && this->stepladder_ != NULL);
    }

    // XXX What does this return value mean ???
    bool call(T1 obj)
    {
        if (this->is_available()) {
            return this->stepladder_(this->pyfunc_, obj);
        } else {
            throw IllegalState("Callback failed.");
        }
    }
    pyfunc_type pyfunc_;
    stepladder_type stepladder_;
};

struct PythonHook_Space
{
    typedef void *pyfunc_type;
    typedef bool (*stepladder_type_space)(pyfunc_type, boost::shared_ptr<Space> space);

    PythonHook_Space(
            stepladder_type_space stepladder, pyfunc_type pyfunc)
        :stepladder_(stepladder), pyfunc_(pyfunc)
    {;}
    ~PythonHook_Space()
    {
        ;
    }
    bool is_available()
    {
        return (this->pyfunc_ != NULL && this->stepladder_ != NULL);
    }

    // XXX What does this return value mean ???
    bool call(const boost::shared_ptr<Space>& space)
    {
        if (this->is_available()) {
            return this->stepladder_(this->pyfunc_, space);
        } else {
            throw IllegalState("Callback failed.");
        }
    }
    pyfunc_type pyfunc_;
    stepladder_type_space stepladder_;
};

struct PythonPositionHooker
{
};

struct PyObjectHandler
{
    typedef void* pylayer_object_pointer_type;
    typedef void (*pylayer_operation_function_type)(pylayer_object_pointer_type);
    PyObjectHandler(
            pylayer_operation_function_type inc_ref,
            pylayer_operation_function_type dec_ref
            ) :
        inc_ref_(inc_ref), dec_ref_(dec_ref)
    {
        ;
    }
    pylayer_operation_function_type inc_ref_, dec_ref_;
};

}
#endif /* __ECELL4_HPP */

