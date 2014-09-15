#ifndef __ECELL4_RATELOW_HPP
#define __ECELL4_RATELOW_HPP

#include "types.hpp"
#include "exceptions.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/variant.hpp>

namespace ecell4
{

class Ratelow
{
public:
    Ratelow(){;}
    virtual ~Ratelow() {;}
    virtual bool is_available() const = 0;
    virtual Real operator()(Real *state_array, Real value) = 0;
};

class RatelowCppImpl : public Ratelow
{
    /** Function object for calculate ratelow called by C++
     *  This class must not be exposed to Cython interface.
     */
public:
    typedef double (*Ratelow_Callback)(double *, double);
public:
    RatelowCppImpl(Ratelow_Callback func) : func_(func) {;}
    RatelowCppImpl() : func_(0) {}
    virtual ~RatelowCppImpl(){;}

    virtual bool is_available() const
    {
        return this->func_ != 0;
    }
    virtual Real operator()(Real *state_array, Real volume)
    {
        if (!is_available() )
        {
            throw IllegalState("Callback Function has not been registerd");
        }
        return this->func_(state_array, volume);
    }
    Ratelow_Callback get_callback() const
    {
        return this->func_;
    }
    Ratelow_Callback set_callback(Ratelow_Callback new_func)
    {
        if (new_func == 0)
        {
            throw std::invalid_argument("Ratelow Callback must not be 0");
        }
        Ratelow_Callback prev = get_callback();
        this->func_ = new_func;
        return prev;
    }
private:
    Ratelow_Callback func_;
};

} // ecell4

#endif  //__ECELL4_RATELOW_HPP
