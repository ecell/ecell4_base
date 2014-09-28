#ifndef __ECELL4_RATELOW_HPP
#define __ECELL4_RATELOW_HPP

#include "types.hpp"
#include "exceptions.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/variant.hpp>
#include <vector>

namespace ecell4
{

class Ratelow
{
public:
    // The order of the species must be the same as
    //  reactants' container of ReactionRule object.
    typedef std::vector<Real> state_container_type;
public:
    Ratelow(){;}
    virtual ~Ratelow() {;}
    virtual bool is_available() const = 0;
    virtual Real operator()(state_container_type const &state_array, Real volume) = 0;
};

class RatelowCppCallback : public Ratelow
{
    /** Function object for calculate ratelow called by C++
     *  This class must not be exposed to Cython interface.
     */
public:
    typedef double (*Ratelow_Callback)(state_container_type const &, double);
public:
    RatelowCppCallback(Ratelow_Callback func) : func_(func) {;}
    RatelowCppCallback() : func_(0) {}
    virtual ~RatelowCppCallback(){;}

    virtual bool is_available() const
    {
        return this->func_ != 0;
    }
    virtual Real operator()(state_container_type const &state_array, Real volume)
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

class RatelowMassAction : public Ratelow
{
public:
    RatelowMassAction(Real k = 0.0, std::size_t num_reactant = 0) 
        : k_(k), num_reactant_(num_reactant) {}
    virtual ~RatelowMassAction(){;}
    virtual bool is_available() const
    {
        return true;    // always true
    }
    virtual Real operator()(state_container_type const &state_array, Real volume)
    {
        Real flux(this->k_ * volume);
        for(int i(0); i < this->num_reactant_; i++) {
            flux *= Real(state_array[i]) / volume;
        }
        return flux;
    }
    void set_k(Real k)
    {
        this->k_ = k;
    }
    Real get_k(Real k) const
    {
        return this->k_;
    }
private:
    Real k_;
    std::size_t num_reactant_;
};

} // ecell4

#endif  //__ECELL4_RATELOW_HPP
