#ifndef __ECELL4_RATELOW_HPP
#define __ECELL4_RATELOW_HPP

#include "types.hpp"
#include "exceptions.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/variant.hpp>
#include <vector>
#include <algorithm>

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
    virtual Real deriv_func(state_container_type const &state_array, Real volume);
    virtual void jacobi_func(state_container_type const &state_array, state_container_type &jacobian, Real volume) = 0;
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
    virtual Real deriv_func(state_container_type const &state_array, Real volume)
    {
        // This is the same as deriv_func
        // Should they be united to this function?
        return (*this)(state_array, volume);
    }
    virtual void jacobi_func(
            state_container_type const &state_array, state_container_type &jacobian, Real volume)
    {
        Real h( 1.0e-8 ); //XXX  Temporary using 1.0e-8. Should be fixed
        std::fill(jacobian.begin(), jacobian.end(), Real(0.0));
        Real flux( this->deriv_func(state_array, volume) );
        for(int i(0); i < jacobian.size(); i++) 
        {
            //XXX For now, we are using FORWARD difference method.
            state_container_type h_shift(state_array);
            h_shift[i] += h;
            jacobian[i] = ((this->deriv_func(h_shift, volume)) - flux) / h;
        }
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
    virtual Real deriv_func(state_container_type const &state_array, Real volume)
    {
        Real flux( (*this)(state_array, volume) );
        return flux;
    }
    virtual void jacobi_func(
            state_container_type const &state_array, state_container_type &jacobian, 
            Real volume)
    {
        // The size of the argument 'jacobian' must be resized to 
        //  the number of reactants.
        std::fill(jacobian.begin(), jacobian.end(), Real(0.0));
        Real flux( this->deriv_func(state_array, volume) );
        if (flux == Real(0.0))
        {
            return;
        }
        for(int i(0); i < state_array.size(); i++) {
            Real partial(flux / state_array[i]);
            jacobian[i] = partial;
        }
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
