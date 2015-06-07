#ifndef __ECELL4_ODE_RATELOW_HPP
#define __ECELL4_ODE_RATELOW_HPP

#include <ecell4/core/types.hpp>
#include <ecell4/core/exceptions.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/variant.hpp>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

// For Jacobi
#include <boost/numeric/ublas/matrix.hpp>


template <typename TInputIterator, typename T>
T cartesian_product(TInputIterator begin, TInputIterator end, T const init)
{
    return std::accumulate(begin, end, init, std::multiplies<T>());
}


namespace ecell4
{

namespace ode
{

class ODEReactionRule;

class ODERatelaw
{
public:

    // The order of the species must be the same as
    // reactants' container of ReactionRule object.
    //
    // state_container_type must be resized when called
    // jacobi_func and deriv_func.
    typedef std::vector<Real> state_container_type;
    typedef boost::numeric::ublas::matrix<double> matrix_type;

public:

    virtual bool is_available() const = 0;

    virtual Real deriv_func(
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, 
        Real const volume, Real const t,
        ODEReactionRule const &reaction) = 0;
};

class ODERatelawCppCallback
    : public ODERatelaw
{
public:
    /** Function object to calculate ratelaw called by C++
     *  This class must not be exposed to Cython interface.
     */

    // reactants_state, products_state, volume
    typedef double (*ODERatelaw_Callback)(
        state_container_type const &, state_container_type const &, double const);

public:

    ODERatelawCppCallback(ODERatelaw_Callback func)
        : func_(func), h_(1.0e-8)
    {
        ;
    }

    ODERatelawCppCallback()
        : func_(0), h_(1.0e-8)
    {
        ;
    }

    virtual ~ODERatelawCppCallback()
    {
        ;
    }

    virtual bool is_available() const
    {
        return this->func_ != 0;
    }

    virtual Real deriv_func(
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, 
        Real const volume, Real const t, ODEReactionRule const &rr);

    ODERatelaw_Callback get_callback() const
    {
        return this->func_;
    }

    ODERatelaw_Callback set_callback(ODERatelaw_Callback new_func)
    {
        if (new_func == 0)
        {
            throw std::invalid_argument("ODERatelaw Callback must not be 0");
        }
        ODERatelaw_Callback prev = get_callback();
        this->func_ = new_func;
        return prev;
    }

private:

    ODERatelaw_Callback func_;
    Real h_;
};

class ODERatelawCythonCallback
    : public ODERatelaw
{
public:
    /** Function object to calculate ratelaw called by Cython
     *  This class must not be used from C++ users' code.
     */

    // reactants_state, products_state, volume
    // typedef double (*ODERatelaw_Callback)(
    //     state_container_type const &, state_container_type const &, double const);
    typedef void* Python_CallbackFunctype;
    typedef double (*Stepladder_Functype)(
        Python_CallbackFunctype, 
        state_container_type, state_container_type, 
        Real volume, Real t, ODEReactionRule *rr);

public:

    ODERatelawCythonCallback(Stepladder_Functype indirect_func, void* pyfunc)
        : indirect_func_(indirect_func), python_func_(pyfunc), h_(1.0e-8)
    {
        ;
    }

    ODERatelawCythonCallback()
        : indirect_func_(0), python_func_(0), h_(1.0e-8) {;}

    virtual ~ODERatelawCythonCallback(){;}

    virtual bool is_available() const
    {
        return (this->indirect_func_ != 0 && this->python_func_ != 0);
    }

    virtual Real deriv_func(
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, 
        Real const volume, Real const t,
        ODEReactionRule const &rr);

    void set_callback_pyfunc(Python_CallbackFunctype new_func)
    {
        if (new_func == 0)
        {
            throw std::invalid_argument("ODERatelaw Callback must not be 0");
        }
        this->python_func_ = new_func;
    }

private:

    Stepladder_Functype indirect_func_;
    Python_CallbackFunctype python_func_;
    Real h_;
};

class ODERatelawMassAction
    : public ODERatelaw
{
public:

    ODERatelawMassAction(Real k = 0.0)
        : k_(k)
    {
        ;
    }

    virtual ~ODERatelawMassAction()
    {
        ;
    }

    virtual bool is_available() const
    {
        return true;    // always true
    }

    virtual Real deriv_func(
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, 
        Real const volume, Real const t, ODEReactionRule const &rr);

    void set_k(Real k)
    {
        this->k_ = k;
    }

    Real get_k() const
    {
        return this->k_;
    }

private:

    Real k_;
};

} // ode

} // ecell4

#endif  //__ECELL4_ODE_RATELOW_HPP
