#ifndef __ECELL4_ODE_RATELOW_HPP
#define __ECELL4_ODE_RATELOW_HPP

#include <ecell4/core/types.hpp>
#include <ecell4/core/exceptions.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/variant.hpp>
#include <vector>
#include <algorithm>

// For Jacobi
#include <boost/numeric/ublas/matrix.hpp>

namespace ecell4
{

namespace ode
{

class Ratelaw
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

    virtual Real operator()(
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, Real const volume) = 0;

    virtual Real deriv_func(
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, Real const volume) = 0;

    virtual void jacobi_func(
        matrix_type &jacobian,
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array,
        Real const volume) = 0;
};

class RatelawCppCallback
    : public Ratelaw
{
public:
    /** Function object to calculate ratelaw called by C++
     *  This class must not be exposed to Cython interface.
     */

    // reactants_state, products_state, volume
    typedef double (*Ratelaw_Callback)(
        state_container_type const &, state_container_type const &, double const);

public:

    RatelawCppCallback(Ratelaw_Callback func)
        : func_(func), h_(1.0e-8)
    {
        ;
    }

    RatelawCppCallback()
        : func_(0), h_(1.0e-8)
    {
        ;
    }

    virtual ~RatelawCppCallback()
    {
        ;
    }

    virtual bool is_available() const
    {
        return this->func_ != 0;
    }

    virtual Real operator()(
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, Real const volume)
    {
        return this->deriv_func(reactants_state_array, products_state_array, volume);
    }

    virtual Real deriv_func(
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, Real volume)
    {
        if (!is_available())
        {
            throw IllegalState("Callback Function has not been registerd");
        }
        return this->func_(reactants_state_array, products_state_array, volume);
    }

    virtual void jacobi_func(
        matrix_type &jacobian,
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array,
        Real const volume)
    {
        Real h(this->h_); //XXX: 1.0e-8. Should be fixed

        std::fill(jacobian.data().begin(), jacobian.data().end(), Real(0.0));
        Real flux(this->deriv_func(reactants_state_array, products_state_array, volume));
        double num_reactants(reactants_state_array.size());
        double num_products(products_state_array.size());

        // Differentiates by Reactants.
        for (int i(0); i < num_reactants; i++)
        {
            //XXX: For now, we are using FORWARD difference method.
            state_container_type h_shift(reactants_state_array);
            h_shift[i] += h;
            double deriv = (
                this->deriv_func(h_shift, products_state_array, volume) - flux) / h;
            for (matrix_type::size_type j(0); j < jacobian.size1(); j++)
            {
                if (j < num_reactants)
                {
                    jacobian(j, i) -= deriv;
                }
                else
                {
                    jacobian(j, i) += deriv;
                }
            }
        }

        // Differentiates by Products.
        for (int i(0); i < num_products; i++)
        {
            state_container_type h_shift(products_state_array);
            h_shift[i] += h;
            double deriv = (
                this->deriv_func(reactants_state_array, h_shift, volume) - flux) / h;
            for (matrix_type::size_type j(0); j < jacobian.size1(); j++)
            {
                if (j < num_reactants)
                {
                    jacobian(j, i + num_reactants) -= deriv;
                }
                else
                {
                    jacobian(j, i + num_reactants) += deriv;
                }
            }
        }
    }

    Ratelaw_Callback get_callback() const
    {
        return this->func_;
    }

    Ratelaw_Callback set_callback(Ratelaw_Callback new_func)
    {
        if (new_func == 0)
        {
            throw std::invalid_argument("Ratelaw Callback must not be 0");
        }
        Ratelaw_Callback prev = get_callback();
        this->func_ = new_func;
        return prev;
    }

private:

    Ratelaw_Callback func_;
    Real h_;
};

class RatelawCythonCallback
    : public Ratelaw
{
public:
    /** Function object to calculate ratelaw called by Cython
     *  This class must not be used from C++ users' code.
     */

    // reactants_state, products_state, volume
    // typedef double (*Ratelaw_Callback)(
    //     state_container_type const &, state_container_type const &, double const);
    typedef void* Python_Functype;
    typedef double (*Indirect_Functype)(
        Python_Functype, state_container_type, state_container_type, Real);

public:

    RatelawCythonCallback(Indirect_Functype indirect_func, void* pyfunc)
        : indirect_func_(indirect_func), python_func_(pyfunc), h_(1.0e-8)
    {
        ;
    }

    RatelawCythonCallback()
        : indirect_func_(0), python_func_(0), h_(1.0e-8) {;}

    virtual ~RatelawCythonCallback(){;}

    virtual bool is_available() const
    {
        return (this->indirect_func_ != 0 && this->python_func_ != 0);
    }

    virtual Real operator()(
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, Real const volume)
    {
        return this->deriv_func(reactants_state_array, products_state_array, volume);
    }

    virtual Real deriv_func(
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, Real volume)
    {
        if (!is_available())
        {
            throw IllegalState("Callback Function has not been registerd");
        }
        return this->indirect_func_(
            this->python_func_, reactants_state_array, products_state_array, volume);
    }

    virtual void jacobi_func(
        matrix_type &jacobian,
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array,
        Real const volume)
    {
        Real h(this->h_); //XXX: 1.0e-8. Should be fixed
        std::fill(jacobian.data().begin(), jacobian.data().end(), Real(0.0));
        Real flux(this->deriv_func(reactants_state_array, products_state_array, volume));
        double num_reactants(reactants_state_array.size());
        double num_products(products_state_array.size());

        // Differentiates by Reactants.
        for (int i(0); i < num_reactants; i++)
        {
            //XXX: For now, we are using FORWARD difference method.
            state_container_type h_shift(reactants_state_array);
            h_shift[i] += h;
            double deriv = (
                this->deriv_func(h_shift, products_state_array, volume) - flux) / h;
            for (matrix_type::size_type j(0); j < jacobian.size1(); j++)
            {
                if (j < num_reactants)
                {
                    jacobian(j, i) -= deriv;
                }
                else
                {
                    jacobian(j, i) += deriv;
                }
            }
        }

        // Differentiates by Products.
        for (int i(0); i < num_products; i++)
        {
            state_container_type h_shift(products_state_array);
            h_shift[i] += h;
            double deriv = (
                this->deriv_func(reactants_state_array, h_shift, volume) - flux) / h;
            for (matrix_type::size_type j(0); j < jacobian.size1(); j++)
            {
                if (j < num_reactants)
                {
                    jacobian(j, i + num_reactants) -= deriv;
                }
                else
                {
                    jacobian(j, i + num_reactants) += deriv;
                }
            }
        }

        return; //XXX:
    }

    void set_callback_pyfunc(Python_Functype new_func)
    {
        if (new_func == 0)
        {
            throw std::invalid_argument("Ratelaw Callback must not be 0");
        }
        this->python_func_ = new_func;
    }

private:

    Python_Functype python_func_;
    Indirect_Functype indirect_func_;
    Real h_;
};

class RatelawMassAction
    : public Ratelaw
{
public:

    RatelawMassAction(Real k = 0.0)
        : k_(k)
    {
        ;
    }

    virtual ~RatelawMassAction()
    {
        ;
    }

    virtual bool is_available() const
    {
        return true;    // always true
    }

    virtual Real operator()(
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, Real const volume)
    {
        // Forward to deriv_func()
        Real flux(this->deriv_func(reactants_state_array, products_state_array, volume));
        return flux;
    }

    virtual Real deriv_func(
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, Real const volume)
    {
        // The 1st argument 'state_array' must be resized when calling.
        Real flux(this->k_ * volume);
        for(state_container_type::const_iterator it(reactants_state_array.begin());
            it != reactants_state_array.end(); it++)
        {
            flux *= (*it) / volume;
        }
        return flux;
    }

    virtual void jacobi_func(
        matrix_type &jacobian,
        state_container_type const &reactants_state_array,
        state_container_type const &products_state_array, Real const volume)
    {
        // The size of the argument 'state_array' must be resized
        // to the number of reactants.
        // The size of the argument 'jacobian' must be resized
        // to the number of (reactants + products)
        std::fill(jacobian.data().begin(), jacobian.data().end(), Real(0.0));
        Real flux(this->deriv_func(reactants_state_array, products_state_array, volume));
        if (flux == Real(0.0))
        {
            return;
        }

        matrix_type::size_type num_reactants(
            static_cast<matrix_type::size_type>(reactants_state_array.size()));
        for(matrix_type::size_type i(0); i < num_reactants; i++)
        {
            Real partial(flux / reactants_state_array[i]);
            for(matrix_type::size_type j(0); j < jacobian.size1(); j++)
            {
                if (j < num_reactants)
                {
                    jacobian(j, i) -= partial;
                }
                else
                {
                    jacobian(j, i) += partial;
                }
            }
        }
    }

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
