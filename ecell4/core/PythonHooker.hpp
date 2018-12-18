#ifndef ECELL4_PYTHON_HOOKER_HPP
#define ECELL4_PYTHON_HOOKER_HPP

#include "Python.h"

#include "types.hpp"
#include "functions.hpp"
#include "Space.hpp"
#include "Species.hpp"
#include "Real3.hpp"
#include "Simulator.hpp"
#include "WorldInterface.hpp"
#include "observers.hpp"

#include <fstream>
#include <boost/format.hpp>
#include <time.h>


namespace ecell4
{

class FixedIntervalPythonHooker
    : public FixedIntervalObserver
{
public:

    typedef FixedIntervalObserver base_type;

    typedef PyObject* pyfunc_type;
    typedef bool (*stepladder_func_type)(pyfunc_type, const boost::shared_ptr<WorldInterface>& world, bool check_reaction);

public:

    FixedIntervalPythonHooker(
        const Real &dt, stepladder_func_type stepladder, pyfunc_type pyfunc)
        : base_type(dt), stepladder_(stepladder), pyfunc_(pyfunc)
    {
        Py_INCREF(this->pyfunc_);
    }

    ~FixedIntervalPythonHooker()
    {
        Py_DECREF(this->pyfunc_);
    }

    virtual void initialize(const boost::shared_ptr<WorldInterface>& world, const boost::shared_ptr<Model>& model);
    virtual bool fire(const Simulator* sim, const boost::shared_ptr<WorldInterface>& world);
    virtual void reset();

protected:

    stepladder_func_type stepladder_;
    pyfunc_type pyfunc_;
};

} // ecell4

#endif /* ECELL4_PYTHON_HOOKER_HPP */
