#include "PythonHooker.hpp"

namespace ecell4
{

void FixedIntervalPythonHooker::initialize(const boost::shared_ptr<WorldInterface>& world, const boost::shared_ptr<Model>& model)
{
    base_type::initialize(world, model);
}

bool FixedIntervalPythonHooker::fire(const Simulator* sim, const boost::shared_ptr<WorldInterface>& world)
{
    const bool ret1 = this->stepladder_(this->pyfunc_, world, sim->check_reaction());
    const bool ret2 = base_type::fire(sim, world);
    return (ret1 & ret2);
}

void FixedIntervalPythonHooker::reset()
{
    base_type::reset();
}

} // ecell4
