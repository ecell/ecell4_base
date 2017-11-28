
#include "ODESimulator_New.hpp"

#include <boost/numeric/odeint.hpp>
#include <algorithm>

namespace odeint = boost::numeric::odeint;

namespace ecell4
{

namespace ode
{


//void ODEWorld_New::bind_to(boost::shared_ptr<Model> model)
//{
//    if (generated_)
//    {
//        std::cerr << "Warning: NetworkModel is already bound to ODEWorld."
//            << std::endl;
//    }
//    else if (model_.expired())
//    {
//        std::cerr << "Warning: ODENetworkModel is already bound to ODEWorld."
//            << std::endl;
//    }
//
//    try
//    {
//        boost::shared_ptr<NetworkModel> tmp(new NetworkModel(model));
//        generated_.swap(tmp);
//        model_.reset();
//    }
//    catch (NotSupported e)
//    {
//        throw NotSupported(
//            "Not supported yet. Either ODENetworkModel or NetworkModel must be given.");
//    }
//}

void ODEWorld_New::save(const std::string& filename) const
{
    throw NotImplemented("ODEWorld_new::save not implemented");
}

void ODEWorld_New::load(const std::string& filename)
{
    throw NotImplemented("ODEWorld_new::load not implemented");
}

void ODEWorld_New::bind_to(boost::shared_ptr<NetworkModel> model)
{
    if (boost::shared_ptr<NetworkModel> bound_model = model_.lock())
    {
        if (bound_model.get() != model.get())
        {
            std::cerr << "Warning: ODENetworkModel is already bound to ODEWorld."
                << std::endl;
        }
    }
    else if (generated_)
    {
        std::cerr << "Warning: NetworkModel is already bound to ODEWorld."
            << std::endl;
    }

    this->model_ = model;
    generated_.reset();
}

bool ODESimulator_New::step(const Real &upto)
{
    if (upto <= t()) {
        return false;
    }
    const Real dt(std::min(dt_,upto - t()));
    const Real ntime(std::min(upto, t() + dt_));
    const std::vector<Species> species(world_->list_species());
    state_type x(species.size());
    state_type::size_type i(0);
    for (NetworkModel::species_container_type::const_iterator it(species.begin());
            it != species.end(); it++) 
    {
        x[i] = static_cast<double>(world_->get_value_exact(*it));
        i++;
    }
    //std::pair<deriv_func, jacobi_func> system(generate_system());
    //StateAndTimeBackInserter::state_container_type x_vec;
    //StateAndTimeBackInserter::time_container_type times;
    return true;
}

} // ode
} // ecell4
