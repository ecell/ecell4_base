
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


} // ode
} // ecell4
