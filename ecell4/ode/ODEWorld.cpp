#include "ODEWorld.hpp"
#include <ecell4/core/extras.hpp>
#include <ecell4/core/exceptions.hpp>

namespace ecell4
{

namespace ode
{

void ODEWorld::bind_to(boost::shared_ptr<Model> model)
{
    if (generated_)
    {
        std::cerr << "Warning: NetworkModel is already bound to ODEWorld."
            << std::endl;
    }
    else if (model_.expired())
    {
        std::cerr << "Warning: ODENetworkModel is already bound to ODEWorld."
            << std::endl;
    }

    try
    {
        boost::shared_ptr<ODENetworkModel> tmp(new ODENetworkModel(model));
        generated_.swap(tmp);
        model_.reset();
    }
    catch (NotSupported e)
    {
        throw NotSupported(
            "Not supported yet. Either ODENetworkModel or NetworkModel must be given.");
    }
}

void ODEWorld::bind_to(boost::shared_ptr<ODENetworkModel> model)
{
    if (boost::shared_ptr<ODENetworkModel> bound_model = model_.lock())
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

void ODEWorld::save(const std::string& filename) const
{
#ifdef WITH_HDF5
    boost::scoped_ptr<H5::H5File>
        fout(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
    boost::scoped_ptr<H5::Group>
        group(new H5::Group(fout->createGroup("CompartmentSpace")));
    save_compartment_space<ODEWorldHDF5Traits<ODEWorld> >(*this, group.get());

    const uint32_t space_type = static_cast<uint32_t>(Space::ELSE);
    group->openAttribute("type").write(H5::PredType::STD_I32LE, &space_type);

    extras::save_version_information(fout.get(), "ecell4-ode-0.0-1");
#else
    throw NotSupported(
        "This method requires HDF5. The HDF5 support is turned off.");
#endif
}

void ODEWorld::load(const std::string& filename)
{
#ifdef WITH_HDF5
    boost::scoped_ptr<H5::H5File>
        fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));
    const H5::Group group(fin->openGroup("CompartmentSpace"));
    load_compartment_space<ODEWorldHDF5Traits<ODEWorld> >(group, this);
#else
    throw NotSupported(
        "This method requires HDF5. The HDF5 support is turned off.");
#endif
}

Real ODEWorld::evaluate(ODEReactionRule& rr)
{
    //XXX: rr cannot be const because ratelaw.deriv_func is not const.
    ODERatelaw::state_container_type::size_type cnt(0);
    const ODEReactionRule::reactant_container_type reactants(rr.reactants()); //XXX: ODEReactionRule should return the const reference of reactants, but not yet.
    ODERatelaw::state_container_type rvalues(reactants.size());
    for (ODEReactionRule::reactant_container_type::const_iterator
         i(reactants.begin()); i != reactants.end(); i++, cnt++)
    {
        rvalues[cnt] = get_value_exact(*i);
    }

    cnt = 0;
    const ODEReactionRule::product_container_type products(rr.products()); //XXX: ODEReactionRule should return the const reference of products, but not yet.
    ODERatelaw::state_container_type pvalues(products.size());
    for (ODEReactionRule::product_container_type::const_iterator
         i(products.begin()); i != products.end(); i++, cnt++)
    {
        pvalues[cnt] = get_value_exact(*i);
    }

    if (rr.has_ratelaw())
    {
        boost::shared_ptr<ODERatelaw> ratelaw(rr.get_ratelaw());
        if (!ratelaw->is_available())
        {
            throw IllegalState("The given rate law is not available.");
        }
        return ratelaw->deriv_func(rvalues, pvalues, volume(), t(), rr);
    }
    else
    {
        return ODERatelawMassAction(rr.k()).deriv_func(
            rvalues, pvalues, volume(), t(), rr);
    }
}

} // ode

} // ecell4
