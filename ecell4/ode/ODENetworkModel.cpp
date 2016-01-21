
#include "ODENetworkModel.hpp"

#include <algorithm>
#include <iterator>

namespace ecell4
{
namespace ode
{

ODENetworkModel::ODENetworkModel()
{
    ;
}

ODENetworkModel::ODENetworkModel(const boost::shared_ptr<ecell4::Model> model)
    :expanded_(model)
{
    this->convert_from_networkmodel(model);
}

void ODENetworkModel::update_model(void)
{
    if (this->has_network_model() == false)
    {
        throw IllegalState("ecell4::Model object has not been registered");
    }
    this->convert_from_networkmodel(this->get_networkmodel());
}

bool ODENetworkModel::convert_from_networkmodel(const boost::shared_ptr<ecell4::Model> model)
{
    if (!model->is_static())
    {
        throw NotSupported("Not supported yet. NetworkModel must be given.");
    }
    this->species_attributes_.clear();
    this->species_attributes_.reserve(model->species_attributes().size());
    std::copy(
            model->species_attributes().begin(), model->species_attributes().end(), 
            std::back_inserter(this->species_attributes_) );

    this->ode_reaction_rules_.clear();
    const network_model_type::reaction_rule_container_type source(model->reaction_rules());
    for(network_model_type::reaction_rule_container_type::const_iterator it(source.begin());
            it != source.end(); it++)
    {
        ode_reaction_rule_type ode_rr(*it);
        this->ode_reaction_rules_.push_back(ode_rr);
    }
    return true;
}

ODENetworkModel::~ODENetworkModel()
{
    ;
}

}   // ode

}   // ecell4
