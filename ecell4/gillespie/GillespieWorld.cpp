#include <vector>
#include <string>
#include <sstream>
#include <iostream>

#include "GillespieWorld.hpp"


// using namespace std;

namespace ecell4
{

namespace gillespie
{

void GillespieWorld::set_t(const Real& t)
{
    this->cs_->set_t(t);
}

const Real& GillespieWorld::t(void) const
{
    return this->cs_->t();
}

std::vector<Species> GillespieWorld::list_species() const
{
    return this->cs_->list_species();
}

Integer GillespieWorld::num_molecules(const Species& sp) const
{
    return this->cs_->num_molecules(sp);
}

Integer GillespieWorld::num_molecules_exact(const Species& sp) const
{
    return this->cs_->num_molecules_exact(sp);
}

void GillespieWorld::add_molecules(const Species& sp, const Integer& num)
{
    this->cs_->add_molecules(sp, num);
}

void GillespieWorld::remove_molecules(const Species& sp, const Integer& num)
{
    this->cs_->remove_molecules(sp, num);
}

bool GillespieWorld::has_species(const Species& sp) const
{
    return this->cs_->has_species(sp);
}

} // gillespie

} // ecell4
