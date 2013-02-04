#include <vector>
#include <string>
#include <sstream>
#include <iostream>

using namespace std;
#include "./GillespieWorld.hpp"

namespace ecell4
{

namespace gillespie
{

void GillespieWorld::set_t(Real const &t)
{
    this->cs_->set_t(t);
}
Real GillespieWorld::t(void)
{
    return this->cs_->t();
}

Integer GillespieWorld::num_species(void)
{
    return this->cs_->num_species();
}

bool GillespieWorld::has_species(Species const &sp)
{
    return this->cs_->has_species(sp);
}

Integer GillespieWorld::num_molecules(Species const &sp)
{
    return this->cs_->num_molecules(sp);
}

void GillespieWorld::add_species(Species const &sp)
{
    this->cs_->add_species(sp);
    return;
}

void GillespieWorld::remove_species(Species const &sp)
{
    this->cs_->remove_species(sp);
    return;
}

void GillespieWorld::add_molecules(Species const &sp, Integer const &num)
{
    this->cs_->add_molecules(sp, num);
    return;
}

void GillespieWorld::remove_molecules(Species const &sp, Integer const &num)
{
    this->cs_->remove_molecules(sp, num);
    return;
}

} //gillespie

} //ecell4
