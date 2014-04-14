#ifndef __ECELL4_UNIT_SPECIES_HPP
#define __ECELL4_UNIT_SPECIES_HPP

#include <string>

#include "config.h"

#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif


namespace ecell4
{

class UnitSpecies
{
public:

    typedef std::string serial_type;

public:

    UnitSpecies(const std::string& name)
        : name_(name)
    {
        ;
    }

    std::string name() const
    {
        return name_;
    }

    serial_type serial() const
    {
        return name_;
    }

    bool operator==(const UnitSpecies& rhs) const
    {
        return (name() == rhs.name());
    }

    bool operator<(const UnitSpecies& rhs) const
    {
        return (name() < rhs.name());
    }

    bool operator>(const UnitSpecies& rhs) const
    {
        return (name() > rhs.name());
    }

protected:

    std::string name_;
};

} // ecell4

#if defined(HAVE_TR1_FUNCTIONAL)
namespace std
{

namespace tr1
{
#elif defined(HAVE_STD_HASH)
namespace std
{
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
namespace boost
{
#endif

template<>
struct hash<ecell4::UnitSpecies>
{
    std::size_t operator()(const ecell4::UnitSpecies& val) const
    {
        return hash<ecell4::UnitSpecies::serial_type>()(val.serial());
    }
};

#if defined(HAVE_TR1_FUNCTIONAL)
} // tr1

} // std
#elif defined(HAVE_STD_HASH)
} // std
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
} // boost
#endif

#endif /* __ECELL4_SPECIES_HPP */

