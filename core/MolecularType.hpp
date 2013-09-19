#ifndef __ECELL4_MOLECULE_TYPE_HPP
#define __ECELL4_MOLECULE_TYPE_HPP

#include <vector>
#include "Species.hpp"
#include "Voxel.hpp"

#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

namespace ecell4
{

class MolecularType
{

public:

    typedef std::vector<Voxel*> voxel_container_type;

public:

    MolecularType(const std::string& name = "")
        : species_(name)
    {
    }
    void addVoxel(Voxel* voxel);
    void removeVoxel(const Voxel& voxel);
    const Species& species() const;
    const voxel_container_type& voxels() const;

protected:

    Species species_;
    voxel_container_type voxels_;

};

}

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
struct hash<ecell4::MolecularType>
{
    std::size_t operator()(const ecell4::MolecularType& val) const
    {
        return hash<ecell4::Species>()(val.species());
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

#endif
