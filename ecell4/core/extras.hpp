#ifndef ECELL4_EXTRAS_HPP
#define ECELL4_EXTRAS_HPP

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include <ecell4/core/config.h>

#ifdef WITH_HDF5
#include <hdf5.h>
#include <H5Cpp.h>
#endif

#include "types.hpp"
#include "Real3.hpp"
#include "Species.hpp"
#include "Particle.hpp"
#include "AABB.hpp"
#include "Model.hpp"


namespace ecell4
{

namespace extras
{

template<typename Tworld_, typename Trng_>
void throw_in_particles(
    Tworld_& world, const Species& sp, const Integer& N,
    const boost::shared_ptr<Shape> shape,
    boost::shared_ptr<Trng_>& rng)
{
    typedef typename Tworld_::molecule_info_type molecule_info_type;
    boost::shared_ptr<RandomNumberGenerator>
        myrng(static_cast<boost::shared_ptr<RandomNumberGenerator> >(rng));

    if (N < 0)
    {
        throw std::invalid_argument("the number of particles must be positive.");
    }

    // const Real3 edge_lengths(world.edge_lengths());
    const molecule_info_type info(world.get_molecule_info(sp));

    for (int i(0); i < N; ++i)
    {
        while (true)
        {
            // const Real3 pos(
            //     rng.uniform(0.0, edge_lengths[0]),
            //     rng.uniform(0.0, edge_lengths[1]),
            //     rng.uniform(0.0, edge_lengths[2]));
            // if (world.list_particles_within_radius(pos, info.radius).size()
            //     == 0)
            // {
            //     world.new_particle(Particle(sp, pos, info.radius, info.D));
            //     break;
            // }
            const Real3 pos(shape->draw_position(myrng));
            if (world.new_particle(Particle(sp, pos, info.radius, info.D)).second)
            {
                break;
            }
        }
    }
}

template<typename Tworld_, typename Trng_>
void throw_in_particles(
    Tworld_& world, const Species& sp, const Integer& N, boost::shared_ptr<Trng_>& rng)
{
    boost::shared_ptr<Shape> shape(new AABB(Real3(0, 0, 0), world.edge_lengths()));
    throw_in_particles(world, sp, N, shape, rng);
}

template<typename Tfactory_>
typename Tfactory_::world_type* generate_world_from_model(
    const Tfactory_& f, const boost::shared_ptr<Model>& m)
{
    typename Tfactory_::world_type* w = f.create_world();
    w->bind_to(m);
    return w;
}

struct VersionInformation
{
    std::string header;
    int majorno, minorno, patchno;

    VersionInformation(
        const std::string& header, const int majorno, const int minorno, const int patchno)
        : header(header), majorno(majorno), minorno(minorno), patchno(patchno)
    {
        ;
    }
};

VersionInformation parse_version_information(const std::string& version);
bool check_version_information(const std::string& version, const std::string& required);

#ifdef WITH_HDF5
void save_version_information(H5::CommonFG* root, const std::string& version);
std::string load_version_information(const H5::CommonFG& root);
#endif
std::string load_version_information(const std::string& filename);

} // extras

} // ecell4

#endif // ECELL4_EXTRAS_HPP
