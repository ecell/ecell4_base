#ifndef ECELL4_GILLESPIE_GILLESPIE_WORLD_HPP
#define ECELL4_GILLESPIE_GILLESPIE_WORLD_HPP

#include <stdexcept>
#include <sstream>
#include <map>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <string>

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/CompartmentSpace.hpp>
#ifdef WITH_HDF5
#include <ecell4/core/CompartmentSpaceHDF5Writer.hpp>
#endif
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Shape.hpp>
#include <ecell4/core/extras.hpp>


namespace ecell4
{

namespace gillespie
{

class GillespieWorld
    : public Space
{
public:

    GillespieWorld(const Real3& edge_lengths,
                   boost::shared_ptr<RandomNumberGenerator> rng)
        : cs_(new CompartmentSpaceVectorImpl(edge_lengths)), rng_(rng)
    {
        ;
    }

    GillespieWorld(const Real3& edge_lengths = Real3(1, 1, 1))
        : cs_(new CompartmentSpaceVectorImpl(edge_lengths))
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    GillespieWorld(const std::string filename)
        : cs_(new CompartmentSpaceVectorImpl(Real3(1, 1, 1)))
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        this->load(filename);
    }

    // SpaceTraits

    const Real t() const;
    void set_t(const Real& t);

    const Real3& edge_lengths() const
    {
        return cs_->edge_lengths();
    }

    void reset(const Real3& edge_lengths)
    {
        cs_->reset(edge_lengths);
    }

    // CompartmentSpaceTraits

    const Real volume() const
    {
        return cs_->volume();
    }

    Integer num_molecules(const Species& sp) const;
    Integer num_molecules_exact(const Species& sp) const;
    Real get_value(const Species& sp) const;
    Real get_value_exact(const Species& sp) const;
    void set_value(const Species& sp, const Real value);
    std::vector<Species> list_species() const;
    bool has_species(const Species& sp) const;

    // CompartmentSpace member functions

    void set_volume(const Real& volume)
    {
        (*cs_).set_volume(volume);
    }

    void add_molecules(const Species& sp, const Integer& num);
    void remove_molecules(const Species& sp, const Integer& num);

    // Optional members

    inline const boost::shared_ptr<RandomNumberGenerator>& rng()
    {
        return rng_;
    }

    void save(const std::string& filename) const
    {
#ifdef WITH_HDF5
        boost::scoped_ptr<H5::H5File>
            fout(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
        rng_->save(fout.get());
        boost::scoped_ptr<H5::Group>
            group(new H5::Group(fout->createGroup("CompartmentSpace")));
        cs_->save_hdf5(group.get());
        extras::save_version_information(fout.get(), std::string("ecell4-gillespie-") + std::string(ECELL4_VERSION));
#else
        throw NotSupported(
            "This method requires HDF5. The HDF5 support is turned off.");
#endif
    }

    void load(const std::string& filename)
    {
#ifdef WITH_HDF5
        boost::scoped_ptr<H5::H5File>
            fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));

        const std::string required = "ecell4-gillespie-4.1.0";
        try
        {
            const std::string version = extras::load_version_information(*fin);
            if (!extras::check_version_information(version, required))
            {
                std::stringstream ss;
                ss << "The version of the given file [" << version
                    << "] is too old. [" << required << "] or later is required.";
                throw NotSupported(ss.str());
            }
        }
        catch(H5::GroupIException not_found_error)
        {
            throw NotFound("No version information was found.");
        }

        rng_->load(*fin);
        const H5::Group group(fin->openGroup("CompartmentSpace"));
        cs_->load_hdf5(group);
#else
        throw NotSupported(
            "This method requires HDF5. The HDF5 support is turned off.");
#endif
    }

    void bind_to(boost::shared_ptr<Model> model)
    {
        if (boost::shared_ptr<Model> bound_model = lock_model())
        {
            if (bound_model.get() != model.get())
            {
                std::cerr << "Warning: Model already bound to GillespieWorld."
                    << std::endl;
            }
        }

        this->model_ = model;
    }

    boost::shared_ptr<Model> lock_model() const
    {
        return model_.lock();
    }

    void add_molecules(const Species& sp, const Integer& num, const boost::shared_ptr<Shape> shape)
    {
        add_molecules(sp, num);
    }

    std::pair<std::pair<ParticleID, Particle>, bool> new_particle(const Particle& p)
    {
        add_molecules(p.species(), 1);
        return std::make_pair(std::make_pair(ParticleID(), p), true);
    }

    std::pair<std::pair<ParticleID, Particle>, bool> new_particle(
        const Species& sp, const Real3& pos)
    {
        add_molecules(sp, 1);
        return std::make_pair(
            std::make_pair(ParticleID(), Particle(sp, pos, 0.0, 0.0)), true);
    }

    std::vector<std::pair<ParticleID, Particle> > list_particles() const;
    std::vector<std::pair<ParticleID, Particle> > list_particles_exact(const Species& sp) const;
    std::vector<std::pair<ParticleID, Particle> > list_particles(const Species& sp) const;

private:

    boost::scoped_ptr<CompartmentSpace> cs_;
    boost::shared_ptr<RandomNumberGenerator> rng_;

    boost::weak_ptr<Model> model_;
};

} // gillespie

} // ecell4

#endif /* ECELL4_GILLESPIE_GILLESPIE_WORLD_HPP */
