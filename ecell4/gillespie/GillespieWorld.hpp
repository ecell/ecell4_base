#ifndef __ECELL4_GILLESPIE_GILLESPIE_WORLD_HPP
#define __ECELL4_GILLESPIE_GILLESPIE_WORLD_HPP

#include <stdexcept>
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

    const Real& t(void) const;
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
        cs_->save(group.get());
#else
        throw NotSupported("not supported yet.");
#endif
    }

    void load(const std::string& filename)
    {
#ifdef WITH_HDF5
        boost::scoped_ptr<H5::H5File>
            fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));
        rng_->load(*fin);
        const H5::Group group(fin->openGroup("CompartmentSpace"));
        cs_->load(group);
#else
        throw NotSupported("not supported yes.");
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
                extras::set_parameters(*model, *this);
            }
        }
        else
        {
            extras::set_parameters(*model, *this);
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

#endif /* __ECELL4_GILLESPIE_GILLESPIE_WORLD_HPP */
