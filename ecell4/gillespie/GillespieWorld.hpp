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
#include <ecell4/core/CompartmentSpaceHDF5Writer.hpp>
#include <ecell4/core/NetworkModel.hpp>


namespace ecell4
{

namespace gillespie
{

class GillespieWorld
{
public:

    GillespieWorld(const Position3& edge_lengths,
                   boost::shared_ptr<RandomNumberGenerator> rng)
        : cs_(new CompartmentSpaceVectorImpl(edge_lengths)), rng_(rng)
    {
        ;
    }

    GillespieWorld(const Position3& edge_lengths)
        : cs_(new CompartmentSpaceVectorImpl(edge_lengths))
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    // SpaceTraits

    const Real& t(void) const;
    void set_t(const Real& t);

    const Position3& edge_lengths() const
    {
        return cs_->edge_lengths();
    }

    void set_edge_lengths(const Position3& edge_lengths)
    {
        cs_->set_edge_lengths(edge_lengths);
    }

    // CompartmentSpaceTraits

    const Real volume() const
    {
        return cs_->volume();
    }

    Integer num_molecules(const Species& sp) const;
    std::vector<Species> list_species() const;

    // CompartmentSpace member functions

    void set_volume(const Real& volume)
    {
        (*cs_).set_volume(volume);
    }

    void add_molecules(const Species& sp, const Integer& num);
    void remove_molecules(const Species& sp, const Integer& num);

    // Optional members

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return rng_;
    }

    void save(const std::string& filename) const
    {
        boost::scoped_ptr<H5::H5File>
            fout(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
        rng_->save(fout.get());
        boost::scoped_ptr<H5::Group>
            group(new H5::Group(fout->createGroup("CompartmentSpace")));
        cs_->save(group.get());
    }

    void load(const std::string& filename)
    {
        boost::scoped_ptr<H5::H5File>
            fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));
        rng_->load(*fin);
        const H5::Group group(fin->openGroup("CompartmentSpace"));
        cs_->load(group);
    }

    void bind_to(boost::shared_ptr<NetworkModel> model)
    {
        if (boost::shared_ptr<NetworkModel> bound_model = lock_model())
        {
            if (bound_model.get() != model.get())
            {
                std::cerr << "Warning: Model already bound to GillespieWorld."
                    << std::endl;
            }
        }
        this->model_ = model;
    }

    boost::shared_ptr<NetworkModel> lock_model() const
    {
        return model_.lock();
    }

private:

    boost::scoped_ptr<CompartmentSpace> cs_;
    boost::shared_ptr<RandomNumberGenerator> rng_;

    boost::weak_ptr<NetworkModel> model_;
};

} // gillespie

} // ecell4

#endif /* __ECELL4_GILLESPIE_GILLESPIE_WORLD_HPP */
