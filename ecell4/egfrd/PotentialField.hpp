#ifndef POTENTIAL_FIELD_HPP
#define POTENTIAL_FIELD_HPP

#include <boost/shared_ptr.hpp>

#include <ecell4/core/get_mapper_mf.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/Shape.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

namespace ecell4
{

template <typename Tcontainer>
class PotentialField
{
public:

    typedef Tcontainer container_type;
    typedef RandomNumberGenerator rng_type;
    typedef std::pair<ParticleID, Particle> pid_particle_pair;

public:

    virtual bool try_move(rng_type& rng, const pid_particle_pair& p, const Real3& newpos, const container_type& space)
    {
        return true;
    }
};

template <typename Tcontainer>
class ShapedPotentialField
    : public PotentialField<Tcontainer>
{
public:

    typedef PotentialField<Tcontainer> base_type;
    typedef typename base_type::container_type container_type;
    typedef typename base_type::pid_particle_pair pid_particle_pair;
    typedef typename base_type::rng_type rng_type;

public:

    ShapedPotentialField(const boost::shared_ptr<Shape>& shape)
        : shape_(shape)
    {
        ;
    }

    bool try_move(rng_type& rng, const pid_particle_pair& p, const Real3& newpos, const container_type& space)
    {
        return (shape_->is_inside(newpos) <= 0.0);
    }

protected:

    boost::shared_ptr<Shape> shape_;
};

template <typename Tcontainer>
class LeashPotentialField
    : public PotentialField<Tcontainer>
{
public:

    typedef PotentialField<Tcontainer> base_type;
    typedef typename base_type::container_type container_type;
    typedef typename base_type::pid_particle_pair pid_particle_pair;
    typedef typename base_type::rng_type rng_type;

    typedef typename ecell4::utils::get_mapper_mf<ParticleID, Real3>::type particle_id_position_map_type;

public:

    LeashPotentialField(const Real radius)
        : radius_(radius)
    {
        ;
    }

    bool try_move(rng_type& rng, const pid_particle_pair& p, const Real3& newpos, const container_type& space)
    {
        particle_id_position_map_type::const_iterator it = centers_.find(p.first);
        if (it == centers_.end())
        {
            centers_.insert(particle_id_position_map_type::value_type(p.first, p.second.position()));
            return (space.distance(p.second.position(), newpos) <= radius_);
        }
        return (space.distance((*it).second, newpos) <= radius_);
    }

protected:

    Real radius_;
    particle_id_position_map_type centers_;
};



} // ecell4

#endif /* POTENTIAL_FIELD_HPP */
