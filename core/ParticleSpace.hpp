#ifndef __PARTICLE_SPACE_HPP
#define __PARTICLE_SPACE_HPP

#include <cmath>
// #include <gsl/gsl_pow_int.h>

#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "functions.hpp"
#include "exceptions.hpp"
#include "Position3.hpp"
#include "Particle.hpp"
#include "Species.hpp"
#include "Space.hpp"

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
    using std::cout;
    using std::endl;
#endif  // H5_NO_STD
#endif

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

const H5std_string FILE_NAME( "hoge.h5" );
const H5std_string DATASET_NAME( "TimePoint" );
const H5std_string MEMBER1( "particle_id" );
const H5std_string MEMBER2( "positions" );

namespace ecell4
{

Real pow_2(Real const& a);

class ParticleSpace
    : public Space
{
public:

    typedef std::vector<std::pair<ParticleID, Particle> >
    particle_container_type;

    virtual Position3 const& edge_lengths() const
    {
        throw NotImplemented("edge_lengths() not implemented");
    }

    virtual Integer num_particles() const
    {
        throw NotImplemented("num_particles() not implemented");
    }

    virtual Integer num_particles(Species const& species) const
    {
        throw NotImplemented("num_particles() not implemented");
    }

    virtual std::vector<std::pair<ParticleID, Particle> >
    list_particles(Species const& species) const
    {
        throw NotImplemented("list_particles() not implemented");
    }

    virtual bool has_particle(ParticleID const& pid) const = 0;
    virtual bool update_particle(ParticleID const& pid, Particle const& p) = 0;
    virtual void remove_particle(ParticleID const& pid) = 0;

    virtual particle_container_type const& particles() const = 0;

    virtual std::pair<ParticleID, Particle>
    get_particle(ParticleID const& pid) const = 0;
    virtual std::vector<std::pair<ParticleID, Particle> >
    list_particles() const = 0;

    virtual void insert_space(particle_container_type particle_container) = 0;

    Position3 periodic_transpose(
        Position3 const& pos1, Position3 const& pos2) const
    {
        Position3 retval(pos1);
        Position3 const& edges(edge_lengths());
        for (Position3::size_type dim(0); dim < 3; ++dim)
        {
            const Real edge_length(edges[dim]);
            const Real diff(pos2[dim] - pos1[dim]), half(edge_length * 0.5);

            if (diff > half)
            {
                retval[dim] += edge_length;
            }
            else if (diff < -half)
            {
                retval[dim] -= edge_length;
            }
        }
        return retval;
    }

    inline Position3 apply_boundary(Position3 const& pos) const
    {
        return modulo(pos, edge_lengths());
    }

    Real distance_sq(
        Position3 const& pos1, Position3 const& pos2) const
    {
        Real retval(0);
        Position3 const& edges(edge_lengths());
        for (Position3::size_type dim(0); dim < 3; ++dim)
        {
            const Real edge_length(edges[dim]);
            const Real diff(pos2[dim] - pos1[dim]), half(edge_length * 0.5);

            if (diff > half)
            {
                retval += pow_2(diff - edge_length);
            }
            else if (diff < -half)
            {
                retval += pow_2(diff + edge_length);
            }
            else
            {
                retval += pow_2(diff);
            }
        }
        return retval;
    }

    inline Real distance(Position3 const& pos1, Position3 const& pos2) const
    {
        return std::sqrt(distance_sq(pos1, pos2));
    }

    virtual std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        Position3 const& pos, Real const& radius) const = 0;
    virtual std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        Position3 const& pos, Real const& radius,
        ParticleID const& ignore) const = 0;
    virtual std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        Position3 const& pos, Real const& radius,
        ParticleID const& ignore1, ParticleID const& ignore2) const = 0;
};

class ParticleSpaceVectorImpl
    : public ParticleSpace
{
public:

    typedef ParticleSpace::particle_container_type particle_container_type;

    ParticleSpaceVectorImpl(Position3 const& edge_lengths)
    {
        set_edge_lengths(edge_lengths);
    }

    Position3 const& edge_lengths() const
    {
        return edge_lengths_;
    }

    Integer num_particles() const;
    Integer num_particles(Species const& species) const;

    bool has_particle(ParticleID const& pid) const;
    bool update_particle(ParticleID const& pid, Particle const& p);
    void remove_particle(ParticleID const& pid);

    particle_container_type const& particles() const
    {
        return particles_;
    }

    void save(std::string const& filename)
    {
    	typedef struct h5_partcles {
    		int h5_particle_id;
    		double h5_particle_position[3];
    	} h5_particles;

    	h5_particles h5_p[particles_.size()];

    	for (int i=0; i<particles_.size(); i++){
    		h5_p[i].h5_particle_id = particles_[i].first;
    		h5_p[i].h5_particle_position[0] = particles_[i].second.position()[0];
    		h5_p[i].h5_particle_position[1] = particles_[i].second.position()[1];
    		h5_p[i].h5_particle_position[2] = particles_[i].second.position()[2];
    	}

    	H5::Exception::dontPrint();

    	H5File* file = new H5File( FILE_NAME, H5F_ACC_RDONLY );
    	CompType mtype( sizeof(h5_particles) );
        mtype.insertMember( MEMBER1, HOFFSET(h5_particles, h5_particle_id), PredType::NATIVE_INT);
        mtype.insertMember( MEMBER2, HOFFSET(h5_particles, h5_particle_position), PredType::NATIVE_DOUBLE);


    }

    std::pair<ParticleID, Particle> get_particle(ParticleID const& pid) const;
    std::vector<std::pair<ParticleID, Particle> > list_particles() const;
    std::vector<std::pair<ParticleID, Particle> >
    list_particles(Species const& species) const;

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        Position3 const& pos, Real const& radius) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        Position3 const& pos, Real const& radius,
        ParticleID const& ignore) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        Position3 const& pos, Real const& radius,
        ParticleID const& ignore1, ParticleID const& ignore2) const;

private:

    void set_edge_lengths(Position3 const& edge_lengths);

protected:

    typedef utils::get_mapper_mf<
        ParticleID, particle_container_type::size_type>::type particle_map_type;

    Position3 edge_lengths_;
    particle_container_type particles_;
    particle_map_type index_map_;
};

} // ecell4

#endif /* __PARTICLE_SPACE_HPP */
