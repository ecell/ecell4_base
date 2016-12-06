#ifndef ECELL4_BD_BD_CONTAINER_2D
#define ECELL4_BD_BD_CONTAINER_2D

#include <ecell4/core/ParticleSpace.hpp>
#include "BDPolygon.hpp"
#include <set>

namespace ecell4
{

namespace bd
{

class ParticleContainer2D : public ParticleSpace
{
public:
    typedef ParticleSpace base_type;
    typedef typename base_type::particle_container_type particle_container_type;
    typedef typename particle_container_type::size_type container_index_type;
    typedef utils::get_mapper_mf<ParticleID, container_index_type>::type 
        pid_to_particle_index_type;
    typedef std::set<ParticleID> particle_id_set;
    typedef std::map<Species::serial_type, particle_id_set>
        per_species_particle_id_set;

    typedef Polygon polygon_type;

public:

    ParticleContainer2D();
    ~ParticleContainer2D(){}

    const Real3& edge_lengths() const {return edge_lengths_;}
    Integer num_particles() const;
    Integer num_particles(const Species& sp) const;
    Integer num_particles_exact(const Species& sp) const;

    Integer num_molecules(const Species& sp) const;
    Integer num_molecules_exact(const Species& sp) const;

    std::vector<std::pair<ParticleID, Particle> >
        list_particles() const;
    std::vector<std::pair<ParticleID, Particle> >
        list_particles(const Species& sp) const;
    std::vector<std::pair<ParticleID, Particle> >
        list_particles_exact(const Species& sp) const;
 
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        list_particles_within_radius(
            const Real3& pos, const Real& radius) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        list_particles_within_radius(
            const Real3& pos, const Real& radius,
            const ParticleID& ignore) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        list_particles_within_radius(
            const Real3& pos, const Real& radius,
            const ParticleID& ignore1, const ParticleID& ignore2) const;

private:

    Real3                       edge_lengths_;
    polygon_type                polygon_;
    pid_to_particle_index_type  rmap_;
    per_species_particle_id_set particle_pool_;
    particle_container_type     particles_;
};


}// bd
}// ecell4
#endif /* ECELL4_BD_BD_CONTAINER_2D */
