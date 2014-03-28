#ifndef __ECELL4_BD_BD_WORLD_HPP
#define __ECELL4_BD_BD_WORLD_HPP

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include <ecell4/core/extras.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/ParticleSpace.hpp>
#include <ecell4/core/Model.hpp>


namespace ecell4
{

namespace bd
{

struct MoleculeInfo
{
    const Real radius;
    const Real D;
};

class BDWorld
{
public:

    typedef MoleculeInfo molecule_info_type;
    typedef ParticleSpace::particle_container_type particle_container_type;

public:

    BDWorld(
        const Position3& edge_lengths,
        boost::shared_ptr<RandomNumberGenerator> rng)
        : ps_(new ParticleSpaceVectorImpl(edge_lengths)), rng_(rng)
    {
        ;
    }

    BDWorld(const Position3& edge_lengths)
        : ps_(new ParticleSpaceVectorImpl(edge_lengths))
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    /**
     * create and add a new particle
     * @param p a particle
     * @return pid a particle id
     */
    ParticleID new_particle(const Particle& p)
    {
        ParticleID pid(pidgen_());
        // if (has_particle(pid))
        // {
        //     throw AlreadyExists("particle already exists");
        // }
        (*ps_).update_particle(pid, p);
        return pid;
    }

    /**
     * draw attributes of species and return it as a molecule info.
     * @param sp a species
     * @return info a molecule info
     */
    MoleculeInfo get_molecule_info(const Species& sp) const
    {
        Real radius, D;
        if (sp.has_attribute("radius") && sp.has_attribute("D")) 
        {
            radius = std::atof(sp.get_attribute("radius").c_str());
            D = std::atof(sp.get_attribute("D").c_str());
        } 
        else if (boost::shared_ptr<Model> bound_model = this->model_.lock())
        {
            if (bound_model->has_species_attribute(sp)) 
            {
                Species sp_model(bound_model->apply_species_attributes(sp));
                if (sp_model.has_attribute("radius") && sp_model.has_attribute("D"))
                {
                    radius = std::atof(sp_model.get_attribute("radius").c_str());
                    D = std::atof(sp_model.get_attribute("D").c_str());
                }
            }
        } 
        MoleculeInfo info = {radius, D};
        return info;
    }

    // SpaceTraits

    const Real& t() const
    {
        return (*ps_).t();
    }

    void set_t(const Real& t)
    {
        (*ps_).set_t(t);
    }

    // ParticleSpaceTraits

    const Position3& edge_lengths() const
    {
        return (*ps_).edge_lengths();
    }

    Integer num_particles() const
    {
        return (*ps_).num_particles();
    }

    Integer num_particles(const Species& species) const
    {
        return (*ps_).num_particles(species);
    }

    bool has_particle(const ParticleID& pid) const
    {
        return (*ps_).has_particle(pid);
    }

    std::vector<std::pair<ParticleID, Particle> > list_particles() const
    {
        return (*ps_).list_particles();
    }

    std::vector<std::pair<ParticleID, Particle> >
    list_particles(const Species& species) const
    {
        return (*ps_).list_particles(species);
    }

    // ParticleSpace member functions

    bool update_particle(const ParticleID& pid, const Particle& p)
    {
        return (*ps_).update_particle(pid, p);
    }

    std::pair<ParticleID, Particle>
    get_particle(const ParticleID& pid) const
    {
        return (*ps_).get_particle(pid);
    }

    void remove_particle(const ParticleID& pid)
    {
        (*ps_).remove_particle(pid);
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Position3& pos, const Real& radius) const
    {
        return (*ps_).list_particles_within_radius(pos, radius);
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Position3& pos, const Real& radius, const ParticleID& ignore) const
    {
        return (*ps_).list_particles_within_radius(pos, radius, ignore);
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Position3& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        return (*ps_).list_particles_within_radius(pos, radius, ignore1, ignore2);
    }

    inline Position3 periodic_transpose(
        const Position3& pos1, const Position3& pos2) const
    {
        return (*ps_).periodic_transpose(pos1, pos2);
    }

    inline Position3 apply_boundary(const Position3& pos) const
    {
        return (*ps_).apply_boundary(pos);
    }

    inline Real distance_sq(const Position3& pos1, const Position3& pos2) const
    {
        return (*ps_).distance_sq(pos1, pos2);
    }

    inline Real distance(const Position3& pos1, const Position3& pos2) const
    {
        return (*ps_).distance(pos1, pos2);
    }

    // CompartmentSpaceTraits

    Integer num_molecules(const Species& sp) const
    {
        return num_particles(sp);
    }

    void add_molecules(const Species& sp, const Integer& num)
    {
        extras::throw_in_particles(*this, sp, num, *rng());
    }

    const Real volume() const
    {
        const Position3& lengths(edge_lengths());
        return lengths[0] * lengths[1] * lengths[2];
    }

    // Optional members

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return rng_;
    }

    const particle_container_type& particles() const
    {
        return (*ps_).particles();
    }

    void save(const std::string& filename) const
    {
        boost::scoped_ptr<H5::H5File>
            fout(new H5::H5File(filename, H5F_ACC_TRUNC));
        rng_->save(fout.get());
        pidgen_.save(fout.get());
        boost::scoped_ptr<H5::Group>
            group(new H5::Group(fout->createGroup("ParticleSpace")));
        ps_->save(group.get());
    }

    void load(const std::string& filename)
    {
        boost::scoped_ptr<H5::H5File>
            fin(new H5::H5File(filename, H5F_ACC_RDONLY));
        const H5::Group group(fin->openGroup("ParticleSpace"));
        ps_->load(group);
        pidgen_.load(*fin);
        rng_->load(*fin);
    }

    void bind_to(boost::shared_ptr<Model> model)
    {
        this->model_ = model;
    }

protected:

    boost::scoped_ptr<ParticleSpace> ps_;
    boost::shared_ptr<RandomNumberGenerator> rng_;
    SerialIDGenerator<ParticleID> pidgen_;
    boost::weak_ptr<Model> model_;
};

} // bd

} // ecell4

#endif /* __ECELL4_BD_BD_WORLD_HPP */
