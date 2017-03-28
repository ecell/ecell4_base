#ifndef __ECELL4_BD_BD_WORLD_HPP
#define __ECELL4_BD_BD_WORLD_HPP

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <sstream>

#include <ecell4/core/extras.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/ParticleSpace.hpp>
#include <ecell4/core/ParticleSpaceCellListImpl.hpp>
#include <ecell4/core/Model.hpp>
#include "BDPolygon.hpp"
#include "BDContainer2D.hpp"


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
    : public Space
{
public:

    typedef MoleculeInfo molecule_info_type;
    typedef ParticleSpaceCellListImpl particle_space_type;
    // typedef ParticleSpaceVectorImpl particle_space_type;
    typedef particle_space_type::particle_container_type particle_container_type;

    typedef BDPolygon::triangle_type    face_type;
    typedef BDPolygon::face_id_type face_id_type;

public:

    BDWorld(const Real3& edge_lengths = Real3(1, 1, 1),
        const Integer3& matrix_sizes = Integer3(3, 3, 3))
        : ps3d_(new particle_space_type(edge_lengths, matrix_sizes)),
          ps2d_(new ParticleContainer2D(edge_lengths))
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    BDWorld(
        const Real3& edge_lengths, const Integer3& matrix_sizes,
        boost::shared_ptr<RandomNumberGenerator> rng)
        : ps3d_(new particle_space_type(edge_lengths, matrix_sizes)),
          ps2d_(new ParticleContainer2D(edge_lengths)), rng_(rng)
    {
        ;
    }

    BDWorld(const std::string& filename)
        : ps3d_(new particle_space_type(Real3(1, 1, 1)))
    {
        throw NotImplemented("2D does not support HDF5 yet");
//         rng_ = boost::shared_ptr<RandomNumberGenerator>(
//             new GSLRandomNumberGenerator());
//         this->load(filename);
    }

    /**
     * @brief create and add a new particle
     * XXX checking overlap MUST be done for both 3D particles and surfaces
     *     using 3D distance (TODO for surfaces)
     * @param p a particle
     * @return a pair of a pair of pid (a particle id) and p (a particle)
     * and bool (if it's succeeded or not)
     */
    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Particle& p)
    {
        ParticleID pid(pidgen_());
        // if (has_particle(pid)) throw AlreadyExists("particle already exists");

        // TODO search surfaces

        if (list_particles_within_radius(p.position(), p.radius()).size() == 0)
        { //XXX: DONOT call this->update_particle
            (*ps3d_).update_particle(pid, p);
            return std::make_pair(std::make_pair(pid, p), true);
        }
        else
        {
            return std::make_pair(std::make_pair(pid, p), false);
        }
    }

    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Species& sp, const Real3& pos)
    {
        const MoleculeInfo info(get_molecule_info(sp));
        return new_particle(Particle(sp, pos, info.radius, info.D));
    }

    /**
     * @brief create and add a new particle on a surface.
     * @param p a particle
     * @param fid ID of the face that you want to put the particle
     * @return a pair of a pair of pid (a particle id) and p (a particle)
     * and bool (if it's succeeded or not)
     */
    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Particle& p, const face_id_type& fid)
    {
        const ParticleID pid(pidgen_());
        // if (has_particle(pid)) throw AlreadyExists("particle already exists");

        if(list_particles_within_radius(
                std::make_pair(p.position(), fid), p.radius()).size() == 0)
        {
            (*ps2d_).update_particle(pid, p, fid);
            return std::make_pair(std::make_pair(pid, p), true);
        }
        else
        {
            return std::make_pair(std::make_pair(pid, p), false);
        }
    }

    std::pair<std::pair<ParticleID, Particle>, bool>
    new_particle(const Species& sp, const Real3& pos, const face_id_type& fid)
    {
        const MoleculeInfo info(get_molecule_info(sp));
        return new_particle(Particle(sp, pos, info.radius, info.D), fid);
    }

    /**
     * draw attributes of species and return it as a molecule info.
     * @param sp a species
     * @return info a molecule info
     */
    MoleculeInfo get_molecule_info(const Species& sp) const
    {
        Real radius(0.0), D(0.0);

        if (sp.has_attribute("radius") && sp.has_attribute("D"))
        {
            radius = std::atof(sp.get_attribute("radius").c_str());
            D = std::atof(sp.get_attribute("D").c_str());
        }
        else if (boost::shared_ptr<Model> bound_model = lock_model())
        {
            Species attributed(bound_model->apply_species_attributes(sp));
            if (attributed.has_attribute("radius")
                && attributed.has_attribute("D"))
            {
                radius = std::atof(
                    attributed.get_attribute("radius").c_str());
                D = std::atof(attributed.get_attribute("D").c_str());
            }
        }

        MoleculeInfo info = {radius, D};
        return info;
    }

    // SpaceTraits

    const Real t() const
    {
        // check synchro
        assert((*ps3d_).t() == (*ps2d_).t());
        return (*ps3d_).t();
    }

    void set_t(const Real& t)
    {
        (*ps2d_).set_t(t);
        (*ps3d_).set_t(t);
    }

    // ParticleSpaceTraits

    const Real3& edge_lengths() const
    {
        // check
//         assert((*ps3d_).edge_lengths() == (*ps2d_).edge_lengths());
        return (*ps3d_).edge_lengths();
    }

    Integer num_particles() const
    {
        return (*ps2d_).num_particles() + (*ps3d_).num_particles();
    }

    Integer num_particles(const Species& species) const
    {
        return (*ps2d_).num_particles(species) + (*ps3d_).num_particles(species);
    }

    Integer num_particles_exact(const Species& species) const
    {
        return (*ps2d_).num_particles_exact(species) +
               (*ps3d_).num_particles_exact(species);
    }

    bool has_particle(const ParticleID& pid) const
    {
        return (*ps2d_).has_particle(pid) || (*ps3d_).has_particle(pid);
    }

    std::vector<std::pair<ParticleID, Particle> > list_particles() const
    {
        return concat((*ps3d_).list_particles(), (*ps2d_).list_particles());
    }

    std::vector<std::pair<ParticleID, Particle> >
    list_particles(const Species& sp) const
    {
        return concat((*ps3d_).list_particles(sp), (*ps2d_).list_particles(sp));
    }

    std::vector<std::pair<ParticleID, Particle> >
    list_particles_exact(const Species& sp) const
    {
        return concat((*ps3d_).list_particles_exact(sp),
                      (*ps2d_).list_particles_exact(sp));
    }

    std::vector<Species> list_species() const
    {
        return concat((*ps3d_).list_species(), (*ps2d_).list_species());
    }

    virtual Real get_value(const Species& sp) const
    {
        return static_cast<Real>(num_molecules(sp));
    }

    virtual Real get_value_exact(const Species& sp) const
    {
        return static_cast<Real>(num_molecules_exact(sp));
    }

    // ParticleSpace member functions

    bool update_particle_without_checking(const ParticleID& pid, const Particle& p)
    {
        return (*ps3d_).update_particle(pid, p);
    }

    bool update_particle_without_checking(
            const ParticleID& pid, const Particle& p, const face_id_type& fid)
    {
        return (*ps2d_).update_particle(pid, p, fid);
    }

    bool update_particle(const ParticleID& pid, const Particle& p)
    {
        if (list_particles_within_radius(p.position(), p.radius(), pid).size()
            == 0)
        {
            return (*ps3d_).update_particle(pid, p);
        }
        else
        {
            return true;
        }
    }

    bool update_particle(
            const ParticleID& pid, const Particle& p, const face_id_type& fid)
    {
        // XXX: checking overlap in 3D spherical region
        if (list_particles_within_radius(std::make_pair(p.position(), fid),
                    p.radius(), pid).size() == 0)
        {
            return (*ps2d_).update_particle(pid, p, fid);
        }
        else
        {
            return true;
        }
    }

    std::pair<ParticleID, Particle>
    get_particle(const ParticleID& pid) const
    {
        if((*ps3d_).has_particle(pid))
            return (*ps3d_).get_particle(pid);
        else if((*ps2d_).has_particle(pid))
            return (*ps2d_).get_particle(pid);
        else
            throw NotFound("BDWorld: no such particle");
    }

    void remove_particle(const ParticleID& pid)
    {
        if((*ps3d_).has_particle(pid))
            (*ps3d_).remove_particle(pid);
        else if((*ps2d_).has_particle(pid))
            (*ps2d_).remove_particle(pid);
        else
            throw NotFound("BDWorld: no such particle");
    }

    // list-up all the particles using 3D distance.
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Real3& pos, const Real& radius) const
    {
        return concat((*ps3d_).list_particles_within_radius(pos, radius),
                      (*ps2d_).list_particles_within_radius(pos, radius));
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Real3& pos, const Real& radius, const ParticleID& ignore) const
    {
        return concat(
            (*ps3d_).list_particles_within_radius(pos, radius, ignore),
            (*ps2d_).list_particles_within_radius(pos, radius, ignore)
            );
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Real3& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        return concat(
            (*ps3d_).list_particles_within_radius(pos, radius, ignore1, ignore2),
            (*ps2d_).list_particles_within_radius(pos, radius, ignore1, ignore2)
            );
    }

    // list-up 2D particles using 2D distance(along surface).
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real& radius) const
    {
        return (*ps2d_).list_particles_within_radius(pos, radius);
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real& radius,
        const ParticleID& ignore) const
    {
        return (*ps2d_).list_particles_within_radius(pos, radius, ignore);
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        return (*ps2d_).list_particles_within_radius(
                pos, radius, ignore1, ignore2);
    }

    inline Real3 periodic_transpose(
        const Real3& pos1, const Real3& pos2) const
    {
        //XXX: only use 3d.
        //     assuming 2d has same space size(assured by a ctor, maybe)
        return (*ps3d_).periodic_transpose(pos1, pos2);
    }

    inline Real3 apply_boundary(const Real3& pos) const
    {
        //XXX: only use 3d.
        //     assuming 2d has same space size(assured by a ctor, maybe)
        return (*ps3d_).apply_boundary(pos);
    }

    inline std::pair<Real3, face_id_type>
    apply_surface(const std::pair<Real3, face_id_type>& position,
                  const Real3& displacement) const
    {
        return (*ps2d_).apply_surface(position, displacement);
    }

    //! @brief square distance in normal 3 dimensional space
    inline Real distance_sq(const Real3& pos1, const Real3& pos2) const
    {
        return (*ps3d_).distance_sq(pos1, pos2);
    }

    //! @brief distance in normal 3 dimensional space
    inline Real distance(const Real3& pos1, const Real3& pos2) const
    {
        return (*ps3d_).distance(pos1, pos2);
    }

    //! @brief square distance along the surface. If inaccessible, return infty.
    inline Real distance_sq(const Real3& pos1, const face_id_type& fid1,
                            const Real3& pos2, const face_id_type& fid2) const
    {
        return (*ps2d_).distance_sq(pos1, fid1, pos2, fid2);
    }

    //! @brief distance along the surface. If inaccessible, return infty.
    inline Real distance(const Real3& pos1, const face_id_type& fid1,
                         const Real3& pos2, const face_id_type& fid2) const
    {
        return (*ps2d_).distance(pos1, fid1, pos2, fid2);
    }

    // CompartmentSpaceTraits

    Integer num_molecules(const Species& sp) const
    {
        return (*ps3d_).num_molecules(sp) + (*ps2d_).num_molecules(sp);
    }

    Integer num_molecules_exact(const Species& sp) const
    {
        return (*ps3d_).num_molecules_exact(sp) +
               (*ps2d_).num_molecules_exact(sp);
    }

    //XXX only 3D-molecules are available yet.
    void add_molecules(const Species& sp, const Integer& num)
    {
        extras::throw_in_particles(*this, sp, num, rng());
    }

    //XXX only 3D-molecules are available yet.
    void add_molecules(const Species& sp, const Integer& num,
                       const boost::shared_ptr<Shape> shape)
    {
        extras::throw_in_particles(*this, sp, num, shape, rng());
    }

    //XXX both 3D and 2D particles are randomly removed.
    void remove_molecules(const Species& sp, const Integer& num)
    {
        if (num < 0)
        {
            throw std::invalid_argument(
                "The number of molecules must be positive.");
        }

        std::vector<std::pair<ParticleID, Particle> >
            particles(list_particles(sp)); // 2D + 3D concatenated.

        const Integer num_particles(particles.size());
        if (num_particles < num)
        {
            throw std::invalid_argument(
                "The number of molecules cannot be negative.");
        }

        shuffle((*rng_), particles);
        for (std::vector<std::pair<ParticleID, Particle> >::const_iterator
            i(particles.begin()); i != particles.begin() + num; ++i)
        {
            remove_particle((*i).first);
        }
    }

    const Real volume() const
    {
        // XXX: assuming 2d has same space size(assured by a ctor, maybe)
        const Real3& lengths(edge_lengths());
        return lengths[0] * lengths[1] * lengths[2];
    }

    // Optional members

    inline boost::shared_ptr<RandomNumberGenerator>& rng()
    {
        return rng_;
    }

    // XXX return 3D particles ONLY!!!
    const particle_container_type& particles() const
    {
        return (*ps3d_).particles();
    }

    const particle_container_type& particles_3d() const
    {
        return (*ps3d_).particles();
    }

    const particle_container_type& particles_2d() const
    {
        return (*ps2d_).particles();
    }

    void save(const std::string& filename) const
    {
#ifdef WITH_HDF5
        boost::scoped_ptr<H5::H5File>
            fout(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
        rng_->save(fout.get());
        pidgen_.save(fout.get());
        boost::scoped_ptr<H5::Group>
            group(new H5::Group(fout->createGroup("ParticleSpace")));
        ps3d_->save_hdf5(group.get());
        extras::save_version_information(fout.get(), std::string("ecell4-bd-") + std::string(ECELL4_VERSION));
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

        const std::string required = "ecell4-bd-4.1.0";
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

        const H5::Group group(fin->openGroup("ParticleSpace"));
        ps3d_->load_hdf5(group);
        pidgen_.load(*fin);
        rng_->load(*fin);
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
                std::cerr << "Warning: Model already bound to BDWorld"
                    << std::endl;
            }
        }

        model_ = model;
    }

    boost::shared_ptr<Model> lock_model() const
    {
        return model_.lock();
    }

    // XXX
    ParticleSpace&       container_3D(){return *ps3d_;}
    ParticleContainer2D& container_2D(){return *ps2d_;}

    // TODO make return type Shape*
//     BDPolygon const& polygon() const {return (*ps2d_).polygon();}

    void set_polygon(const BDPolygon& poly)
    {
        (*ps2d_).polygon() = poly;
        (*ps2d_).setup_polygon();
        return;
    }

    const face_type& belonging_face(const ParticleID& pid)
    {
        return (*ps2d_).belonging_face(pid);
    }

    const face_id_type& belonging_faceid(const ParticleID& pid)
    {
        return (*ps2d_).belonging_faceid(pid);
    }

    face_type const& get_face(const std::size_t i)
    {
        return (*ps2d_).polygon().triangle_at(face_id_type(i));
    }

    Real3 get_inter_position_vector(const Real3& lhs, const Real3& rhs) const
    {
        return rhs - lhs;
    }

    Real3 get_inter_position_vector(
            const std::pair<Real3, face_id_type>& lhs,
            const std::pair<Real3, face_id_type>& rhs) const
    {
        return (*ps2d_).get_inter_position_vector(lhs, rhs);
    }

private:

    template<typename T>
    std::vector<T> concat(std::vector<T> lhs, std::vector<T> rhs) const
    {
        if(rhs.empty()) return lhs;

        lhs.reserve(lhs.size() + rhs.size());
        std::copy(rhs.begin(), rhs.end(), std::back_inserter(lhs));
        return lhs;
    }

protected:

    boost::scoped_ptr<ParticleSpace>       ps3d_;
    boost::scoped_ptr<ParticleContainer2D> ps2d_;

    boost::shared_ptr<RandomNumberGenerator> rng_;
    SerialIDGenerator<ParticleID> pidgen_;

    boost::weak_ptr<Model> model_;
};

} // bd

} // ecell4

#endif /* __ECELL4_BD_BD_WORLD_HPP */
