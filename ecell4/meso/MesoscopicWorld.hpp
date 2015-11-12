#ifndef __ECELL4_MESO_MESOSCOPIC_WORLD_HPP
#define __ECELL4_MESO_MESOSCOPIC_WORLD_HPP

#include <numeric>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/SubvolumeSpace.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/Shape.hpp>
#include <ecell4/core/extras.hpp>

namespace ecell4
{

namespace meso
{

struct MoleculeInfo
{
    const Real D;
    const std::string loc;
};

class MesoscopicWorld
    : public Space
{
public:

    typedef SubvolumeSpace::coordinate_type coordinate_type;
    typedef MoleculeInfo molecule_info_type;

public:

    MesoscopicWorld(const std::string& filename)
        : cs_(new SubvolumeSpaceVectorImpl(Real3(1, 1, 1), Integer3(1, 1, 1)))
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        this->load(filename);
    }

    MesoscopicWorld(const Real3& edge_lengths = Real3(1, 1, 1))
        : cs_(new SubvolumeSpaceVectorImpl(edge_lengths, Integer3(1, 1, 1)))
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    MesoscopicWorld(const Real3& edge_lengths, const Integer3& matrix_sizes)
        : cs_(new SubvolumeSpaceVectorImpl(edge_lengths, matrix_sizes))
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    MesoscopicWorld(const Real3& edge_lengths,
        const Integer3& matrix_sizes, boost::shared_ptr<RandomNumberGenerator> rng)
        : cs_(new SubvolumeSpaceVectorImpl(edge_lengths, matrix_sizes)), rng_(rng)
    {
        ;
    }

    MesoscopicWorld(const Real3& edge_lengths, const Real subvolume_length);
    MesoscopicWorld(
        const Real3& edge_lengths, const Real subvolume_length,
        boost::shared_ptr<RandomNumberGenerator> rng);

    virtual ~MesoscopicWorld()
    {
        ;
    }

    void bind_to(boost::shared_ptr<Model> model)
    {
        if (boost::shared_ptr<Model> bound_model = lock_model())
        {
            if (bound_model.get() != model.get())
            {
                std::cerr << "Warning: Model already bound to MesoscopicWorld."
                    << std::endl;
            }
        }

        this->model_ = model;
    }

    void save(const std::string& filename) const
    {
#ifdef WITH_HDF5
        boost::scoped_ptr<H5::H5File>
            fout(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
        rng_->save(fout.get());
        boost::scoped_ptr<H5::Group>
            group(new H5::Group(fout->createGroup("SubvolumeSpace")));
        cs_->save(group.get());
        extras::save_version_information(fout.get(), "ecell4-meso-0.0-1");
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
        const H5::Group group(fin->openGroup("SubvolumeSpace"));
        cs_->load(group);
#else
        throw NotSupported("not supported yet.");
#endif
    }

    boost::shared_ptr<Model> lock_model() const
    {
        return model_.lock();
    }

    inline const boost::shared_ptr<RandomNumberGenerator>& rng()
    {
        return rng_;
    }

    MoleculeInfo get_molecule_info(const Species& sp) const;

    const Real& t() const;
    void set_t(const Real& t);
    const Integer num_subvolumes() const;
    const Real subvolume() const;
    const Real volume() const;
    const Real3 subvolume_edge_lengths() const;
    const Real3& edge_lengths() const;

    const Integer num_subvolumes(const Species& sp) const
    {
        return cs_->num_subvolumes(sp);
    }

    const Integer3 matrix_sizes() const
    {
        return cs_->matrix_sizes();
    }

    coordinate_type global2coord(const Integer3& g) const;
    Integer3 coord2global(const coordinate_type& c) const;
    Integer3 position2global(const Real3& pos) const;
    Integer position2coordinate(const Real3& pos) const;

    coordinate_type get_neighbor(const coordinate_type& c, const Integer rnd) const
    {
        return cs_->get_neighbor(c, rnd);
    }

    Real get_value(const Species& sp) const;
    Real get_value_exact(const Species& sp) const;
    Integer num_molecules(const Species& sp) const;
    Integer num_molecules_exact(const Species& sp) const;
    Integer num_molecules(const Species& sp, const coordinate_type& c) const;
    Integer num_molecules_exact(const Species& sp, const coordinate_type& c) const;
    void add_molecules(const Species& sp, const Integer& num, const coordinate_type& c);
    void remove_molecules(const Species& sp, const Integer& num, const coordinate_type& c);

    Integer num_molecules(const Species& sp, const Integer3& g) const
    {
        return cs_->num_molecules(sp, g);
    }

    Integer num_molecules_exact(const Species& sp, const Integer3& g) const
    {
        return cs_->num_molecules_exact(sp, g);
    }

    void add_molecules(const Species& sp, const Integer& num, const Integer3& g)
    {
        cs_->add_molecules(sp, num, g);
    }

    void remove_molecules(const Species& sp, const Integer& num, const Integer3& g)
    {
        cs_->remove_molecules(sp, num, g);
    }

    std::vector<coordinate_type> list_coordinates(const Species& sp) const;
    std::vector<coordinate_type> list_coordinates_exact(const Species& sp) const;

    void add_molecules(const Species& sp, const Integer& num)
    {
        const molecule_info_type minfo(get_molecule_info(sp));
        if (minfo.loc == "")
        {
            for (Integer i(0); i < num; ++i)
            {
                cs_->add_molecules(sp, 1, rng_->uniform_int(0, num_subvolumes() - 1));
            }
            return;
        }

        const Species st(minfo.loc);
        if (!cs_->has_structure(st))
        {
            throw NotFound("no space to throw-in.");
        }

        Integer i(0);
        while (i < num)
        {
            const coordinate_type j(rng_->uniform_int(0, num_subvolumes() - 1));
            if (cs_->check_structure(minfo.loc, j))
            {
                cs_->add_molecules(sp, 1, j);
                i++;
            }
        }
    }

    void add_molecules(const Species& sp, const Integer& num,
        const boost::shared_ptr<Shape> shape)
    {
        const molecule_info_type minfo(get_molecule_info(sp));
        if (minfo.loc == "")
        {
            for (Integer i(0); i < num; ++i)
            {
                const Real3 pos(shape->draw_position(rng_));
                const coordinate_type& coord(
                    cs_->global2coord(cs_->position2global(pos)));
                cs_->add_molecules(sp, 1, coord);
            }
            return;
        }

        const Species st(minfo.loc);
        if (!cs_->has_structure(st))
        {
            throw NotFound("no space to throw-in.");
        }

        Integer i(0);
        while (i < num)
        {
            const Real3 pos(shape->draw_position(rng_));
            const Integer3 g(cs_->position2global(pos));
            const coordinate_type j(cs_->global2coord(g));
            if (cs_->check_structure(minfo.loc, j))
            {
                cs_->add_molecules(sp, 1, j);
                i++;
            }
        }
    }

    void remove_molecules(const Species& sp, const Integer& num)
    {
        std::vector<Integer> a(num_subvolumes());
        for (coordinate_type c(0); c < num_subvolumes(); ++c)
        {
            a[c] = num_molecules_exact(sp, c);
        }

        Integer num_tot(std::accumulate(a.begin(), a.end(), 0));
        if (num_tot < num)
        {
            std::ostringstream message;
            message << "The number of molecules cannot be negative. [" << sp.serial() << "]";
            throw std::invalid_argument(message.str());
        }

        for (Integer i(0); i < num; ++i)
        {
            const Integer rnd1(rng_->uniform_int(0, num_tot - 1));
            Integer acct(0);
            for (coordinate_type c(0); c < num_subvolumes(); ++c)
            {
                acct += a[c];
                if (acct > rnd1)
                {
                    cs_->remove_molecules(sp, 1, c);
                    a[c] -= 1;
                    --num_tot;
                    break;
                }
            }
        }
    }

    std::pair<std::pair<ParticleID, Particle>, bool> new_particle(const Particle& p)
    {
        add_molecules(p.species(), 1, position2coordinate(p.position()));
        return std::make_pair(std::make_pair(ParticleID(), p), true);
    }

    std::pair<std::pair<ParticleID, Particle>, bool> new_particle(
        const Species& sp, const Real3& pos)
    {
        add_molecules(sp, 1, position2coordinate(pos));
        const MoleculeInfo info(get_molecule_info(sp));
        return std::make_pair(
            std::make_pair(ParticleID(), Particle(sp, pos, 0.0, info.D)), true);
    }

    void add_structure(const Species& sp, const boost::shared_ptr<const Shape>& shape);
    bool on_structure(const Species& sp, const coordinate_type& coord) const;

    bool has_structure(const Species& sp) const
    {
        return cs_->has_structure(sp);
    }

    Real get_occupancy(const Species& sp, const coordinate_type& coord) const
    {
        return cs_->get_occupancy(sp, coord);
    }

    Real get_occupancy(const Species& sp, const Integer3& g) const
    {
        return cs_->get_occupancy(sp, g);
    }

    // Shape::dimension_kind get_dimension(const Species& sp) const
    // {
    //     return cs_->get_dimension(sp);
    // }

    inline bool on_structure(const Species& sp, const Integer3& g) const
    {
        return on_structure(sp, global2coord(g));
    }

    bool check_structure(const Species& sp, const Integer3& g) const
    {
        return cs_->check_structure(sp, g);
    }

    Real get_volume(const Species& sp) const;

    const std::vector<Species>& species() const;
    std::vector<Species> list_species() const;

    std::vector<std::pair<ParticleID, Particle> > list_particles() const;
    std::vector<std::pair<ParticleID, Particle> > list_particles_exact(const Species& sp) const;
    std::vector<std::pair<ParticleID, Particle> > list_particles(const Species& sp) const;

private:

    boost::scoped_ptr<SubvolumeSpace> cs_;
    boost::shared_ptr<RandomNumberGenerator> rng_;

    boost::weak_ptr<Model> model_;
};

} // meso

} // ecell4

#endif /* __ECELL4_MESO_MESOSCOPIC_WORLD_HPP */
