#ifndef __ECELL4_MESO_MESOSCOPIC_WORLD_HPP
#define __ECELL4_MESO_MESOSCOPIC_WORLD_HPP

#include <numeric>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/SubvolumeSpace.hpp>
#include <ecell4/core/Model.hpp>

namespace ecell4
{

namespace meso
{

struct MoleculeInfo
{
    const Real D;
};

class MesoscopicWorld
    : public Space
{
public:

    typedef SubvolumeSpace::coordinate_type coordinate_type;

public:

    MesoscopicWorld(const std::string& filename)
        : cs_(new SubvolumeSpaceVectorImpl(Position3(1, 1, 1), Global(1, 1, 1)))
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        this->load(filename);
    }

    MesoscopicWorld(const Position3& edge_lengths = Position3(1, 1, 1))
        : cs_(new SubvolumeSpaceVectorImpl(edge_lengths, Global(1, 1, 1)))
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    MesoscopicWorld(const Position3& edge_lengths,
        const Global& matrix_sizes, boost::shared_ptr<RandomNumberGenerator> rng)
        : cs_(new SubvolumeSpaceVectorImpl(edge_lengths, matrix_sizes)), rng_(rng)
    {
        ;
    }

    MesoscopicWorld(const Position3& edge_lengths, const Global& matrix_sizes)
        : cs_(new SubvolumeSpaceVectorImpl(edge_lengths, matrix_sizes))
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(
            new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

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
                std::cerr << "Warning: Model already bound to GillespieWorld."
                    << std::endl;
            }
        }
        this->model_ = model;
    }

    void save(const std::string& filename) const
    {
        boost::scoped_ptr<H5::H5File>
            fout(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
        rng_->save(fout.get());
        boost::scoped_ptr<H5::Group>
            group(new H5::Group(fout->createGroup("SubvolumeSpace")));
        cs_->save(group.get());
    }

    void load(const std::string& filename)
    {
        boost::scoped_ptr<H5::H5File>
            fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));
        rng_->load(*fin);
        const H5::Group group(fin->openGroup("SubvolumeSpace"));
        cs_->load(group);
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
    const Position3 subvolume_edge_lengths() const;
    const Position3& edge_lengths() const;

    const Global matrix_sizes() const
    {
        return cs_->matrix_sizes();
    }

    coordinate_type global2coord(const Global& g) const;
    Global coord2global(const coordinate_type& c) const;
    coordinate_type get_neighbor(const coordinate_type& c, const Integer rnd) const;

    Real get_value(const Species& sp) const;
    Real get_value_exact(const Species& sp) const;
    Integer num_molecules(const Species& sp) const;
    Integer num_molecules_exact(const Species& sp) const;
    Integer num_molecules(const Species& sp, const coordinate_type& c) const;
    Integer num_molecules_exact(const Species& sp, const coordinate_type& c) const;
    void add_molecules(const Species& sp, const Integer& num, const coordinate_type& c);
    void remove_molecules(const Species& sp, const Integer& num, const coordinate_type& c);

    Integer num_molecules(const Species& sp, const Global& g) const
    {
        return cs_->num_molecules(sp, g);
    }

    Integer num_molecules_exact(const Species& sp, const Global& g) const
    {
        return cs_->num_molecules_exact(sp, g);
    }

    void add_molecules(const Species& sp, const Integer& num, const Global& g)
    {
        cs_->add_molecules(sp, num, g);
    }

    void remove_molecules(const Species& sp, const Integer& num, const Global& g)
    {
        cs_->remove_molecules(sp, num, g);
    }

    void add_molecules(const Species& sp, const Integer& num)
    {
        for (Integer i(0); i < num; ++i)
        {
            cs_->add_molecules(sp, 1, rng_->uniform_int(0, num_subvolumes() - 1));
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
