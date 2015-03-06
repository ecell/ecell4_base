#ifndef __ECELL4_ODE_ODE_WORLD_HPP
#define __ECELL4_ODE_ODE_WORLD_HPP

#include <ecell4/core/Species.hpp>
#include <ecell4/core/Context.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Space.hpp>
#include <ecell4/core/Model.hpp>
#ifdef WITH_HDF5
#include <ecell4/core/CompartmentSpaceHDF5Writer.hpp>
#endif
#include <ecell4/core/Shape.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

namespace ecell4
{

namespace ode
{

template<typename Tspace_>
struct ODEWorldHDF5Traits
    : public CompartmentSpaceHDF5TraitsBase<Tspace_, H5DataTypeTraits_double>
{
    typedef CompartmentSpaceHDF5TraitsBase<Tspace_, H5DataTypeTraits_double> base_type;
    typedef typename base_type::num_molecules_type num_molecules_type;
    typedef typename base_type::space_type space_type;

    num_molecules_type getter(const space_type& space, const Species& sp) const
    {
        return space.get_value_exact(sp);
    }

    void setter(
        Tspace_& space, const Species& sp, const num_molecules_type& value) const
    {
        space.set_value(sp, value);
    }
};

class ODEWorld
    : public Space
{
protected:

    typedef std::vector<Real> num_molecules_container_type;
    typedef std::vector<Species> species_container_type;
    typedef utils::get_mapper_mf<
        Species, num_molecules_container_type::size_type>::type species_map_type;

public:

    ODEWorld(const Real3& edge_lengths = Real3(1, 1, 1))
        : t_(0.0)
    {
        reset(edge_lengths);
    }

    ODEWorld(const std::string& filename)
        : t_(0.0)
    {
        reset(Real3(1, 1, 1));
        this->load(filename);
    }

    // SpaceTraits

    const Real& t() const
    {
        return t_;
    }

    void set_t(const Real& t)
    {
        if (t < 0.0)
        {
            throw std::invalid_argument("the time must be positive.");
        }
        t_ = t;
    }

    const Real3& edge_lengths() const
    {
        return edge_lengths_;
    }

    void reset(const Real3& edge_lengths)
    {
        t_ = 0.0;
        index_map_.clear();
        num_molecules_.clear();
        species_.clear();

        for (Real3::size_type dim(0); dim < 3; ++dim)
        {
            if (edge_lengths[dim] <= 0)
            {
                throw std::invalid_argument("the edge length must be positive.");
            }
        }

        edge_lengths_ = edge_lengths;
        volume_ = edge_lengths[0] * edge_lengths[1] * edge_lengths[2];
    }

    const Real volume() const
    {
        return volume_;
    }

    void set_volume(const Real& volume)
    {
        if (volume <= 0.0)
        {
            throw std::invalid_argument("The volume must be positive.");
        }

        volume_ = volume;
        const Real L(cbrt(volume));
        edge_lengths_ = Real3(L, L, L);
    }

    // CompartmentSpaceTraits

    Integer num_molecules(const Species& sp) const
    {
        return static_cast<Integer>(get_value(sp));
    }

    Integer num_molecules_exact(const Species& sp) const
    {
        return static_cast<Integer>(get_value_exact(sp));
    }

    std::vector<Species> list_species() const
    {
        return species_;
    }

    // CompartmentSpace member functions

    void add_molecules(const Species& sp, const Real& num)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            reserve_species(sp);
            i = index_map_.find(sp);
        }

        num_molecules_[(*i).second] += num;
    }

    void remove_molecules(const Species& sp, const Real& num)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            throw NotFound("Species not found");
        }

        num_molecules_[(*i).second] -= num;
    }

    // Optional members

    Real get_value(const Species& sp) const
    {
        SpeciesExpressionMatcher sexp(sp);
        Real retval(0);
        for (species_map_type::const_iterator i(index_map_.begin());
            i != index_map_.end(); ++i)
        {
            if (sexp.match((*i).first))
            {
                do
                {
                    retval += num_molecules_[(*i).second];
                } while (sexp.next());
            }
        }
        return retval;
    }

    Real get_value_exact(const Species& sp) const
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            // throw NotFound("Species not found");
            return 0.0;
        }

        return num_molecules_[(*i).second];
    }

    void set_value(const Species& sp, const Real& num)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            throw NotFound("Species not found");
        }

        num_molecules_[(*i).second] = num;
    }

    void save(const std::string& filename) const;
    void load(const std::string& filename);

    bool has_species(const Species& sp)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        return (i != index_map_.end());
    }

    void reserve_species(const Species& sp)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i != index_map_.end())
        {
            throw AlreadyExists("Species already exists");
        }

        index_map_.insert(std::make_pair(sp, num_molecules_.size()));
        species_.push_back(sp);
        num_molecules_.push_back(0);
    }

    void release_species(const Species& sp)
    {
        species_map_type::iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            throw NotFound("Species not found");
        }

        species_map_type::mapped_type
            idx((*i).second), last_idx(num_molecules_.size() - 1);
        if (idx != last_idx)
        {
            const species_container_type::size_type
                idx_(static_cast<species_container_type::size_type>(idx)),
                last_idx_(
                    static_cast<species_container_type::size_type>(last_idx));
            const Species& last_sp(species_[last_idx_]);
            species_[idx_] = last_sp;
            num_molecules_[idx] = num_molecules_[last_idx];
            index_map_[last_sp] = idx;
        }

        species_.pop_back();
        num_molecules_.pop_back();
        index_map_.erase(sp);
    }

    void bind_to(boost::shared_ptr<Model> model)
    {
        if (boost::shared_ptr<Model> bound_model = lock_model())
        {
            if (bound_model.get() != model.get())
            {
                std::cerr << "Warning: Model already bound to ODEWorld."
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

protected:

    Real3 edge_lengths_;
    Real volume_;
    Real t_;

    num_molecules_container_type num_molecules_;
    species_container_type species_;
    species_map_type index_map_;

    boost::weak_ptr<Model> model_;
    bool is_netfree_;
};

} // ode

} // ecell4

#endif /* __ECELL4_ODE_ODE_WORLD_HPP */
