#ifndef ECELL4_PARTICLE_SPACE_CELL_LIST_IMPL_HPP
#define ECELL4_PARTICLE_SPACE_CELL_LIST_IMPL_HPP

#include <set>
#include <boost/multi_array.hpp>

#include "ParticleSpace.hpp"

#ifdef WITH_HDF5
#include "ParticleSpaceHDF5Writer.hpp"
#endif

#include "Integer3.hpp"


namespace ecell4
{

class ParticleSpaceCellListImpl
    : public ParticleSpace
{
public:

    typedef ParticleSpace base_type;
    typedef ParticleSpace::particle_container_type particle_container_type;

    typedef utils::get_mapper_mf<ParticleID, particle_container_type::size_type>::type
        key_to_value_map_type;

    typedef std::set<ParticleID> particle_id_set;
    typedef std::map<Species::serial_type, particle_id_set> per_species_particle_id_set;

    typedef std::vector<particle_container_type::size_type> cell_type; // sorted
    typedef boost::multi_array<cell_type, 3> matrix_type;
    typedef boost::array<matrix_type::size_type, 3> cell_index_type;
    typedef boost::array<matrix_type::difference_type, 3> cell_offset_type;

public:

    ParticleSpaceCellListImpl(const Real3& edge_lengths)
        : base_type(), edge_lengths_(edge_lengths), matrix_(boost::extents[3][3][3])
    {
        cell_sizes_[0] = edge_lengths_[0] / matrix_.shape()[0];
        cell_sizes_[1] = edge_lengths_[1] / matrix_.shape()[1];
        cell_sizes_[2] = edge_lengths_[2] / matrix_.shape()[2];
    }

    ParticleSpaceCellListImpl(
        const Real3& edge_lengths, const Integer3& matrix_sizes)
        : base_type(), edge_lengths_(edge_lengths),
        matrix_(boost::extents[matrix_sizes.col][matrix_sizes.row][matrix_sizes.layer])
    {
        cell_sizes_[0] = edge_lengths_[0] / matrix_.shape()[0];
        cell_sizes_[1] = edge_lengths_[1] / matrix_.shape()[1];
        cell_sizes_[2] = edge_lengths_[2] / matrix_.shape()[2];
    }

    void diagnosis() const
    {
        for (matrix_type::size_type i(0); i < matrix_.shape()[0]; ++i)
        {
            for (matrix_type::size_type j(0); j < matrix_.shape()[1]; ++j)
            {
                for (matrix_type::size_type k(0); k < matrix_.shape()[2]; ++k)
                {
                    const cell_type& c = matrix_[i][j][k];
                    for (cell_type::const_iterator it(c.begin()); it != c.end(); ++it)
                    {
                        if (*it >= particles_.size())
                        {
                            throw IllegalState("out of bounds.");
                        }
                    }
                }
            }
        }
    }

    // Space

    virtual Integer num_species() const
    {
        return particle_pool_.size();
    }
    virtual bool has_species(const Species& sp) const
    {
        return (particle_pool_.find(sp.serial()) != particle_pool_.end());
    }

    virtual std::vector<Species> list_species() const
    {
        std::vector<Species> retval;
        for (per_species_particle_id_set::const_iterator
            i(particle_pool_.begin()); i != particle_pool_.end(); ++i)
        {
            retval.push_back(Species((*i).first));
        }
        return retval;
    }

    // ParticleSpaceTraits

    const Real3& edge_lengths() const
    {
        return edge_lengths_;
    }

    const Real3& cell_sizes() const
    {
        return cell_sizes_;
    }

    const Integer3 matrix_sizes() const
    {
        return Integer3(matrix_.shape()[0], matrix_.shape()[1], matrix_.shape()[2]);
    }

    void reset(const Real3& edge_lengths);

    bool update_particle(const ParticleID& pid, const Particle& p);

    const particle_container_type& particles() const
    {
        return particles_;
    }

    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const;
    bool has_particle(const ParticleID& pid) const;
    void remove_particle(const ParticleID& pid);

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

    virtual void save(const std::string& filename) const
    {
        throw NotSupported(
            "save(const std::string) is not supported by this space class");
    }

#ifdef WITH_HDF5
    void save_hdf5(H5::Group* root) const
    {
        save_particle_space(*this, root);
    }

    void load_hdf5(const H5::Group& root)
    {
        load_particle_space(root, this);
    }
#endif

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

protected:

    // inline cell_index_type index(const Real3& pos, double t = 1e-10) const
    inline cell_index_type index(const Real3& pos) const
    {
        cell_index_type retval = {{
            static_cast<matrix_type::size_type>(
                pos[0] / cell_sizes_[0]) % matrix_.shape()[0],
            static_cast<matrix_type::size_type>(
                pos[1] / cell_sizes_[1]) % matrix_.shape()[1],
            static_cast<matrix_type::size_type>(
                pos[2] / cell_sizes_[2]) % matrix_.shape()[2]
            }}; // boost::array<matrix_type::size_type, 3>
        return retval;
    }

    inline Real3 offset_index_cyclic(
        cell_index_type& i, const cell_offset_type& o) const
    {
        Real3 retval;

        if (o[0] < 0 &&
            static_cast<matrix_type::size_type>(-o[0]) > i[0])
        {
            matrix_type::size_type t(
                (i[0] + matrix_.shape()[0] - (-o[0] % matrix_.shape()[0]))
                % matrix_.shape()[0]);
            retval[0] = (o[0] - static_cast<matrix_type::difference_type>(t - i[0]))
                * cell_sizes_[0];
            i[0] = t;
        }
        else if (matrix_.shape()[0] - o[0] <= i[0])
        {
            matrix_type::size_type
                t((i[0] + (o[0] % matrix_.shape()[0])) % matrix_.shape()[0]);
            retval[0] = (o[0] - static_cast<matrix_type::difference_type>(t - i[0]))
                * cell_sizes_[0];
            i[0] = t;
        }
        else
        {
            i[0] += o[0];
        }

        if (o[1] < 0 &&
            static_cast<matrix_type::size_type>(-o[1]) > i[1])
        {
            matrix_type::size_type t(
                (i[1] + matrix_.shape()[1] - (-o[1] % matrix_.shape()[1]))
                % matrix_.shape()[1]);
            retval[1] = (o[1] - static_cast<matrix_type::difference_type>(t - i[1]))
                * cell_sizes_[1];
            i[1] = t;
        }
        else if (matrix_.shape()[1] - o[1] <= i[1])
        {
            matrix_type::size_type
                t((i[1] + (o[1] % matrix_.shape()[1])) % matrix_.shape()[1]);
            retval[1] = (o[1] - static_cast<matrix_type::difference_type>(t - i[1]))
                * cell_sizes_[1];
            i[1] = t;
        }
        else
        {
            i[1] += o[1];
        }

        if (o[2] < 0 &&
            static_cast<matrix_type::size_type>(-o[2]) > i[2])
        {
            matrix_type::size_type
                t((i[2] + matrix_.shape()[2] - (-o[2] % matrix_.shape()[2]))
                % matrix_.shape()[2]);
            retval[2] = (o[2] - static_cast<matrix_type::difference_type>(t - i[2]))
                * cell_sizes_[2];
            i[2] = t;
        }
        else if (matrix_.shape()[2] - o[2] <= i[2])
        {
            matrix_type::size_type t(
                (i[2] + (o[2] % matrix_.shape()[2])) % matrix_.shape()[2]);
            retval[2] = (o[2] - static_cast<matrix_type::difference_type>(t - i[2]))
                * cell_sizes_[2];
            i[2] = t;
        }
        else
        {
            i[2] += o[2];
        }

        return retval;
    }

    inline const cell_type& cell(const cell_index_type& i) const
    {
        return matrix_[i[0]][i[1]][i[2]];
    }

    inline cell_type& cell(const cell_index_type& i)
    {
        return matrix_[i[0]][i[1]][i[2]];
    }

    inline particle_container_type::iterator find(const ParticleID& k)
    {
        key_to_value_map_type::const_iterator p(rmap_.find(k));
        if (rmap_.end() == p)
        {
            return particles_.end();
        }
        return particles_.begin() + (*p).second;
    }

    inline particle_container_type::const_iterator find(const ParticleID& k) const
    {
        key_to_value_map_type::const_iterator p(rmap_.find(k));
        if (rmap_.end() == p)
        {
            return particles_.end();
        }
        return particles_.begin() + (*p).second;
    }

    inline particle_container_type::iterator update(
        particle_container_type::iterator const& old_value,
        const std::pair<ParticleID, Particle>& v)
    {
        cell_type* new_cell(&cell(index(v.second.position())));
        cell_type* old_cell(0);

        if (old_value != particles_.end())
        {
            old_cell = &cell(index((*old_value).second.position()));
        }

        if (new_cell == old_cell)
        {
            // reinterpret_cast<nonconst_value_type&>(*old_value) = v;
            *old_value = v;
            return old_value;
        }
        else
        {
            particle_container_type::size_type idx(0);

            if (old_cell)
            {
                // reinterpret_cast<nonconst_value_type&>(*old_value) = v;
                *old_value = v;

                cell_type::iterator
                    i(find_in_cell(old_cell, old_value - particles_.begin()));
                idx = *i;
                erase_from_cell(old_cell, i);
                push_into_cell(new_cell, idx);
            }
            else
            {
                idx = particles_.size();
                particles_.push_back(v);
                push_into_cell(new_cell, idx);
                rmap_[v.first] = idx;
            }
            return particles_.begin() + idx;
        }
    }

    inline std::pair<particle_container_type::iterator, bool> update(
        const std::pair<ParticleID, Particle>& v)
    {
        cell_type* new_cell(&cell(index(v.second.position())));
        particle_container_type::iterator old_value(particles_.end());
        cell_type* old_cell(0);

        {
            key_to_value_map_type::const_iterator i(rmap_.find(v.first));
            if (i != rmap_.end())
            {
                old_value = particles_.begin() + (*i).second;
                old_cell = &cell(index(old_value->second.position()));
            }
        }

        if (new_cell == old_cell)
        {
            // reinterpret_cast<nonconst_value_type&>(*old_value) = v;
            *old_value = v;
            // return std::pair<particle_container_type::iterator, bool>(old_value, false);
            return std::make_pair(old_value, false);
        }
        else
        {
            particle_container_type::size_type idx(0);

            if (old_cell)
            {
                // reinterpret_cast<nonconst_value_type&>(*old_value) = v;
                *old_value = v;

                cell_type::iterator
                    i(find_in_cell(old_cell, old_value - particles_.begin()));
                idx = *i;
                erase_from_cell(old_cell, i);
                push_into_cell(new_cell, idx);
                return std::pair<particle_container_type::iterator, bool>(
                    particles_.begin() + idx, false);
            }
            else
            {
                idx = particles_.size();
                particles_.push_back(v);
                push_into_cell(new_cell, idx);
                rmap_[v.first] = idx;
                return std::pair<particle_container_type::iterator, bool>(
                    particles_.begin() + idx, true);
            }
        }
    }

    inline bool erase(particle_container_type::iterator const& i)
    {
        if (particles_.end() == i)
        {
            return false;
        }

        particle_container_type::size_type old_idx(i - particles_.begin());
        cell_type& old_cell(cell(index((*i).second.position())));
        const bool succeeded(erase_from_cell(&old_cell, old_idx));
        assert(succeeded);
        // BOOST_ASSERT(succeeded);
        rmap_.erase((*i).first);

        particle_container_type::size_type const last_idx(particles_.size() - 1);

        if (old_idx < last_idx)
        {
            const std::pair<ParticleID, Particle>& last(particles_[last_idx]);
            cell_type& last_cell(cell(index(last.second.position())));
            const bool tmp(erase_from_cell(&last_cell, last_idx));
            // BOOST_ASSERT(tmp);
            assert(succeeded);
            push_into_cell(&last_cell, old_idx);
            rmap_[last.first] = old_idx;
            // reinterpret_cast<nonconst_value_type&>(*i) = last;
            (*i) = last;
        }
        particles_.pop_back();
        return true;
    }

    inline bool erase(const ParticleID& k)
    {
        key_to_value_map_type::const_iterator p(rmap_.find(k));
        if (rmap_.end() == p)
        {
            return false;
        }
        return erase(particles_.begin() + (*p).second);
    }

    inline void erase_from_cell(cell_type* c, const cell_type::iterator& i)
    {
        c->erase(i);
    }

    inline cell_type::size_type erase_from_cell(
        cell_type* c, const particle_container_type::size_type& v)
    {
        cell_type::iterator e(c->end());
        std::pair<cell_type::iterator, cell_type::iterator>
            i(std::equal_range(c->begin(), e, v));
        const cell_type::size_type retval(i.second - i.first);
        c->erase(i.first, i.second);
        return retval;
    }

    inline void push_into_cell(
        cell_type* c, const particle_container_type::size_type& v)
    {
        cell_type::iterator i(std::upper_bound(c->begin(), c->end(), v));
        c->insert(i, v);
    }

    inline cell_type::iterator find_in_cell(
        cell_type* c, const particle_container_type::size_type& v)
    {
        cell_type::iterator i(std::lower_bound(c->begin(), c->end(), v));
        if (i != c->end() && *i == v)
        {
            return i;
        }
        else
        {
            return c->end();
        }
    }

    inline cell_type::const_iterator find_in_cell(
        cell_type* c, const particle_container_type::size_type& v) const
    {
        cell_type::iterator i(std::lower_bound(c->begin(), c->end(), v));
        if (i != c->end() && *i == v)
        {
            return i;
        }
        else
        {
            return c->end();
        }
    }

protected:

    Real3 edge_lengths_;

    particle_container_type particles_;
    key_to_value_map_type rmap_;
    per_species_particle_id_set particle_pool_;

    matrix_type matrix_;
    Real3 cell_sizes_;
};

}; // ecell4

#endif /* ECELL4_PARTICLE_SPACE_CELL_LIST_IMPL_HPP */
