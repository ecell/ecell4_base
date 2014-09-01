#ifndef __ECELL4_SUBVOLUME_SPACE_HPP
#define __ECELL4_SUBVOLUME_SPACE_HPP

#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "exceptions.hpp"
#include "Species.hpp"
#include "Space.hpp"
#include "Global.hpp"

namespace ecell4
{

class SubvolumeSpace
    : public Space
{
public:

    typedef Integer coordinate_type;

public:

    SubvolumeSpace()
        : t_(0.0)
    {
    }

    virtual ~SubvolumeSpace()
    {
        ;
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

    virtual Integer num_molecules(const Species& sp) const
    {
        return num_molecules_exact(sp);
    }

    virtual Integer num_molecules_exact(const Species& sp) const
    {
        throw NotImplemented("num_molecules_exact(const Species&) not implemented");
    }

    virtual Real get_value(const Species& sp) const
    {
        return static_cast<Real>(num_molecules(sp));
    }

    virtual Real get_value_exact(const Species& sp) const
    {
        return static_cast<Real>(num_molecules_exact(sp));
    }

    virtual const Integer num_subvolumes() const = 0;
    virtual const Real subvolume() const = 0;
    virtual coordinate_type global2coord(const Global& g) const = 0;
    virtual Global coord2global(const coordinate_type& c) const = 0;
    virtual Integer num_molecules(const Species& sp, const coordinate_type& c) const = 0;
    virtual Integer num_molecules_exact(const Species& sp, const coordinate_type& c) const = 0;
    virtual void add_molecules(const Species& sp, const Integer& num, const coordinate_type& c) = 0;
    virtual void remove_molecules(const Species& sp, const Integer& num, const coordinate_type& c) = 0;
    virtual std::vector<Species> list_species() const = 0;

protected:

    double t_;
};

class SubvolumeSpaceVectorImpl
    : public SubvolumeSpace
{
public:

    typedef SubvolumeSpace base_type;
    typedef base_type::coordinate_type coordinate_type;

    typedef std::vector<Integer> cell_type;
    typedef utils::get_mapper_mf<Species, cell_type>::type matrix_type;

public:

    SubvolumeSpaceVectorImpl(const Position3& edge_lengths,
        const Integer& cx, const Integer& cy, const Integer& cz)
        : base_type()
    {
        cell_sizes_[0] = cx;
        cell_sizes_[1] = cy;
        cell_sizes_[2] = cz;

        set_edge_lengths(edge_lengths);
    }

    virtual ~SubvolumeSpaceVectorImpl()
    {
        ;
    }

    const Position3& edge_lengths() const
    {
        return edge_lengths_;
    }

    void set_edge_lengths(const Position3& edge_lengths)
    {
        for (Position3::size_type dim(0); dim < 3; ++dim)
        {
            if (edge_lengths[dim] <= 0)
            {
                throw std::invalid_argument("the edge length must be positive.");
            }
        }

        edge_lengths_ = edge_lengths;
    }

    const Real volume() const
    {
        return edge_lengths_[0] * edge_lengths_[1] * edge_lengths_[2];
    }

    const Integer num_subvolumes() const
    {
        return cell_sizes_[0] * cell_sizes_[1] * cell_sizes_[2];
    }

    const Real subvolume() const
    {
        return volume() / num_subvolumes();
    }

    coordinate_type global2coord(const Global& g) const
    {
        return g.col + cell_sizes_[0] * (g.row + cell_sizes_[1] * g.layer);
    }

    Global coord2global(const coordinate_type& c) const
    {
        const Integer rowcol(cell_sizes_[0] * cell_sizes_[1]);
        const Integer layer(static_cast<Integer>(c / rowcol));
        const Integer surplus(c - layer * rowcol);
        const Integer row(surplus / cell_sizes_[0]);

        return Global(surplus - row * cell_sizes_[0], row, layer);
    }

    Integer num_molecules(const Species& sp) const;
    Integer num_molecules_exact(const Species& sp) const;
    Integer num_molecules(const Species& sp, const coordinate_type& c) const;
    Integer num_molecules_exact(const Species& sp, const coordinate_type& c) const;
    void add_molecules(const Species& sp, const Integer& num, const coordinate_type& c);
    void remove_molecules(const Species& sp, const Integer& num, const coordinate_type& c);

    std::vector<Species> list_species() const
    {
        std::vector<Species> retval;
        for (matrix_type::const_iterator i(matrix_.begin()); i != matrix_.end(); ++i)
        {
            retval.push_back((*i).first);
        }
        return retval;
    }

protected:

    void reserve_species(const Species& sp, const coordinate_type& c);

protected:

    Position3 edge_lengths_;
    boost::array<Integer, 3> cell_sizes_;
    matrix_type matrix_;
    std::vector<Species> species_;
};

} // ecell4

#endif /* __ECELL4_SUBVOLUME_SPACE_HPP */
