#ifndef __ECELL4_SUBVOLUME_SPACE_HPP
#define __ECELL4_SUBVOLUME_SPACE_HPP

#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "exceptions.hpp"
#include "Species.hpp"
#include "Space.hpp"
#include "Global.hpp"
#include "SubvolumeSpaceHDF5Writer.hpp"

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

    virtual const Global matrix_sizes() const = 0;
    virtual const Position3 subvolume_edge_lengths() const = 0;
    virtual const Integer num_subvolumes() const = 0;
    virtual const Real subvolume() const = 0;
    virtual coordinate_type global2coord(const Global& g) const = 0;
    virtual Global coord2global(const coordinate_type& c) const = 0;
    virtual Integer num_molecules(
        const Species& sp, const coordinate_type& c) const = 0;
    virtual Integer num_molecules_exact(
        const Species& sp, const coordinate_type& c) const = 0;
    virtual void add_molecules(
        const Species& sp, const Integer& num, const coordinate_type& c) = 0;
    virtual void remove_molecules(
        const Species& sp, const Integer& num, const coordinate_type& c) = 0;
    virtual const std::vector<Species>& species() const = 0;
    virtual std::vector<Species> list_species() const = 0;
    virtual coordinate_type get_neighbor(
        const coordinate_type& c, const Integer rnd) const = 0;

    virtual Integer num_molecules(const Species& sp, const Global& g) const
    {
        return num_molecules(sp, global2coord(g));
    }

    virtual Integer num_molecules_exact(const Species& sp, const Global& g) const
    {
        return num_molecules_exact(sp, global2coord(g));
    }

    virtual void add_molecules(const Species& sp, const Integer& num, const Global& g)
    {
        add_molecules(sp, num, global2coord(g));
    }

    virtual void remove_molecules(const Species& sp, const Integer& num, const Global& g)
    {
        remove_molecules(sp, num, global2coord(g));
    }

    virtual void cleanup(const Position3& edge_lengths, const Global& matrix_sizes) = 0;
    virtual void save(H5::Group* root) const = 0;
    virtual void load(const H5::Group& root) = 0;

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

    SubvolumeSpaceVectorImpl(
        const Position3& edge_lengths, const Global matrix_sizes)
        : base_type()
    {
        matrix_sizes_[0] = matrix_sizes.col;
        matrix_sizes_[1] = matrix_sizes.row;
        matrix_sizes_[2] = matrix_sizes.layer;

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

    const Global matrix_sizes() const
    {
        return Global(matrix_sizes_[0], matrix_sizes_[1], matrix_sizes_[2]);
    }

    const Position3 subvolume_edge_lengths() const
    {
        return Position3(
            edge_lengths_[0] / matrix_sizes_[0],
            edge_lengths_[1] / matrix_sizes_[1],
            edge_lengths_[2] / matrix_sizes_[2]);
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
        return matrix_sizes_[0] * matrix_sizes_[1] * matrix_sizes_[2];
    }

    const Real subvolume() const
    {
        return volume() / num_subvolumes();
    }

    coordinate_type global2coord(const Global& g) const
    {
        const coordinate_type coord(
            modulo(g.col, matrix_sizes_[0])
                + matrix_sizes_[0] * (modulo(g.row, matrix_sizes_[1])
                    + matrix_sizes_[1] * modulo(g.layer, matrix_sizes_[2])));
        return coord;
    }

    Global coord2global(const coordinate_type& c) const
    {
        const Integer rowcol(matrix_sizes_[0] * matrix_sizes_[1]);
        const Integer layer(static_cast<Integer>(c / rowcol));
        const Integer surplus(c - layer * rowcol);
        const Integer row(static_cast<Integer>(surplus / matrix_sizes_[0]));
        return Global(surplus - row * matrix_sizes_[0], row, layer);
    }

    Integer num_molecules(const Species& sp) const;
    Integer num_molecules_exact(const Species& sp) const;

    Integer num_molecules(const Species& sp, const coordinate_type& c) const;
    Integer num_molecules_exact(const Species& sp, const coordinate_type& c) const;
    void add_molecules(const Species& sp, const Integer& num, const coordinate_type& c);
    void remove_molecules(const Species& sp, const Integer& num, const coordinate_type& c);

    coordinate_type get_neighbor(const coordinate_type& c, const Integer rnd) const;

    const std::vector<Species>& species() const
    {
        return species_;
    }

    std::vector<Species> list_species() const
    {
        // std::vector<Species> retval;
        // for (matrix_type::const_iterator i(matrix_.begin()); i != matrix_.end(); ++i)
        // {
        //     retval.push_back((*i).first);
        // }
        // return retval;
        return species_;
    }

    void save(H5::Group* root) const
    {
        save_subvolume_space(*this, root);
    }

    void load(const H5::Group& root)
    {
        load_subvolume_space(root, this);
    }

    void cleanup(const Position3& edge_lengths, const Global& matrix_sizes)
    {
        set_edge_lengths(edge_lengths);
        matrix_sizes_[0] = matrix_sizes.col;
        matrix_sizes_[1] = matrix_sizes.row;
        matrix_sizes_[2] = matrix_sizes.layer;
        matrix_.clear();
        species_.clear();
    }

protected:

    void reserve_species(const Species& sp, const coordinate_type& c);

protected:

    Position3 edge_lengths_;
    boost::array<Integer, 3> matrix_sizes_;
    matrix_type matrix_;
    std::vector<Species> species_;
};

} // ecell4

#endif /* __ECELL4_SUBVOLUME_SPACE_HPP */
