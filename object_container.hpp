#ifndef OBJECT_CONTAINER_HPP
#define OBJECT_CONTAINER_HPP

#include <cstddef>
#include <algorithm>
#include <iterator>
#include <boost/multi_array.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/add_const.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/list/at.hpp>
#include "vector3.hpp"
#include "array_helper.hpp"
#include "utils.hpp"

template<typename Tobj_, typename Tkey_,
        template<typename, typename> class MFget_mapper_ =
            get_default_impl::std::template map>
class object_container
{
public:
    typedef typename Tobj_::length_type length_type;
    typedef Tkey_ key_type;
    typedef Tobj_ mapped_type;
    typedef std::pair<const key_type, mapped_type> value_type;
    typedef vector3<length_type> position_type;
    typedef typename MFget_mapper_<Tkey_, mapped_type>::type cell_type;
    typedef boost::multi_array<cell_type, 3> matrix_type;
    typedef typename cell_type::size_type size_type;
    typedef boost::array<typename matrix_type::size_type, 3>
            cell_index_type;
    typedef boost::array<typename matrix_type::difference_type, 3>
            cell_offset_type;
    typedef typename cell_type::reference reference;
    typedef typename cell_type::const_reference const_reference;
    typedef typename MFget_mapper_<key_type, cell_type*>::type
            id_to_cell_mapper_type;
    typedef id_to_cell_mapper_type key_to_cell_mapper_type;

    template<typename Thost_, typename Treftype_, typename Tciter_>
    class iterator_base
    {
    private:
        struct reference_is_const:
            public boost::is_const<
                    typename boost::remove_reference<Treftype_>::type> {};
    public:
        typedef typename object_container::value_type value_type;
        typedef Treftype_ reference;
        typedef const Treftype_ const_reference;
        typedef std::bidirectional_iterator_tag iterator_category;
        typedef typename Tciter_::difference_type difference_type;
        typedef typename boost::mpl::if_<reference_is_const,
                typename boost::add_const<value_type>::type,
                value_type>::type* pointer;

    protected:
        typedef typename boost::mpl::if_<reference_is_const,
                typename boost::add_const<cell_type>::type,
                cell_type>::type* cell_ptr_type;
        typedef std::pair<cell_ptr_type, cell_ptr_type> cell_range_type;

    public:
        bool operator==(const Thost_& rhs) const
        {
            return cell_p_ == rhs.cell_p_ &&
                cell_iter_ == rhs.cell_iter_;
        }

        bool operator!=(const Thost_& rhs) const
        {
            return !operator==(rhs);
        }

        const_reference operator*() const
        {
            return *cell_iter_;
        }

        reference operator*()
        {
            return *cell_iter_;
        }

        Thost_& operator++()
        {
            if (cell_iter_ != cell_p_->end())
            {
                ++cell_iter_;
            }

            while (cell_iter_ == cell_p_->end())
            {
                ++cell_p_;
                if (cell_p_ == cell_r_.second)
                {
                    --cell_p_;
                    cell_iter_ = cell_p_->end();
                    break;
                }
                cell_iter_ = cell_p_->begin();
            }
            return *static_cast<Thost_*>(this);
        }

        Thost_ operator++(int)
        {
            Thost_ retval(static_cast<Thost_ const&>(*this));
            operator++();
            return retval;
        }

        Thost_& operator--()
        {
            while (cell_iter_ == cell_p_->begin() && cell_p_ > cell_r_.first)
            {
                --cell_p_;
                cell_iter_ = cell_p_->end();
            }
            --cell_iter_;
            return *static_cast<Thost_*>(this);
        }

        Thost_ operator--(int)
        {
            Thost_ retval(static_cast<Thost_ const&>(*this));
            operator--();
            return retval;
        }

        Thost_& operator=(const Thost_& rhs)
        {
            cell_p_ = rhs.cell_p_;
            cell_iter_ = rhs.cell_iter_;
            return *static_cast<Thost_*>(this);
        }

        iterator_base(cell_ptr_type cell_p, cell_range_type const& cell_r,
                const Tciter_& cell_iter)
                : cell_p_(cell_p), cell_r_(cell_r), cell_iter_(cell_iter) {}

    public:
        cell_ptr_type cell_p_;
        cell_range_type cell_r_;
        Tciter_ cell_iter_;
    };

    class iterator: public iterator_base<iterator, reference,
            typename cell_type::iterator>
    {
    public:
        typedef iterator_base<iterator, reference,
                typename cell_type::iterator> base_type;

    public:
        iterator(typename base_type::cell_ptr_type cell_p,
                typename base_type::cell_range_type cell_r,
                const typename cell_type::iterator& cell_iter)
                : base_type(cell_p, cell_r, cell_iter) {}
    };

    class const_iterator: public iterator_base<const_iterator,
            const_reference, typename cell_type::const_iterator>
    {
    public:
        typedef iterator_base<const_iterator, const_reference,
                typename cell_type::const_iterator> base_type;

    public:
        const_iterator(typename base_type::cell_ptr_type cell_p,
                typename base_type::cell_range_type const& cell_r,
                const typename cell_type::const_iterator& cell_iter)
                : base_type(cell_p, cell_r, cell_iter) {}

        const_iterator(iterator const& that)
                : base_type(that.cell_p_, that.cell_r_, that.cell_iter_) {}
    };

public:
    object_container(length_type world_size = 1.0,
            typename matrix_type::size_type size = 1)
        : world_size_(world_size),
          cell_size_(world_size / size),
          matrix_(boost::extents[size][size][size]),
          size_(0)
    {
    }

    inline cell_index_type index(const position_type& pos,
            double t = 1e-10) const
    {
        return array_gen<typename matrix_type::size_type>(
            static_cast<typename matrix_type::size_type>(
                pos[0] / cell_size_ ) % matrix_.shape()[0],
            static_cast<typename matrix_type::size_type>(
                pos[1] / cell_size_ ) % matrix_.shape()[1],
            static_cast<typename matrix_type::size_type>(
                pos[2] / cell_size_ ) % matrix_.shape()[2] );
    }

    inline cell_offset_type offset(const position_type& pos,
            double t = 1e-10) const
    {
        return array_gen<typename matrix_type::difference_type>(
            pos[0] / cell_size, pos[1] / cell_size, pos[2] / cell_size);
    }

    inline bool offset_index(
            cell_index_type& i,
            const cell_offset_type& o) const
    {
        if ((o[0] < 0 && static_cast<size_type>(-o[0]) > i[0])
                || (matrix_.shape()[0] - o[0] <= i[0])
                || (o[1] < 0 && static_cast<size_type>(-o[1]) > i[1])
                || (matrix_.shape()[1] - o[1] <= i[1])
                || (o[2] < 0 && static_cast<size_type>(-o[2]) > i[2])
                || (matrix_.shape()[2] - o[2] <= i[2]))
        {
            return false;
        }
        i[0] += o[0];
        i[1] += o[1];
        i[2] += o[2];
        return true;
    }

    inline position_type offset_index_cyclic(cell_index_type& i,
                                             const cell_offset_type& o)
    {
        position_type retval;

        if (o[0] < 0 &&
            static_cast<typename matrix_type::size_type>(-o[0]) > i[0])
        {
            typename matrix_type::size_type t(
                (i[0] + matrix_.shape()[0] - (-o[0] % matrix_.shape()[0])) %
                matrix_.shape()[0]);
            retval[0] 
                = (o[0] - 
                   static_cast<typename matrix_type::difference_type>
                   (t - i[0])) * cell_size_;
            i[0] = t;
        }
        else if (matrix_.shape()[0] - o[0] <= i[0])
        {
            typename matrix_type::size_type t(
                    (i[0] + (o[0] % matrix_.shape()[0])) % matrix_.shape()[0]);
            retval[0] 
                = (o[0] - 
                   static_cast<typename matrix_type::difference_type>
                   (t - i[0])) * cell_size_;
            i[0] = t;
        }
        else
        {
            i[0] += o[0];
        }

        if (o[1] < 0 &&
                static_cast<typename matrix_type::size_type>(-o[1]) > i[1])
        {
            typename matrix_type::size_type t(
                    (i[1] + matrix_.shape()[1] - (-o[1] % matrix_.shape()[1])) %
                        matrix_.shape()[1]);
            retval[1] = (o[1] - static_cast<typename matrix_type::difference_type>(t - i[1])) * cell_size_;
            i[1] = t;
        }
        else if (matrix_.shape()[1] - o[1] <= i[1])
        {
            typename matrix_type::size_type t(
                    (i[1] + (o[1] % matrix_.shape()[1])) % matrix_.shape()[1]);
            retval[1] = (o[1] - static_cast<typename matrix_type::difference_type>(t - i[1])) * cell_size_;
            i[1] = t;
        }
        else
        {
            i[1] += o[1];
        }

        if (o[2] < 0 &&
                static_cast<typename matrix_type::size_type>(-o[2]) > i[2])
        {
            typename matrix_type::size_type t(
                    (i[2] + matrix_.shape()[2] - (-o[2] % matrix_.shape()[2])) %
                        matrix_.shape()[2]);
            retval[2] = (o[2] - static_cast<typename matrix_type::difference_type>(t - i[2])) * cell_size_;
            i[2] = t;
        }
        else if (matrix_.shape()[2] - o[2] <= i[2])
        {
            typename matrix_type::size_type t(
                    (i[2] + (o[2] % matrix_.shape()[2])) % matrix_.shape()[2]);
            retval[2] = (o[2] - static_cast<typename matrix_type::difference_type>(t - i[2])) * cell_size_;
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

    inline length_type world_size() const
    {
        return world_size_;
    }

    inline length_type cell_size() const
    {
        return cell_size_;
    }

    inline typename matrix_type::size_type matrix_size() const
    {
        return matrix_.shape()[0];
    }

    inline size_type size() const
    {
        return size_;
    }

    inline std::pair<iterator, bool> update(const value_type& v)
    {
        cell_type& c(cell(index(v.second.position())));
        typename key_to_cell_mapper_type::iterator kci(rmap_.find(v.first));
        std::pair<typename cell_type::iterator, bool> ir;
        if (rmap_.end() != kci)
        {
            if (&c != (*kci).second)
            {
                (*kci).second->erase(v.first);
                rmap_.erase(v.first);
                rmap_.insert(std::make_pair(v.first, &c));
                ir = c.insert(v);
            }
            else
            {
                ir.first = c.find(v.first);
                BOOST_ASSERT(c.end() != ir.first);
                if (v.first == (*ir.first).first) {
                    ir.first->second = v.second;
                } else {
                    c.erase(ir.first);
                    ir.first = c.insert(v).first;
                }
            }
            return std::make_pair(iterator(&c, cell_range(), ir.first), false);
        }
        else
        {
            ir = c.insert(v);
            rmap_.insert(std::make_pair(v.first, &c));
            ++size_;
            return std::make_pair(iterator(&c, cell_range(), ir.first), true);
        }
    }


    inline bool erase(const key_type& k)
    {
        typename key_to_cell_mapper_type::iterator p(rmap_.find(k));
        if (rmap_.end() == p)
        {
            return false;
        }
        (*p).second->erase(k);
        rmap_.erase(p);
        --size_;
        return true;
    }

    inline iterator begin()
    {
        BOOST_ASSERT(matrix_.num_elements() > 0);
        std::pair<cell_type*, cell_type*> r(cell_range());
        cell_type* cell_p(r.first);

        while (cell_p->begin() == cell_p->end())
        {
            ++cell_p;
            if (cell_p == r.second)
            {
                --cell_p;
                return iterator(cell_p, r, cell_p->end());
            }
        }

        return iterator(cell_p, r, cell_p->begin());
    }

    inline const_iterator begin() const
    {
        BOOST_ASSERT(matrix_.num_elements() > 0);
        std::pair<cell_type const*, cell_type const*> r(cell_range());
        cell_type const* cell_p(r.first);

        while (cell_p->begin() == cell_p->end())
        {
            ++cell_p;
            if (cell_p == r.second)
            {
                --cell_p;
                return const_iterator(cell_p, r, cell_p->end());
            }
        }

        return const_iterator(cell_p, r, cell_p->begin());
    }

    inline iterator end()
    {
        BOOST_ASSERT(matrix_.num_elements() > 0);
        std::pair<cell_type*, cell_type*> r(cell_range());
        return iterator(r.second - 1, r, (r.second - 1)->end());
    }

    inline const_iterator end() const
    {
        BOOST_ASSERT(matrix_.num_elements() > 0);
        std::pair<cell_type const*, cell_type const*> r(cell_range());
        return const_iterator(r.second - 1, r, (r.second - 1)->end());
    }

    inline iterator find(const key_type& k)
    {
        typename key_to_cell_mapper_type::iterator p(rmap_.find(k));
        if (rmap_.end() == p)
        {
            return end();
        }
        typename cell_type::iterator i((*p).second->find(k));
        if ((*p).second->end() == i)
        {
            return end();
        }
        return iterator((*p).second, cell_range(), i);
    }

    inline const_iterator find(const key_type& k) const
    {
        typename key_to_cell_mapper_type::const_iterator p(rmap_.find(k));
        if (rmap_.end() == p)
        {
            return end();
        }
        typename cell_type::const_iterator i((*p).second->find(k));
        if ((*p).second->end() == i)
        {
            return end();
        }
        return const_iterator((*p).second, cell_range(), i);
    }

    template<typename Tcollect_>
    inline void each_neighbor(const cell_index_type& idx, Tcollect_& collector)
    {
        cell_offset_type _off;
        each_neighbor_loops<Tcollect_>(3, _off, idx, collector);
    }

    template<typename Tcollect_>
    inline void each_neighbor(const cell_index_type& idx, Tcollect_ const& collector)
    {
        cell_offset_type _off;
        each_neighbor_loops<Tcollect_ const>(3, _off, idx, collector);
    }

    template<typename Tcollect_>
    inline void each_neighbor(const cell_index_type& idx, Tcollect_& collector) const
    {
        cell_offset_type _off;
        each_neighbor_loops<Tcollect_>(3, _off, idx, collector);
    }

    template<typename Tcollect_>
    inline void each_neighbor(const cell_index_type& idx, Tcollect_ const& collector) const
    {
        cell_offset_type _off;
        each_neighbor_loops<Tcollect_ const>(3, _off, idx, collector);
    }

    template<typename Tcollect_>
    inline void each_neighbor_cyclic(const cell_index_type& idx,
            Tcollect_& collector)
    {
        cell_offset_type _off;
        each_neighbor_cyclic_loops<Tcollect_>(3, _off, idx, collector);
    }

    template<typename Tcollect_>
    inline void each_neighbor_cyclic(const cell_index_type& idx,
            Tcollect_ const& collector)
    {
        cell_offset_type _off;
        each_neighbor_cyclic_loops<Tcollect_ const>(3, _off, idx, collector);
    }

    template<typename Tcollect_>
    inline void each_neighbor_cyclic(const cell_index_type& idx,
            Tcollect_& collector) const
    {
        cell_offset_type _off;
        each_neighbor_cyclic_loops<Tcollect_>(3, _off, idx, collector);
    }

    template<typename Tcollect_>
    inline void each_neighbor_cyclic(const cell_index_type& idx,
            Tcollect_ const& collector) const
    {
        cell_offset_type _off;
        each_neighbor_cyclic_loops<Tcollect_>(3, _off, idx, collector);
    }

private:
    std::pair<cell_type*, cell_type*> cell_range()
    {
        return std::make_pair(
            matrix_.origin(), matrix_.origin() + matrix_.num_elements());
    }

    std::pair<cell_type const*, cell_type const*> cell_range() const
    {
        return std::make_pair(
            matrix_.origin(), matrix_.origin() + matrix_.num_elements());
    }

    template<typename Tcollect_>
    inline void each_neighbor_loops(const std::size_t depth,
            cell_offset_type& off, const cell_index_type& idx,
            Tcollect_& collector) const
    {
        if (depth > 0)
        {
            for (off[depth - 1] = -1; off[depth - 1] <= 1; ++off[depth - 1])
            {
                each_neighbor_loops(depth - 1, off, idx, collector);
            }
        }
        else
        {
            cell_index_type _idx(idx);
            if (!offset_index(_idx, off)) {
                return;
            }
            cell_type& c(cell(_idx));
            for (typename cell_type::const_iterator i(c.begin()); i != c.end(); ++i) 
            {
                collector(const_iterator(&c, cell_range(), i));
            }
        }
    }

    template<typename Tcollect_>
    inline void each_neighbor_loops(const std::size_t depth,
            cell_offset_type& off, const cell_index_type& idx,
            Tcollect_& collector)
    {
        if (depth > 0)
        {
            for (off[depth - 1] = -1; off[depth - 1] <= 1; ++off[depth - 1])
            {
                each_neighbor_loops(depth - 1, off, idx, collector);
            }
        }
        else
        {
            cell_index_type _idx(idx);
            if (!offset_index(_idx, off)) {
                return;
            }
            cell_type& c(cell(_idx));
            for (typename cell_type::iterator i(c.begin()); i != c.end(); ++i) 
            {
                collector(iterator(&c, cell_range(), i));
            }
        }
    }

    template<typename Tcollect_>
    inline void each_neighbor_cyclic_loops(const std::size_t depth,
            cell_offset_type& off, const cell_index_type& idx,
            Tcollect_& collector) const
    {
        if (depth > 0)
        {
            for (off[depth - 1] = -1; off[depth - 1] <= 1; ++off[depth - 1])
            {
                each_neighbor_cyclic_loops(depth - 1, off, idx, collector);
            }
        }
        else
        {
            cell_index_type _idx(idx);
            const position_type pos_off(offset_index_cyclic(_idx, off));
            cell_type& c(cell(_idx));
            for (typename cell_type::const_iterator i(c.begin()); i != c.end(); ++i) 
            {
                collector(const_iterator(&c, cell_range(), i), pos_off);
            }
        }
    }

    template<typename Tcollect_>
    inline void each_neighbor_cyclic_loops(const std::size_t depth,
            cell_offset_type& off, const cell_index_type& idx,
            Tcollect_& collector)
    {
        if (depth > 0)
        {
            for (off[depth - 1] = -1; off[depth - 1] <= 1; ++off[depth - 1])
            {
                each_neighbor_cyclic_loops(depth - 1, off, idx, collector);
            }
        }
        else
        {
            cell_index_type _idx(idx);
            const position_type pos_off(offset_index_cyclic(_idx, off));
            cell_type& c(cell(_idx));
            for (typename cell_type::iterator i(c.begin()); i != c.end(); ++i) 
            {
                collector(iterator(&c, cell_range(), i), pos_off);
            }
        }
    }

private:
    const length_type world_size_;
    const length_type cell_size_;
    matrix_type matrix_;
    id_to_cell_mapper_type rmap_;
    size_type size_;
};

template<typename T_, typename Tkey_,
        template<typename, typename> class MFget_mapper_>
static inline typename object_container<T_, Tkey_, MFget_mapper_>::cell_index_type&
operator+=(
       typename object_container<T_,
                Tkey_, MFget_mapper_>::cell_index_type& lhs,
       const typename object_container<T_,
                Tkey_, MFget_mapper_>::cell_offset_type& rhs)
{
    rhs[0] += lhs[0];
    rhs[1] += lhs[1];
    rhs[2] += lhs[2];
    return rhs;
}


#endif /* OBJECT_CONTAINER_HPP */
