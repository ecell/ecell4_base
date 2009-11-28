#ifndef TWOFOLD_CONTAINER_HPP
#define TWOFOLD_CONTAINER_HPP

#include <boost/array.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_categories.hpp>
#include "utils.hpp"

template<typename T_>
class twofold_container
{
public:
    typedef T_ value_type;
private:
    typedef boost::array<value_type, 2> containing_type;
public:
    typedef typename containing_type::reference reference;
    typedef typename containing_type::const_reference const_reference;
    typedef typename containing_type::size_type size_type;
    typedef typename containing_type::difference_type difference_type;

    class const_iterator;
    class iterator
        : public boost::iterator_facade<
            iterator, value_type, boost::forward_traversal_tag>
    {
        friend class const_iterator;
        friend class boost::iterator_core_access;

        std::ptrdiff_t distance_to(iterator const& that) const
        {
            return that.idx_ - idx_;
        }

        bool equal(iterator const& that) const
        {
            return &cntnr_ == &that.cntnr_ && idx_ == that.idx_;
        }

        void increment()
        {
            ++idx_;
        }

        value_type& dereference() const
        {
            return cntnr_[idx_];
        }

    public:
        iterator(twofold_container& cntnr, size_type idx)
            : cntnr_(cntnr), idx_(idx) {}

        iterator(const_iterator const&);

    private:
        twofold_container& cntnr_;
        size_type idx_;
    };

    class const_iterator
        : public boost::iterator_facade<
            const_iterator, const value_type, boost::forward_traversal_tag>
    {
        friend class iterator;
        friend class boost::iterator_core_access;

        std::ptrdiff_t distance_to(const_iterator const& that) const
        {
            return that.idx_ - idx_;
        }

        bool equal(const_iterator const& that) const
        {
            return &cntnr_ == &that.cntnr_ && idx_ == that.idx_;
        }

        void increment()
        {
            ++idx_;
        }

        value_type const& dereference() const
        {
            return cntnr_[idx_];
        }

    public:
        const_iterator(twofold_container const& cntnr, size_type idx)
            : cntnr_(cntnr), idx_(idx) {}

        const_iterator(iterator const& that)
            : cntnr_(that.cntnr_), idx_(that.idx_) {}

    private:
        twofold_container const& cntnr_;
        size_type idx_;
    };

public:
    twofold_container()
    {
        items_[0] = value_type();
        items_[1] = value_type();
    }

    twofold_container(value_type const& one)
    {
        BOOST_ASSERT(one);
        items_[0] = one;
        items_[1] = value_type();
    }

    twofold_container(value_type const& one, value_type const& two)
    {
        BOOST_ASSERT(one);
        BOOST_ASSERT(two);
        if (one <= two)
        {
            items_[0] = one;
            items_[1] = two;
        }
        else
        {
            items_[0] = two;
            items_[1] = one;
        }
    }

    size_type size() const
    {
        return items_[0] ? items_[1] ? 2: 1: 0;
    }

    iterator begin()
    {
        return iterator(*this, 0);
    }

    iterator end()
    {
        return iterator(*this, size());
    }

    const_iterator begin() const
    {
        return const_iterator(*this, 0);
    }

    const_iterator end() const
    {
        return const_iterator(*this, size());
    }

    void push_back(value_type const& item)
    {
        if (!items_[0])
        {
            items_[0] = item;
        }
        else if (!items_[1])
        {
            items_[1] = item;
        }
        else
        {
            BOOST_ASSERT(false);
        }
    }

    value_type& operator[](std::size_t idx)
    {
        return items_[idx];
    }

    value_type const& operator[](std::size_t idx) const
    {
        return items_[idx];
    }

    bool operator<(twofold_container const& rhs) const
    {
        return memberwise_compare(*this, rhs) < 0;
    }

    bool operator>=(twofold_container const& rhs) const
    {
        return !operator<(rhs);
    }

    bool operator>(twofold_container const& rhs) const
    {
        return memberwise_compare(*this, rhs) > 0;
    }

    bool operator<=(twofold_container const& rhs) const
    {
        return !operator>(rhs);
    }

    bool operator==(twofold_container const& rhs) const
    {
        if (rhs.size() != size())
            return false;
        switch (size())
        {
        case 0:
            return true;
        case 1:
            return items_[0] == rhs[0];
        case 2:
            return items_[0] == rhs[0] && items_[1] == rhs[1];
        }
        /* never get here */
        return false;
    }

    bool operator!=(twofold_container const& rhs) const
    {
        return !operator==(rhs);
    }

protected:
    containing_type items_;
};

template<typename T_>
inline twofold_container<T_>::iterator::iterator(
        typename twofold_container<T_>::const_iterator const& that)
    : cntnr_(const_cast<twofold_container&>(that.cntnr_)), idx_(that.idx_)
{
}

#endif /* TWOFOLD_CONTAINER_HPP */
