#ifndef REACTION_RULE_HPP
#define REACTION_RULE_HPP

#include <set>
#include <ostream>
#include <algorithm>

#include <boost/array.hpp>
#include <boost/bind.hpp>
#include <boost/range/size.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_categories.hpp>

#include "Defs.hpp"
#include "utils.hpp"
#include "SpeciesType.hpp"

class ReactionRule
{
private:
    typedef std::vector<SpeciesType const*> species_type_vector;

public:
    class Reactants
    {
    private:
        typedef boost::array<SpeciesType const*, 2> containing_type;
    public:
        typedef containing_type::value_type value_type;
        typedef containing_type::reference reference;
        typedef containing_type::const_reference const_reference;
        typedef containing_type::size_type size_type;
        typedef containing_type::difference_type difference_type;
 
        class const_iterator;
        class iterator
            : public boost::iterator_facade<
                iterator, value_type, boost::forward_traversal_tag>
        {
            friend class const_iterator;
            friend class boost::iterator_core_access;

            bool equal(iterator const& that) const
            {
                return &cntnr_ == &that.cntnr_ && idx_ == that.idx_;
            }

            void increment()
            {
                ++idx_;
            }

            SpeciesType const*& dereference() const
            {
                return cntnr_[idx_];
            }

        public:
            iterator(Reactants& cntnr, size_type idx)
                : cntnr_(cntnr), idx_(idx) {}

            iterator(const_iterator const&);

        private:
            Reactants& cntnr_;
            size_type idx_;
        };
 
        class const_iterator
            : public boost::iterator_facade<
                const_iterator, const value_type, boost::forward_traversal_tag>
        {
            friend class iterator;
            friend class boost::iterator_core_access;

            bool equal(const_iterator const& that) const
            {
                return &cntnr_ == &that.cntnr_ && idx_ == that.idx_;
            }

            void increment()
            {
                ++idx_;
            }

            SpeciesType const* const& dereference() const
            {
                return cntnr_[idx_];
            }

        public:
            const_iterator(Reactants const& cntnr, size_type idx)
                : cntnr_(cntnr), idx_(idx) {}

            const_iterator(iterator const& that)
                : cntnr_(that.cntnr_), idx_(that.idx_) {}

        private:
            Reactants const& cntnr_;
            size_type idx_;
        };

    public:
        Reactants()
        {
            items_[0] = 0;
            items_[1] = 0;
        }

        Reactants(SpeciesType const* one)
        {
            BOOST_ASSERT(one);
            items_[0] = one;
            items_[1] = 0;
        }

        Reactants(SpeciesType const* one, SpeciesType const* two)
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

        SpeciesType const*& operator[](std::size_t idx)
        {
            return items_[idx];
        }

        SpeciesType const* const& operator[](std::size_t idx) const
        {
            return items_[idx];
        }

        bool operator==(Reactants const& rhs) const
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

        bool operator!=(Reactants const& rhs) const
        {
            return !operator==(rhs);
        }

    protected:
        containing_type items_;
    };

    typedef species_type_vector::const_iterator species_type_iterator;
    typedef boost::iterator_range<species_type_iterator> species_type_range;

public:
    Reactants const& get_reactants() const
    {
        return reactants_;
    }

    void add_product(SpeciesType const* s)
    {
        products_.insert(
            std::lower_bound(products_.begin(), products_.end(), s),
            s);
    }

    species_type_range get_products() const
    {
        return products_;
    }

    Real k() const
    {
        return k_;
    }

    Real& k()
    {
        return k_;
    }

    explicit ReactionRule(Reactants const& _reactants, Real k = .0)
        : reactants_(_reactants), k_(k) {}

    template<typename Trange_>
    ReactionRule(Reactants const& _reactants, Trange_ const& products, Real k = .0)
        : reactants_(_reactants), k_(k)
    {
        std::for_each(boost::begin(products), boost::end(products),
                boost::bind(&ReactionRule::add_product, this, _1));
        std::stable_sort(products_.begin(), products_.end());
    }

private:
    Reactants reactants_;
    species_type_vector products_;
    Real k_;
};

inline ReactionRule::Reactants::iterator::iterator(const_iterator const& that)
    : cntnr_(const_cast<Reactants&>(that.cntnr_)), idx_(that.idx_)
{
}

inline bool operator<(ReactionRule::Reactants const& lhs, ReactionRule::Reactants const& rhs)
{
    return memberwise_compare(lhs, rhs) < 0;
}

inline bool operator<(ReactionRule const& lhs, ReactionRule const& rhs)
{
    int tmp = memberwise_compare(lhs.get_reactants(), rhs.get_reactants());
    if (tmp > 0)
    {
        return false;
    }
    else if (tmp < 0)
    {
        return true;
    }
    return memberwise_compare(lhs.get_products(), rhs.get_products()) < 0;
}

inline bool operator==(ReactionRule const& lhs, ReactionRule const& rhs)
{
    return lhs.get_reactants() == rhs.get_reactants() &&
            memberwise_compare(lhs.get_products(), rhs.get_products()) == 0;
}

inline bool operator!=(ReactionRule const& lhs, ReactionRule const& rhs)
{
    return !(lhs == rhs);
}

template<typename T2_>
inline ReactionRule new_reaction_rule(SpeciesType const* r1, T2_ const& products, Real k)
{
    ReactionRule retval((ReactionRule::Reactants(r1)));
    retval.k() = k;
    std::for_each(boost::begin(products), boost::end(products),
            boost::bind(&ReactionRule::add_product, &retval, _1));
    return retval;
}

template<typename T2_>
inline ReactionRule new_reaction_rule(SpeciesType const* r1, SpeciesType const* r2, T2_ const& products, Real k)
{
    ReactionRule retval(ReactionRule::Reactants(r1, r2));
    retval.k() = k;
    std::for_each(boost::begin(products), boost::end(products),
            boost::bind(&ReactionRule::add_product, &retval, _1));
    return retval;
}

template<typename Tchar_, typename Ttraits_>
inline std::basic_ostream<Tchar_, Ttraits_>&
operator<<(std::basic_ostream<Tchar_, Ttraits_>& out, ReactionRule const& r)
{
    bool first;
    out << "ReactionRule(reactants={";
    first = true;
    BOOST_FOREACH (SpeciesType const* s, r.get_reactants())
    {
        if (!first)
        {
            out << ", ";
        }
        out << *s;
        first = false;
    }
    out << "}, products={";
    first = true;
    BOOST_FOREACH (SpeciesType const* s, r.get_products())
    {
        if (!first)
        {
            out << ", ";
        }
        out << *s;
        first = false;
    }
    out << "})";
    return out;
}

#endif /* REACTION_RULE_HPP */
