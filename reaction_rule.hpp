#ifndef REACTION_RULE_HPP
#define REACTION_RULE_HPP

#include <set>
#include <ostream>
#include <iostream>

#include <boost/bind.hpp>
#include <boost/range/size.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_categories.hpp>

#include "Defs.hpp"
#include "utils.hpp"
#include "species_type.hpp"

class reaction_rule
{
private:
    typedef std::set<species_type*> species_type_set;

public:
    class reactants
    {
    private:
        typedef boost::array<species_type*, 2> containing_type;
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

            species_type*& dereference() const
            {
                return cntnr_[idx_];
            }

        public:
            iterator(reactants& cntnr, size_type idx)
                : cntnr_(cntnr), idx_(idx) {}

            iterator(const_iterator const&);

        private:
            reactants& cntnr_;
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

            species_type* const& dereference() const
            {
                return cntnr_[idx_];
            }

        public:
            const_iterator(reactants const& cntnr, size_type idx)
                : cntnr_(cntnr), idx_(idx) {}

            const_iterator(iterator const& that)
                : cntnr_(that.cntnr_), idx_(that.idx_) {}

        private:
            reactants const& cntnr_;
            size_type idx_;
        };

    public:
        reactants()
        {
            items_[0] = 0;
            items_[1] = 0;
        }

        reactants(species_type* one)
        {
            BOOST_ASSERT(one);
            items_[0] = one;
            items_[1] = 0;
        }

        reactants(species_type* one, species_type* two)
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

        species_type*& operator[](std::size_t idx)
        {
            return items_[idx];
        }

        species_type* const& operator[](std::size_t idx) const
        {
            return items_[idx];
        }

    protected:
        containing_type items_;
    };

    typedef species_type_set::const_iterator species_type_iterator;
    typedef boost::iterator_range<species_type_iterator> species_type_range;

public:
    reactants const& get_reactants() const
    {
        return reactants_;
    }

    void add_product(species_type* s)
    {
        if (!products_.insert(s).second)
            throw already_exists(boost::lexical_cast<std::string>(*s));
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

    explicit reaction_rule(class reactants const& _reactants)
        : reactants_(_reactants) {}

private:
    reactants reactants_;
    species_type_set products_;
    Real k_;
};

inline reaction_rule::reactants::iterator::iterator(const_iterator const& that)
    : cntnr_(const_cast<reactants&>(that.cntnr_)), idx_(that.idx_)
{
}

inline bool operator<(reaction_rule::reactants const& lhs, reaction_rule::reactants const& rhs)
{
    return memberwise_compare(lhs, rhs) < 0;
}

inline bool operator<(reaction_rule const& lhs, reaction_rule const& rhs)
{
    if (lhs.k() < rhs.k())
        return true;

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

template<typename T2_>
inline reaction_rule new_reaction_rule(species_type* r1, T2_ const& products, Real k)
{
    reaction_rule retval((reaction_rule::reactants(r1)));
    retval.k() = k;
    std::for_each(boost::begin(products), boost::end(products),
            boost::bind(&reaction_rule::add_product, &retval, _1));
    return retval;
}

template<typename T2_>
inline reaction_rule new_reaction_rule(species_type* r1, species_type* r2, T2_ const& products, Real k)
{
    reaction_rule retval(reaction_rule::reactants(r1, r2));
    retval.k() = k;
    std::for_each(boost::begin(products), boost::end(products),
            boost::bind(&reaction_rule::add_product, &retval, _1));
    return retval;
}

template<typename Tchar_, typename Ttraits_>
inline std::basic_ostream<Tchar_, Ttraits_>&
operator<<(std::basic_ostream<Tchar_, Ttraits_>& out, reaction_rule const& r)
{
    bool first;
    out << "reaction_rule(reactants={";
    first = true;
    BOOST_FOREACH (species_type* s, r.get_reactants())
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
    BOOST_FOREACH (species_type* s, r.get_products())
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
