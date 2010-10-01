#ifndef REACTION_RECORD_HPP
#define REACTION_RECORD_HPP

#include <vector>
#include "twofold_container.hpp"
#include "utils/memberwise_compare.hpp"

template<typename Tpid_, typename Trid_>
class ReactionRecord
{
public:
    typedef Tpid_ particle_id_type;
    typedef Trid_ reaction_rule_id_type;
    typedef std::vector<particle_id_type> products_type;
    typedef twofold_container<particle_id_type> reactants_type;

public:
    ReactionRecord()
        : reaction_rule_id_(), reactants_(), products_() {}

    template<typename Tset>
    ReactionRecord(reaction_rule_id_type const& rid,
                   Tset const& products, 
                   particle_id_type const& p1)
        : reaction_rule_id_(rid), reactants_(p1),
          products_(boost::begin(products), boost::end(products)) {}

    template<typename Tset>
    ReactionRecord(reaction_rule_id_type const& rid,
                   Tset const& products,
                   particle_id_type const& p1, particle_id_type const& p2)
        : reaction_rule_id_(rid), reactants_(p1, p2),
          products_(boost::begin(products), boost::end(products)) {}

    // HEADS UP: move constructor!
    ReactionRecord(ReactionRecord const& that)
    {
        swap(const_cast<ReactionRecord&>(that));
    }

    reaction_rule_id_type const& reaction_rule_id() const
    {
        return reaction_rule_id_;
    }

    reactants_type const& reactants() const
    {
        return reactants_;
    }

    products_type const& products() const
    {
        return products_;
    }

    operator bool() const
    {
        return reactants_.size() != 0;
    }

    bool operator==(ReactionRecord const& rhs) const
    {
        return reaction_rule_id_ == rhs.reaction_rule_id() &&
               memberwise_compare(reactants_, rhs.reactants_) == 0 &&
               memberwise_compare(products_, rhs.products_) == 0;
    }

    bool operator!=(ReactionRecord const& rhs) const
    {
        return !operator==(rhs);
    }

    void swap(ReactionRecord& that)
    {
        std::swap(reaction_rule_id_, that.reaction_rule_id_);
        reactants_.swap(that.reactants_);
        products_.swap(that.products_);
    }

protected:
    reaction_rule_id_type reaction_rule_id_;
    reactants_type reactants_;
    products_type products_;
};

template<typename Tchar, typename Ttraits, typename Tpid, typename Trid>
inline std::basic_ostream<Tchar, Ttraits>&
operator<<(std::basic_ostream<Tchar, Ttraits>& out,
           ReactionRecord<Tpid, Trid> const& r)
{
    bool first;
    out << "ReactionRecord(reaction_rule_id=" << r.reaction_rule_id() << ", ";
    out << "reactants={";
    typedef typename ReactionRecord<Tpid, Trid>::reactants_type reactants_type;
    typedef typename ReactionRecord<Tpid, Trid>::products_type products_type;
    reactants_type const& reactants(r.reactants());
    for (typename boost::range_const_iterator<reactants_type>::type
            i(boost::begin(reactants)), e(boost::end(reactants));
         i != e; ++i)
    {
        if (!first)
        {
            out << ", ";
        }
        out << *i;
        first = false;
    }
    out << "}, products={";
    first = true;
    products_type const& products(r.products());
    for (typename boost::range_const_iterator<products_type>::type
                i(boost::begin(products)), e(boost::end(products));
         i != e; ++i)
    {
        if (!first)
        {
            out << ", ";
        }
        out << *i;
        first = false;
    }
    out << "})";
    return out;
}

#endif /* REACTION_RECORD_HPP */
