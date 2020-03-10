#ifndef ECELL4_EGFRD_FILTERS_HPP
#define ECELL4_EGFRD_FILTERS_HPP

#include "geometry.hpp"

namespace ecell4
{
namespace egfrd
{

template<typename Toc_, typename Tfun_, typename Tsphere_>
class neighbor_filter
{
    using Iterator = typename Toc_::iterator;
    using Position = typename Toc_::position_type;

public:
    neighbor_filter(Tfun_& next, const Tsphere_& cmp)
        : next_(next), cmp_(cmp)
    {}

    inline void operator()(Iterator i, const Position& off) const
    {
        const auto dist(distance(shape(offset(i->second, off)), cmp_.position()));
        if (dist < cmp_.radius())
        {
            next_(i, dist);
        }
    }

private:
    Tfun_&   next_;
    Tsphere_ cmp_;
};

template<typename Toc_, typename Tfun_, typename Tsphere_>
inline void take_neighbor(Toc_& oc, Tfun_& fun, const Tsphere_& cmp)
{
    oc.each_neighbor(oc.index(cmp.position()),
                     neighbor_filter<Toc_, Tfun_, Tsphere_>(fun, cmp));
}

template<typename Toc_, typename Tfun_, typename Tsphere_>
inline void take_neighbor(Toc_ const& oc, Tfun_& fun, const Tsphere_& cmp)
{
    oc.each_neighbor(oc.index(cmp.position()),
                     neighbor_filter<Toc_ const, Tfun_, Tsphere_>(fun, cmp));
}

template<typename Toc_, typename Tfun_, typename Tsphere_>
inline void take_neighbor_cyclic(Toc_& oc, Tfun_& fun, const Tsphere_& cmp)
{
    oc.each_neighbor_cyclic(oc.index(cmp.position()),
            neighbor_filter<Toc_, Tfun_, Tsphere_>(fun, cmp));
}

template<typename Toc_, typename Tfun_, typename Tsphere_>
inline void take_neighbor_cyclic(Toc_ const& oc, Tfun_& fun, const Tsphere_& cmp)
{
    oc.each_neighbor_cyclic(oc.index(cmp.position()),
            neighbor_filter<Toc_ const, Tfun_, Tsphere_>(fun, cmp));
}

} // egfrd
} // ecell4
#endif /* ECELL4_EGFRD_FILTERS_HPP */
