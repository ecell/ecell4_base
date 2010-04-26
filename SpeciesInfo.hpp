#ifndef SPECIES_INFO_HPP
#define SPECIES_INFO_HPP

#include <set>
#include <string>
#include <ostream>
#include "Defs.hpp"

template<typename Tid_, typename TD_, typename Tlen_, typename Tsurface_id_>
struct SpeciesInfo
{
    typedef Tid_ identifier_type;
    typedef TD_ D_type;
    typedef Tlen_ length_type;
    typedef Tsurface_id_ surface_id_type;

    identifier_type const& id() const
    {
        return id_;
    }

    length_type const& radius() const
    {
        return radius_;
    }

    length_type& radius()
    {
        return radius_;
    }

    surface_id_type const& surface_id() const
    {
        return surface_id_;
    }

    surface_id_type& surface_id()
    {
        return surface_id_;
    }

    D_type const& D() const
    {
        return diffusion_coef_;
    }

    D_type& D()
    {
        return diffusion_coef_;
    }

    bool operator==(SpeciesInfo const& rhs) const
    {
        return id_ == rhs.id() && diffusion_coef_ == rhs.D() &&
                radius_ == rhs.radius() && surface_id_ == rhs.surface_id();
    }

    bool operator!=(SpeciesInfo const& rhs) const
    {
        return !operator==(rhs);
    }

    SpeciesInfo(identifier_type const& id, D_type const& D = 0.,
                length_type const& r = 0., surface_id_type const& s = "") 
        : id_(id), diffusion_coef_(D), radius_(r), surface_id_(s) {}

    SpeciesInfo() {}

private:
    identifier_type id_;
    D_type diffusion_coef_;
    length_type radius_;
    surface_id_type surface_id_;
};

template<typename Tchar_, typename Ttraits_, typename Tid_, typename TD_, typename Tlen_, typename Tsurface_id_>
inline std::basic_ostream<Tchar_, Ttraits_>&
operator<<(std::basic_ostream<Tchar_, Ttraits_>& strm, const SpeciesInfo<Tid_, TD_, Tlen_, Tsurface_id_>& v)
{
    strm << "SpeciesInfo(id=" << v.id() << ", D=" << v.D() << ", radius=" << v.radius() << ", surface=" << v.surface_id() << ")";
    return strm;
}

#endif /* SPECIES_INFO_HPP */
