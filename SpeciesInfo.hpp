#ifndef SPECIES_INFO_HPP
#define SPECIES_INFO_HPP

#include <set>
#include <string>
#include <ostream>
#include "Defs.hpp"

template<typename Tid_, typename TD_, typename Tlen_>
struct SpeciesInfo
{
    typedef Tid_ identifier_type;
    typedef TD_ D_type;
    typedef Tlen_ length_type;
    typedef std::string surface_type;

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

    surface_type const& surface() const
    {
        return surface_;
    }

    surface_type& surface()
    {
        return surface_;
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
                radius_ == rhs.radius() && surface_ == rhs.surface();
    }

    bool operator!=(SpeciesInfo const& rhs) const
    {
        return !operator==(rhs);
    }

    SpeciesInfo(identifier_type const& id, D_type const& D = 0.,
                length_type const& r = 0., surface_type const& s = "") 
        : id_(id), diffusion_coef_(D), radius_(r), surface_(s) {}

    SpeciesInfo() {}

private:
    identifier_type id_;
    D_type diffusion_coef_;
    length_type radius_;
    surface_type surface_;
};

template<typename Tchar_, typename Ttraits_, typename Tid_, typename TD_, typename Tlen_>
inline std::basic_ostream<Tchar_, Ttraits_>&
operator<<(std::basic_ostream<Tchar_, Ttraits_>& strm, const SpeciesInfo<Tid_, TD_, Tlen_>& v)
{
    strm << "SpeciesInfo(id=" << v.id() << ", SpeciesTypeID=" << v.type_id() <<  ", D=" << v.D() << ", radius=" << v.radius() << ", surface=" << v.surface() << ")";
    return strm;
}

#endif /* SPECIES_INFO_HPP */
