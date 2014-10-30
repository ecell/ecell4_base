#ifndef SPECIES_INFO_HPP
#define SPECIES_INFO_HPP

#include <set>
#include <string>
#include <ostream>
#include "Defs.hpp"

template<typename Tid_, typename TD_, typename Tlen_, typename Tstructure_id_>
struct SpeciesInfo
{
    typedef Tid_ identifier_type;
    typedef TD_ D_type;
    typedef TD_ v_type;
    typedef Tlen_ length_type;
    typedef Tstructure_id_ structure_id_type;

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

    structure_id_type const& structure_id() const
    {
        return structure_id_;
    }

    structure_id_type& structure_id()
    {
        return structure_id_;
    }

    D_type const& D() const
    {
        return diffusion_coef_;
    }

    D_type& D()
    {
        return diffusion_coef_;
    }
    
    v_type const& v() const
    {
        return drift_velocity_;
    }

    v_type& v()
    {
        return drift_velocity_;
    }

    bool operator==(SpeciesInfo const& rhs) const
    {
        return id_ == rhs.id() && diffusion_coef_ == rhs.D() && drift_velocity_ == rhs.v() &&
                radius_ == rhs.radius() && structure_id_ == rhs.structure_id();
    }

    bool operator!=(SpeciesInfo const& rhs) const
    {
        return !operator==(rhs);
    }

    SpeciesInfo() {}

    SpeciesInfo(identifier_type const& id, D_type const& D = 0., 
                length_type const& r = 0., structure_id_type const& s = "", v_type const& v = 0.) 
        : id_(id), diffusion_coef_(D), drift_velocity_(v), radius_(r), structure_id_(s) {}

private:

    identifier_type id_;
    D_type diffusion_coef_;
    v_type drift_velocity_;
    length_type radius_;
    structure_id_type structure_id_;
};

template<typename Tchar_, typename Ttraits_, typename Tid_, typename TD_, typename Tlen_, typename Tstructure_id_>
inline std::basic_ostream<Tchar_, Ttraits_>&
operator<<(std::basic_ostream<Tchar_, Ttraits_>& strm, const SpeciesInfo<Tid_, TD_, Tlen_, Tstructure_id_>& s)
{
    strm << "SpeciesInfo(id=" << s.id() << ", D=" << s.D() << ", v=" << s.v() << ", radius=" << s.radius() << ", surface=" << s.structure_id() << ")";
    return strm;
}

#endif /* SPECIES_INFO_HPP */
