#ifndef SPECIES_HPP
#define SPECIES_HPP

#include <set>
#include <string>
#include <ostream>
#include "Defs.hpp"
#include "species_type_id.hpp"

template<typename T_>
struct species_info
{
    typedef T_ length_type;

    species_type const* type() const
    {
        return type_;
    }

    length_type const& radius() const
    {
        return radius_;
    }

    length_type& radius()
    {
        return radius_;
    }

    Real const& D() const
    {
        return diffusion_constant_;
    }

    Real& D()
    {
        return diffusion_constant_;
    }

    bool operator==(species const& rhs) const
    {
        return id_ == rhs.id() && name_ == rhs.name() &&
            diffusion_constant_ == rhs.D() && radius_ == rhs.radius();
    }

    bool operator!=(species const& rhs) const
    {
        return !operator==(rhs);
    }

    species_info(species_type const* type)
        : type_(type), diffusion_constant_(
            boost::lexical_cast<Real>((*type)["D"])),
          radius_(boost::lexical_cast<Real>((*type)["radius"])) {}

private:
    const species_type const* type_;
    double diffusion_constant_;
    length_type radius_;
};

template<typename Tchar_, typename Ttraits_ typename T_>
inline std::basic_ostream<Tchar_, Ttraits_>&
operator<<(std::basic_ostream<Tchar_, Ttraits_>& strm, const species<T_>& v)
{
    strm << "species_info(type=" << v.type() << ", D=" << v.D() << ", radius=" << v.radius() << ")";
    return strm;
}

#endif /* SPECIES_HPP */
