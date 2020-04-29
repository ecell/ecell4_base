#ifndef ECELL4_SGFRD_SHELL_VISITORS
#define ECELL4_SGFRD_SHELL_VISITORS
#include "Shell.hpp"
#include "SGFRDEvent.hpp"
#include <type_traits>

namespace ecell4
{
namespace sgfrd
{

struct minimal_eval_or : std::true_type
{
    static bool is_resolved(const bool v) {return v;}
};

struct minimal_eval_and : std::false_type
{
    static bool is_resolved(const bool v) {return !v;}
};

struct domain_id_setter : boost::static_visitor<void>
{
    DomainID did_;
    domain_id_setter(DomainID did) : did_(did){}

    template<typename shapeT, typename stridT>
    void operator()(Shell<shapeT, stridT>& shell) const
    {
        shell.domain_id() = did_;
    }
};

struct domain_id_getter : boost::static_visitor<DomainID>
{
    template<typename shapeT, typename stridT>
    DomainID operator()(const Shell<shapeT, stridT>& shell) const
    {
        return shell.domain_id();
    }
};

struct shell_size_getter : boost::static_visitor<Real>
{
    template<typename shapeT, typename stridT>
    Real operator()(const Shell<shapeT, stridT>& shell) const
    {
        return shell.size();
    }
};

struct shell_position_getter : boost::static_visitor<Real3>
{
    template<typename shapeT, typename stridT>
    Real3 operator()(const Shell<shapeT, stridT>& shell) const
    {
        return shell.position();
    }
};

struct inside_checker : boost::static_visitor<bool>
{
    typedef minimal_eval_or eval_manner;
    typedef ecell4::Polygon polygon_type;
    typedef polygon_type::FaceID FaceID;


    inside_checker(Real3 pos, Real rad, FaceID f, polygon_type const& p)
        : radius(rad), position(pos), fid(f), poly(p)
    {}

    //XXX: dispatch using shapeT::dimension to use 3D shells
    template<typename shapeT, typename stridT>
    bool operator()(const Shell<shapeT, stridT>& shell) const
    {
        return ecell4::polygon::distance(poly,
                std::make_pair(shell.position(), shell.structure_id()),
                std::make_pair(position, fid)) < shell.size() - radius;
    }

  private:

    Real   radius;
    Real3  position;
    FaceID fid;
    polygon_type const& poly;
};


} // sgfrd
} // ecell4
#endif// ECELL4_SGFRD_SHELL_VISITORS
