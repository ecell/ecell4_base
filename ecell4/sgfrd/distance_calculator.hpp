#ifndef ECELL4_SGFRD_DISTANCE_CALCULATOR
#define ECELL4_SGFRD_DISTANCE_CALCULATOR
#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/Circle.hpp>
#include <ecell4/core/Cone.hpp>
#include <boost/variant.hpp>
#include "distance.hpp"
#include "Shell.hpp"

namespace ecell4
{
namespace sgfrd
{

struct distance_calculator : public boost::static_visitor<Real>
{
    distance_calculator(const Real3& pos): pos_(pos){}

    template<typename shapeT, typename stridT>
    Real operator()(const Shell<shapeT, stridT>& sh)
    {
        return ecell4::sgfrd::distance<Real3, shapeT>()(pos_, sh.shape());
    }

  private:
    Real3 pos_;
};

template<typename T_polygon_traits, typename strID>
struct distance_calculator_on_surface : public boost::static_visitor<Real>
{
public:
    typedef T_polygon_traits polygon_traits;
    typedef strID structure_id_type;
    typedef ecell4::Polygon<polygon_traits> polygon_type;
    typedef typename polygon_type::face_id_type   face_id_type;
    typedef typename polygon_type::vertex_id_type vertex_id_type;

private:

    template<typename shapeT, typename stridT>
    struct is_2d_shell : boost::mpl::or_<
        boost::mpl::and_<
            boost::is_same<shapeT, ecell4::Circle>,
            boost::is_same<stridT, face_id_type>
            >,
        boost::mpl::and_<
            boost::is_same<shapeT, ecell4::ConicalSurface>,
            boost::is_same<stridT, vertex_id_type>
            >
        >
    {};

public:

    distance_calculator_on_surface(
            const std::pair<Real3, structure_id_type>& pos,
            const polygon_type& poly)
        : poly_(poly), pos_(pos)
    {}

    template<typename shapeT, typename stridT>
    typename boost::enable_if<is_2d_shell<shapeT, stridT>, Real>::type
    operator()(const Shell<shapeT, stridT>& sh) const
    {
        return poly_.distance(pos_, sh.get_surface_position()) - sh.size();
    }

private:

    polygon_type const& poly_;
    std::pair<Real3, structure_id_type> pos_;
};

} // sgfrd
} // ecell4
#endif// ECELL4_SGFRD_DISTANCE_CALCULATOR
