#ifndef ECELL4_DISTANCE_MEASURER
#define ECELL4_DISTANCE_MEASURER
#include <boost/variant.hpp>
#include <boost/mpl/bool.hpp>
#include <ecell4/Circle.hpp>
#include <ecell4/Cone.hpp>

namespace ecell4
{
namespace sgfrd
{

template<typename T_polygon_traits>
class DistanceMeasurer : public boost::static_visitor<Real>
{
public:
    typedef T_polygon_traits traits_type;
    typedef Polygon<traits_type> polygon_type;
    typedef typename polygon_type::face_id_type   face_id_type;
    typedef typename polygon_type::edge_id_type   edge_id_type;
    typedef typename polygon_type::vertex_id_type vertex_id_type;

    template<typename T>
    struct is_1d_shell : boost::mpl::false_{};
    template<typename T>
    struct is_2d_shell : boost::mpl::false_{};
    template<typename T>
    struct is_3d_shell : boost::mpl::false_{};

public:

    DistanceMeasurer(const polygon_type& poly)
        : polygon_(poly)
    {}

    template<typename shellT, typename sidT>
    typename boost::enable_if<is_2d_shell<shellT>, Real>::type
    operator()(const std::pair<Real3, sidT>& pos, const shellT& sh) const
    {
        return polygon_.distance(pos, sh.get_surface_position()) - sh.size();
    }

    template<typename T_shell1, typename T_shell2>
    typename boost::enable_if<boost::mpl::and_<
        is_2d_shell<T_shell1>, is_2d_shell<T_shell2> >,
    Real>::type
    operator()(const shellT& sh1, const shellT& sh2) const
    {
        return polygon_.distance(sh1.get_surface_position(),
                                 sh2.get_surface_position()) -
               (sh1.size() + sh2.size());
    }

private:

    const polygon_type& polygon_;
};

template<typename T_pt>
template<typename shapeT, typename sidT>
struct DistanceMeasurer<T_pt>::is_2d_shell<Shell<shapeT, sidT> >
    : boost::mpl::or_<
        boost::mpl::and_<
            boost::is_same<shapeT, ecell4::Circle>,
            boost::is_same<sidT,   face_id_type> >,
        boost::mpl::and_<
            boost::is_same<shapeT, ecell4::ConicalSurface>,
            boost::is_same<sidT,   vertex_id_type> >
    >
{};

} // sgfrd
} // ecell4
#endif// ECELL4_DISTANCE_MEASURER
