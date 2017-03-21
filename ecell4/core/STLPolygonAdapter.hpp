#ifndef ECELL4_STL_POLYGON_ADAPTER
#define ECELL4_STL_POLYGON_ADAPTER
#include "Polygon.hpp"
#include "STLFileReader.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <set>

namespace ecell4
{

template<typename T_traits>
class STLPolygonAdapter
{
  public:
    typedef T_traits traits_type;
    typedef Polygon<traits_type> polygon_type;
    typedef typename polygon_type::index_type        index_type;
    typedef typename polygon_type::triangle_type     triangle_type;
    typedef typename polygon_type::face_id_type      face_id_type;
    typedef typename polygon_type::local_idx_type    local_idx_type;
    typedef typename polygon_type::vertex_id_type    vertex_id_type;
    typedef typename polygon_type::edge_id_type      edge_id_type;
    typedef typename polygon_type::vertex_descripter vertex_descripter;
    typedef typename polygon_type::edge_descripter   edge_descripter;
    typedef typename polygon_type::face_descripter   face_descripter;

  public:
    STLPolygonAdapter(const Real tol = 1e-7, const std::size_t max_face = 100)
        : tolerance(tol), max_faces_sharing_vertex(max_face)
    {}
    ~STLPolygonAdapter(){}

    boost::shared_ptr<polygon_type>
    make_polygon(const std::vector<StlTriangle>& triangles) const;

    void detect_edge_connections(  polygon_type& poly) const;
    void detect_vertex_connections(polygon_type& poly) const;

  private:

    const Real tolerance;
    const std::size_t max_faces_sharing_vertex;
};

template<typename T_traits>
boost::shared_ptr<typename STLPolygonAdapter<T_traits>::polygon_type>
STLPolygonAdapter<T_traits>::make_polygon(
        const std::vector<StlTriangle>& triangles) const
{
    boost::shared_ptr<polygon_type> polygon = boost::make_shared<polygon_type>();

    for(typename std::vector<StlTriangle>::const_iterator
            iter = triangles.begin(); iter != triangles.end(); ++iter)
        polygon->add_face(triangle_type(iter->vertices));

    this->detece_edge_connections(*polygon);
    this->detece_vertex_connections(*polygon);
    return polygon;
}

template<typename T_traits>
void STLPolygonAdapter<T_traits>::detect_edge_connections(polygon_type& poly) const
{
    std::set<edge_id_type> is_detected;
    for(std::size_t fidx = 0; fidx < poly.num_triangles(); ++fidx)
    {
        const face_id_type currentf(fidx);
        for(std::size_t eidx = 0; eidx < 3; ++eidx)
        {
            const edge_id_type current_edge =
                polygon_type::make_eid(currentf, local_idx_type(eidx));
            if(is_detected.count(current_edge) == 1) continue;

            const std::size_t start_v = eidx;
            const std::size_t end_v   = (eidx == 2) ? 0 : eidx + 1;

            const Real3 start_pos = poly.triangle_at(currentf).vertex_at(start_v);
            const Real3 end_pos = poly.triangle_at(currentf).vertex_at(end_v);

            bool found = false;
            for(std::size_t f = currentf+1; f < poly.num_triangles(); ++f)
            {
                for(std::size_t e = 0; e < 3; ++e)
                {
                    const Real start_pos_dist =
                        length(start_pos - poly.triangle_at(f).vertex_at(e));
                    if(start_pos_dist > tolerance)
                        continue;

                    const Real end_pos_dist =
                        length(end_pos - poly.triangle_at(f).vertex_at(
                                    e == 0 ? 2 : e-1));

                    if(end_pos_dist <= tolerance)
                    {
                        const edge_id_type detected = polygon_type::make_eid(
                            face_id_type(f), local_idx_type(e == 0 ? 2 : e-1));
                        poly.connect_edges(current_edge, detected);

                        is_detected.insert(detected);
                        found = true;
                        break;
                    }
                }
                if(found) break;
            }
            if(not found) throw std::logic_error("polygon is not closed");
        }
    }
    return;
}

template<typename T_traits>
void STLPolygonAdapter<T_traits>::detect_vertex_connections(polygon_type& poly) const
{
    std::set<vertex_id_type> is_detected;
    for(std::size_t fidx = 0; fidx < poly.num_faces(); ++fidx)
    {
        const face_id_type current_fid(fidx);
        for(std::size_t vidx = 0; vidx < 3; ++vidx)
        {
            const vertex_id_type current_vtx =
                polygon_type::make_vid(current_fid, local_idx_type(vidx));
            if(is_detected.count(current_vtx) == 1)
                continue;

            const Real3 vtx_pos = poly.triangle_at(fidx).vertex_at(vidx);

            std::vector<vertex_id_type> vertices_list;
            vertices_list.push_back(current_vtx);

            edge_id_type lookup =
                polygon_type::make_eid(current_fid, local_idx_type(vidx));

            std::size_t i=0;
            for(; i<max_faces_sharing_vertex; ++i)
            {
                const face_id_type   f(polygon_type::get_face_id(lookup));
                const local_idx_type e(polygon_type::get_local_index(lookup));
                const edge_id_type eid(
                    polygon_type::make_eid(f, local_idx_type(e==0?2:e-1)));
                lookup = poly.connecting_edge(eid); // update lookup face

                if(polygon_type::get_face_id(lookup) == current_fid)
                    break;

                const vertex_id_type detected =
                    polygon_type::make_vid(polygon_type::get_face_id(lookup),
                                           polygon_type::get_local_index(lookup));

                is_detected.insert(detected);
                vertices_list.push_back(detected);
            }
            if(i == max_faces_sharing_vertex)
                throw std::runtime_error("too many faces share a certain vertex");

            poly.connect_vertices(vertices_list);
        }
    }
    return ;
}

} // ecell4
#endif /* ECELL4_STL_POLYGON_ADAPTER */
