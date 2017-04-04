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
    typedef typename polygon_type::vertex_id_type    vertex_id_type;
    typedef typename polygon_type::edge_id_type      edge_id_type;
    typedef typename polygon_type::local_index_type  local_index_type;
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

    void detect_edge_connections(polygon_type& poly,
            const std::vector<face_id_type>& fid) const;

    void detect_vertex_connections(polygon_type& poly,
            const std::vector<face_id_type>& fid) const;

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

    std::vector<face_id_type> fids;
    fids.reserve(triangles.size());
    for(typename std::vector<StlTriangle>::const_iterator
            iter = triangles.begin(); iter != triangles.end(); ++iter)
    {
        fids.push_back(polygon->add_face(triangle_type(iter->vertices)));
    }

    this->detect_edge_connections(*polygon, fids);
    this->detect_vertex_connections(*polygon, fids);
    return polygon;
}

template<typename T_traits>
void STLPolygonAdapter<T_traits>::detect_edge_connections(polygon_type& poly,
        const std::vector<face_id_type>& fids) const
{
    std::set<local_index_type> is_detected;
    for(std::size_t fidx = 0; fidx < poly.num_triangles(); ++fidx)
    {
        const face_id_type currentf(fids.at(fidx));
        for(std::size_t eidx = 0; eidx < 3; ++eidx)
        {
            const local_index_type current_edge = std::make_pair(currentf, eidx);
            if(is_detected.count(current_edge) == 1) continue;

            const std::size_t start_v = eidx;
            const std::size_t end_v   = (eidx == 2) ? 0 : eidx + 1;

            const Real3 start_pos = poly.triangle_at(currentf).vertex_at(start_v);
            const Real3 end_pos = poly.triangle_at(currentf).vertex_at(end_v);

            bool found = false;
            for(std::size_t f = fidx+1; f < poly.num_triangles(); ++f)
            {
                for(std::size_t e = 0; e < 3; ++e)
                {
                    const Real start_pos_dist =
                        length(start_pos - poly.triangle_at(fids.at(f)).vertex_at(e));
                    if(start_pos_dist > tolerance)
                        continue;

                    const Real end_pos_dist =
                        length(end_pos - poly.triangle_at(fids.at(f)).vertex_at(
                                    e == 0 ? 2 : e-1));

                    if(end_pos_dist <= tolerance)
                    {
                        const local_index_type detected = std::make_pair(
                            fids.at(f), (e == 0 ? 2 : e-1));
                        poly.connect_edges(current_edge, detected);

                        is_detected.insert(detected);
                        found = true;
                        break;
                    }
                }
                if(found) break;
            }
            if(!found) throw std::logic_error("polygon is not closed");
        }
    }
    return;
}

template<typename T_traits>
void STLPolygonAdapter<T_traits>::detect_vertex_connections(polygon_type& poly,
        const std::vector<face_id_type>& fids) const
{
    std::set<local_index_type> is_detected;
    for(std::size_t fidx = 0; fidx < poly.num_faces(); ++fidx)
    {
        const face_id_type current_fid(fids.at(fidx));
        for(std::size_t vidx = 0; vidx < 3; ++vidx)
        {
            const local_index_type current_vtx = std::make_pair(current_fid, vidx);
            if(is_detected.count(current_vtx) == 1)
                continue;

            const Real3 vtx_pos = poly.triangle_at(current_fid).vertex_at(vidx);

            std::vector<local_index_type> vertices_list;
            vertices_list.push_back(current_vtx);

            local_index_type lookup = std::make_pair(current_fid, vidx);

            std::size_t i=0;
            for(; i<max_faces_sharing_vertex; ++i)
            {
                const face_id_type f(lookup.first);
                const index_type   e(lookup.second == 0 ? 2 : lookup.second-1);
                const edge_id_type eid(poly.get_edge_id(std::make_pair(f, e)));

                const face_id_type next_f = poly.adjacent_faces(f)[e];
                const boost::array<face_id_type, 3> adj(poly.adjacent_faces(next_f));
                const index_type next_e = std::distance(adj.begin(),
                        std::find(adj.begin(), adj.end(), f));
                lookup = std::make_pair(next_f, next_e);

                if(lookup.first == current_fid) break;

                const local_index_type detected = lookup;

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
