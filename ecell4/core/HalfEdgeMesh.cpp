#include <ecell4/core/HalfEdgeMesh.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/algorithm/cxx11/all_of.hpp>

namespace ecell4
{

const Real HalfEdgePolygon::absolute_tolerance = 1e-12;
const Real HalfEdgePolygon::relative_tolerance = 1e-8;

void HalfEdgePolygon::assign(const std::vector<Triangle>& ts)
{
    const Real tol_abs2 = absolute_tolerance * absolute_tolerance;
    const Real tol_rel2 = relative_tolerance * relative_tolerance;

    vertices_.clear();
       faces_.clear();
       edges_.clear();
    this->total_area_ = 0.0;

    // prepair temporal data storage
    typedef std::pair<face_id_type, std::size_t>          fid_vidx_pair;
    typedef std::pair<Real3, std::vector<fid_vidx_pair> > tmp_vtx_type;
    typedef boost::container::flat_map<vertex_id_type, tmp_vtx_type> tmp_vertex_map;
    tmp_vertex_map tmp_vtxs;

    // first, generate (FaceIDs for all triangles) and (EdgeIDs for all Edges).
    // and collect vertices that are at the same position.
    for(typename std::vector<Triangle>::const_iterator
            t_iter(ts.begin()), t_end(ts.end()); t_iter != t_end; ++t_iter)
    {
        const Triangle& triangle = *t_iter;
        this->total_area_ += triangle.area();

        const face_id_type fid = faces_.size();
        face_data fd;
        fd.triangle = triangle;

        for(std::size_t i=0; i<3; ++i)
        {
            const Real3& v1 = triangle.vertices()[i];
            boost::optional<vertex_id_type> found_vtx = boost::none;

            // find near vertex
            for(typename tmp_vertex_map::iterator
                    vi(tmp_vtxs.begin()), ve(tmp_vtxs.end()); vi != ve; ++vi)
            {
                const Real3&  v2 = vi->second.first;
                const Real dist2 =
                    length_sq(this->periodic_transpose(v1, v2) - v2);

                if(dist2 < tol_abs2 || dist2 < tol_rel2 * length_sq(v1))
                {
                    // vertex that locates near the vertex found
                    found_vtx = vi->first;

                    // calculating mean position...
                    vi->second.first = (v2 * vi->second.second.size() +
                                        this->apply_boundary(v1)) /
                                       (vi->second.second.size() + 1);
                    // assign face-id to the vertex
                    vi->second.second.push_back(std::make_pair(fid, i));
                    break;
                }
            }
            if(!found_vtx) // new vertices! add VertexID.
            {
                const vertex_id_type new_vid = tmp_vtxs.size();
                tmp_vtxs[new_vid] = std::make_pair(v1,
                        std::vector<fid_vidx_pair>(1, std::make_pair(fid, i)));
                found_vtx = new_vid;
            }
            fd.vertices[i] = *found_vtx; // store vertex id to face data
        }

        // make 3 edges around the face
        for(std::size_t i=0; i<3; ++i)
        {
            // in this point, edge length and direction are not fixed (because
            // vertex positions are corrected after all the faces are assigned).
            const edge_id_type eid = edges_.size();
            edge_data ed;
            ed.face   = fid;
            ed.target = fd.vertices[i==2?0:i+1];
            this->edges_.push_back(ed);

            fd.edges[i] = eid;
        }
        // set `next` of these 3 edges
        for(std::size_t i=0; i<3; ++i)
        {
            this->edges_.at(fd.edges[i]).next = fd.edges[i==2?0:i+1];
        }
        faces_.push_back(fd);
    }

    // * assign tmp_vtxs to this->vertices_
    // * set outgoing_edges
    for(typename tmp_vertex_map::const_iterator
            vi(tmp_vtxs.begin()), ve(tmp_vtxs.end()); vi != ve; ++vi)
    {
        const vertex_id_type vid = vi->first;
        const Real3          pos = vi->second.first;
        const std::vector<fid_vidx_pair>& face_pos = vi->second.second;

        vertex_data vd;
        vd.position = pos;

        // * set vertex.outgoing_edges
        for(typename std::vector<fid_vidx_pair>::const_iterator
                i(face_pos.begin()), e(face_pos.end()); i!=e; ++i)
        {
            const face_id_type fid = i->first;
            const std::size_t  idx = i->second;
            face_data& fd = this->faces_.at(fid);

            assert(vid == fd.vertices[idx]);
            vd.outgoing_edges.push_back(fd.edges[idx]);
        }
        this->vertices_.push_back(vd);
    }

    // * refine vertex positions
    for(typename face_container_type::iterator
            fi(this->faces_.begin()), fe(this->faces_.end()); fi != fe; ++fi)
    {
        face_data& fd = *fi;

        boost::array<Real3, 3> vs = fd.triangle.vertices();
        vs[0] = this->periodic_transpose(
                this->vertices_.at(fd.vertices[0]).position, vs[0]);
        vs[1] = this->periodic_transpose(
                this->vertices_.at(fd.vertices[1]).position, vs[1]);
        vs[2] = this->periodic_transpose(
                this->vertices_.at(fd.vertices[2]).position, vs[2]);
        fd.triangle = Triangle(vs);
    }

    // set edge.length, edge.direction by using face.traingle
    for(typename face_container_type::const_iterator
            fi(this->faces_.begin()), fe(this->faces_.end()); fi != fe; ++fi)
    {
        const face_data& fd = *fi;
        for(std::size_t i=0; i<3; ++i)
        {
            const edge_id_type eid = fd.edges[i];
            this->edges_.at(eid).length    = fd.triangle.length_of_edge_at(i);
            this->edges_.at(eid).direction = fd.triangle.edge_at(i);
        }
    }

    // search pairs of opposite edges & calculate edge.tilt.
    for(edge_id_type i=0; i<edges_.size(); ++i)
    {
        const vertex_id_type start  = this->target_of(i);
        const vertex_id_type target = this->target_of(
                this->next_of(this->next_of(i)));

        bool opposite_found = false;
        const std::vector<edge_id_type>& vd =
            this->vertices_.at(start).outgoing_edges;

        for(typename std::vector<edge_id_type>::const_iterator
                iter(vd.begin()), iend(vd.end()); iter != iend; ++iter)
        {
            const edge_id_type outgoing = *iter;
            if(this->target_of(outgoing) == target)
            {
                // found opposite edge! calculate tilt...
                this->edges_.at(i).opposite_edge = outgoing;

                const face_id_type fid1 = face_of(i);
                const face_id_type fid2 = face_of(outgoing);
                const Real3 n1 = this->faces_.at(fid1).triangle.normal();
                const Real3 n2 = this->faces_.at(fid2).triangle.normal();
                const Real3 cr = cross_product(edges_.at(i).direction, n1);
                const Real  sg = dot_product(cr, n2);
                const Real ang = angle(n1, n2) * (sg > 0 ? 1 : -1);

                this->edges_.at(i       ).tilt = ang;
                this->edges_.at(outgoing).tilt = ang;

                opposite_found = true;
                break;
            }
        }
        assert(opposite_found);
    }

    // set vertex_data.angle by traversing edges.
    for(std::size_t i=0; i<vertices_.size(); ++i)
    {
        const vertex_data& vtx = this->vertices_[i];
        std::vector<bool> is_found(vtx.outgoing_edges.size(), false);

        Real total_angle = 0.0;
        const edge_id_type start = vtx.outgoing_edges.front();
        edge_id_type current = start;
        do
        {
            for(std::size_t idx=0; idx<vtx.outgoing_edges.size(); ++idx)
            {
                if(vtx.outgoing_edges[idx] == current)
                {
                    is_found[idx] = true;
                    break;
                }
            }

            const face_data& f = this->faces_.at(this->face_of(current));
            bool angle_found = false;
            for(std::size_t idx=0; idx<3; ++idx)
            {
                if(f.edges[idx] == current)
                {
                    angle_found = true;
                    total_angle += f.triangle.angle_at(idx);
                    break;
                }
            }
            assert(angle_found);

            current = this->opposite_of(this->next_of(this->next_of(current)));
        }
        while(current != start);

        this->vertices_[i].apex_angle = total_angle;
        assert(boost::algorithm::all_of(
                    is_found.begin(), is_found.end(), boost::lambda::_1));
    }

    // make neighbor list for faces!
    // TODO
    return;
}


} // ecell4
