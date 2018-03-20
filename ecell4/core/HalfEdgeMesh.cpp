#include <ecell4/core/Polygon.hpp>

const Real Polygon<T_fid, T_vid, T_eid>::absolute_tolerance = 1e-12;
const Real Polygon<T_fid, T_vid, T_eid>::relative_tolerance = 1e-8;

void Polygon::assign(const std::vector<Triangle>& ts)
{
    const Real tol_abs2 = absolute_tolerance * absolute_tolerance;
    const Real tol_rel2 = relative_tolerance * relative_tolerance;

    vertices_.clear();
       faces_.clear();
       edges_.clear();
    this->total_area = 0.0;

    // prepair temporal data storage
    typedef std::pair<face_id_type, std::size_t>          fid_vidx_pair;
    typedef std::pair<Real3, std::vector<fid_vidx_pair> > tmp_vtx_type;
    typedef boost::flat_map<vertex_id_type, tmp_vtx_type> tmp_vertex_map;
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

        for(std::size_t i=0; i<3; ++i)
        {
            const Real3& v1 = triangle.vertices()[i];
            boost::optional<vertex_id_type> found_vtx = boost::none;

            // find near vertex
            for(typename tmp_vertex_map::const_iterator
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
                    vi->second.first = (v2 * vi->second.second.size() + v1) /
                                       (vi->second.second.size() + 1);
                    // assign face-id to the vertex
                    vi->second.second.push_back(std::make_pair(fid, i));
                    break;
                }
            }
            if(!found_vtx) // new vertices! add VertexID.
            {
                const vertex_id_type new_vid = vertices_.size();
                tmp_vtxs.push_back(std::make_pair(v1,
                        std::vector<fid_vidx_pair>(1, std::make_pair(fid, i))));
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
            this->edges_.update(eid, ed);

            fd.edges[i] = eid;
        }
        // set next of these 3 edges
        for(std::size_t i=0; i<3; ++i)
        {
            this->edges_.at(fd.edges[i]).next = fd.edges[i==2?0:i+1];
        }
        faces_.update(fid, fd);
    }

    // assign tmp_vtxs to this->vertices_.
    // by using tmp_vtxs, correct positions of the vertices.
    // here, face_datas are completed.
    for(typename tmp_vertex_map::const_iterator
            vi(tmp_vtxs.begin()), ve(tmp_vtxs.end()); vi != ve; ++vi)
    {
        const vertex_id_type vid = vi->first;
        const Real3          pos = vi->second.first;
        const std::vector<fid_vidx_pair>& face_pos = vi->second.second;

        vertex_data vd;
        vd.position = pos;
        this->vertices_.update(vid, vd);

        // * correct faces_ by using mean position
        // * set vertex.outgoing_edges
        for(typename std::vector<fid_vidx_pair>::const_iterator
                i(face_pos.begin()), e(face_pos.end()); i!=e; ++i)
        {
            const face_id_type fid = i->first;
            const std::size_t  idx = i->second;
            face_data& fd = this->faces_.at(fid);

            assert(vid == fd.vertices[idx]);

            vd.outgoing_edges.push_back(fd.edges[idx]);

            boost::array<Real3, 3> vs = this->faces_.at(fid).triangle.vertices();
            vs[idx] = pos; // update coordinate of Triangle
            this->faces_.at(fid).triangle = Triangle(vs);
        }
    }

    // set edge.length, edge.direction by using face.traingle
    for(typename face_container_type::const_iterator
            fi(this->faces_.begin()), fe(this->faces_.end()); fi != fe; ++fi)
    {
        const face_data& fd = fi->second;
        for(std::size_t i=0; i<3; ++i)
        {
            const edge_id_type eid = fd.edges[i];
            this->edges_.at(eid).length    = fd.triangle.length_of_edge_at(i);
            this->edges_.at(eid).direction = fd.triangle.edge_at(i);
        }
    }

    // search pairs of opposite edges & calculate edge.tilt.
    // here, edge_datas are completed.
    for(edge_id_type i=0; i<edges_.size(); ++i)
    {
        const vertex_id_type start  = target_of(i);
        const vertex_id_type target = target_of(next_of(next_of(i)));

        const std::vector<edge_id_type>& vd =
            this->vertices_[start].outgoing_edges;
        for(typename std::vector<edge_id_type>::const_iterator
                iter(vd.begin()), iend(vd.end()); iter != iend; ++iter)
        {
            const edge_id_type outgoing = *iter;
            if(target_of(outgoing) == target)
            {
                // found opposite edge! calculate tilt...
                this->edges_.at(i).opposite = outgoing;

                const face_id_type fid1 = face_of(i);
                const face_id_type fid2 = face_of(outgoing);
                const Real3 n1 = this->faces_.at(fid1).triangle.normal();
                const Real3 n2 = this->faces_.at(fid2).triangle.normal();
                const Real3 cr = cross_product(edges_.at(i).direction, n1);
                const Real  sg = dot_product(cr, n2);
                const Real ang = angle(n1, n2) * (sg > 0 ? 1 : -1);

                this->edges_.at(i       ).tilt = ang;
                this->edges_.at(outgoing).tilt = ang;

                break;
            }
        }
    }

    // set vertex_data.is_contiguous, angle by traversing edges.
    // here, vertex_data are completed.
    for(std::size_t i=0; i<vertices_.size(); ++i)
    {
        //TODO;
    }

    // make neighbor list for faces!
    // TODO
    return;
}




