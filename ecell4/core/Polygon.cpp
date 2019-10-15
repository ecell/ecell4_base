#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/exceptions.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/container/static_vector.hpp>
#include <boost/format.hpp>

namespace ecell4
{

const Real Polygon::absolute_tolerance = 1e-12;
const Real Polygon::relative_tolerance = 1e-8;

void Polygon::assign(const std::vector<Triangle>& ts)
{
    constexpr Real pi = boost::math::constants::pi<Real>();
    const Real tol_abs2 = absolute_tolerance * absolute_tolerance;
    const Real tol_rel2 = relative_tolerance * relative_tolerance;

    vertices_.clear();
       faces_.clear();
       edges_.clear();
    this->total_area_ = 0.0;

    // prepair temporal data storage
    typedef std::pair<FaceID, std::size_t>                     fid_vidx_pair;
    typedef std::pair<Real3, std::vector<fid_vidx_pair>>       tmp_vtx_type;
    typedef boost::container::flat_map<VertexID, tmp_vtx_type> tmp_vertex_map;
    tmp_vertex_map tmp_vtxs;

    // first, generate (FaceIDs for all triangles) and (EdgeIDs for all Edges).
    // and collect vertices that are at the same position.
    for(const Triangle& triangle : ts)
    {
        this->total_area_ += triangle.area();

        const FaceID fid = FaceID(faces_.size());
        face_data fd;
        fd.triangle = triangle;

        for(std::size_t i=0; i<3; ++i)
        {
            const Real3& v1 = triangle.vertices()[i];
            boost::optional<VertexID> found_vtx = boost::none;

            // find near vertex
            for(auto& vid_vtx : tmp_vtxs)
            {
                const Real3&  v2 = vid_vtx.second.first;
                const Real dist2 =
                    length_sq(this->periodic_transpose(v1, v2) - v2);

                if(dist2 < tol_abs2 || dist2 < tol_rel2 * length_sq(v1))
                {
                    // vertex that locates near the vertex found
                    found_vtx = vid_vtx.first;
                    auto& vtx = vid_vtx.second;

                    // calculating mean position...
                    vtx.first = (v2 * vtx.second.size() + this->apply_boundary(v1)) /
                                (vtx.second.size() + 1);
                    // assign face-id to the vertex
                    vtx.second.push_back(std::make_pair(fid, i));
                    break;
                }
            }
            if(!found_vtx) // new vertices! add VertexID.
            {
                const VertexID new_vid = VertexID(tmp_vtxs.size());
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
            const EdgeID eid = EdgeID(edges_.size());
            edge_data ed;
            ed.face   = fid;
            ed.target = fd.vertices[i==2?0:i+1];
            this->edges_.push_back(ed);

            fd.edges[i] = eid;
        }
        // set `next` of these 3 edges
        for(std::size_t i=0; i<3; ++i)
        {
            this->edge_at(fd.edges[i]).next = fd.edges[i==2?0:i+1];
        }
        faces_.push_back(fd);
    }

    // * assign tmp_vtxs to this->vertices_
    // * set outgoing_edges without order
    for(const auto& vid_vtx : tmp_vtxs)
    {
        const VertexID                         vid = vid_vtx.first;
        const Real3                            pos = vid_vtx.second.first;
        const std::vector<fid_vidx_pair>& face_pos = vid_vtx.second.second;

        vertex_data vd;
        vd.position = pos;

        // * set vertex.outgoing_edges, but not sorted.
        for(const auto& fid_vidx : face_pos)
        {
            const FaceID      fid = fid_vidx.first;
            const std::size_t idx = fid_vidx.second;
            face_data& fd = this->face_at(fid);

            assert(vid == fd.vertices[idx]);
            vd.outgoing_edges.push_back(std::make_pair(fd.edges[idx], 0.0));
        }
        this->vertices_.push_back(vd);
    }

    // * refine vertex positions
    for(face_data& fd : this->faces_)
    {
        boost::array<Real3, 3> vs = fd.triangle.vertices();
        vs[0] = this->periodic_transpose(
                this->vertex_at(fd.vertices[0]).position, vs[0]);
        vs[1] = this->periodic_transpose(
                this->vertex_at(fd.vertices[1]).position, vs[1]);
        vs[2] = this->periodic_transpose(
                this->vertex_at(fd.vertices[2]).position, vs[2]);
        fd.triangle = Triangle(vs);
    }

    // set edge.length, edge.direction by using face.traingle
    for(const face_data& fd : this->faces_)
    {
        for(std::size_t i=0; i<3; ++i)
        {
            const EdgeID eid = fd.edges[i];
            this->edge_at(eid).length    = fd.triangle.length_of_edge_at(i);
            this->edge_at(eid).direction = fd.triangle.edge_at(i);
        }
    }

    // search pairs of opposite edges & calculate edge.tilt.
    for(std::size_t i=0; i<edges_.size(); ++i)
    {
        const EdgeID   eid = EdgeID(i);
        const VertexID start  = this->target_of(eid);
        const VertexID target = this->target_of(
                this->next_of(this->next_of(eid)));

        bool opposite_found = false;
        for(const auto& eid_angle : this->vertex_at(start).outgoing_edges)
        {
            const EdgeID outgoing = eid_angle.first;
            if(this->target_of(outgoing) == target)
            {
                // found opposite edge! calculate tilt...
                this->edge_at(eid).opposite_edge = outgoing;

                const FaceID fid1 = face_of(eid);
                const FaceID fid2 = face_of(outgoing);
                const Real3 n1 = this->face_at(fid1).triangle.normal();
                const Real3 n2 = this->face_at(fid2).triangle.normal();
                const Real3 cr = cross_product(this->edge_at(eid).direction, n1);
                const Real  sg = dot_product(cr, n2);
                const Real ang = calc_angle(n1, n2) * (sg > 0 ? 1 : -1);

                this->edges_.at(i     ).tilt = ang;
                this->edge_at(outgoing).tilt = ang;

                opposite_found = true;
                break;
            }
        }
        if(!opposite_found)
        {
            throw ecell4::NotSupported("The given polygon is not closed.");
        }
    }

    // set vertex_data.angle by traversing edges.
    for(vertex_data& vtx : this->vertices_)
    {
        const std::size_t num_edges = vtx.outgoing_edges.size();
        std::vector<EdgeID> outgoing_edges_tmp(vtx.outgoing_edges.size());
        for(std::size_t idx=0; idx<vtx.outgoing_edges.size(); ++idx)
        {
            outgoing_edges_tmp[idx] = vtx.outgoing_edges[idx].first;
        }
        vtx.outgoing_edges.clear();

        Real total_angle = 0.0;
        const EdgeID start = outgoing_edges_tmp.front();
        EdgeID current = start;
        do
        {
            {
                const auto found = std::find(outgoing_edges_tmp.begin(),
                        outgoing_edges_tmp.end(), current);
                assert(found != outgoing_edges_tmp.end());
                outgoing_edges_tmp.erase(found);
            }
            const face_data& f = this->face_at(this->face_of(current));
            Real angle = std::numeric_limits<Real>::max();
            for(std::size_t idx=0; idx<3; ++idx)
            {
                if(f.edges[idx] == current)
                {
                    angle = f.triangle.angle_at(idx);
                    break;
                }
            }
            assert(angle != std::numeric_limits<Real>::max());

            total_angle += angle;
            vtx.outgoing_edges.push_back(std::make_pair(current, angle));
            current = this->opposite_of(this->next_of(this->next_of(current)));
        }
        while(current != start);

        vtx.apex_angle = total_angle;

        if(!outgoing_edges_tmp.empty())
        {
            std::cout << outgoing_edges_tmp.size() << std::endl;
            for(const auto& oet: outgoing_edges_tmp)
            {
                const auto fid = this->face_of(oet);
                std::cout << oet << ": on " << fid << ", {";
                const face_data& f = this->face_at(fid);
                for(std::size_t idx=0; idx<3; ++idx)
                {
                    std::cout << f.triangle.vertex_at(idx);
                    if(idx != 2) {std::cout << ", ";}
                }
                std::cout << "}, ";
                for(std::size_t idx=0; idx<3; ++idx)
                {
                    if(f.edges[idx] == oet)
                    {
                        std::cout << f.triangle.vertex_at(idx) << " -> ";
                        std::cout << f.triangle.vertex_at(idx==2?0:idx+1);
                        std::cout << std::endl;
                        break;
                    }
                }
            }
        }
        if(!outgoing_edges_tmp.empty())
        {
            throw std::runtime_error("Polygon::assign: internal error: "
                    "cannot traverse all the outgoing edges from a vertex");
        }
        if(vtx.outgoing_edges.size() != num_edges)
        {
            throw std::runtime_error("Polygon::assign: internal error: "
                    "inconsistent number of outgoing edges");
        }
    }

    // make neighbor list for faces!
    for(face_data& face : this->faces_)
    {
        for(std::size_t i=0; i<3; ++i)
        {
            const VertexID vid = face.vertices[i];
            const Real3 v_pos  = face.triangle.vertex_at(i);
            const Real3 normal = face.triangle.normal();

            // counter clock wise
            {
                const Real3 ref_edge = face.triangle.edge_at(i==0?2:i-1) /
                    (-length(face.triangle.edge_at(i==0?2:i-1)));

                face.neighbor_ccw[i].clear();
                const auto start_edge  = face.edges[i];
                EdgeID   current_edge  = opposite_of(next_of(next_of(start_edge)));
                Real     current_angle = 0.0;
                while(current_edge != start_edge)
                {
                    const FaceID fid  = face_of(current_edge);

                    const std::size_t vidx0 = this->face_at(fid).index_of(vid);
                    const std::size_t vidx1 = this->face_at(fid).index_of(
                                              target_of(current_edge));
                    const std::size_t vidx2 = this->face_at(fid).index_of(
                                              target_of(next_of(current_edge)));
                    assert(vidx0 < 3);
                    assert(vidx1 < 3);
                    assert(vidx2 < 3);

                    const Real next_angle = current_angle +
                        this->face_at(fid).triangle.angle_at(vidx0);

                    boost::array<Real3, 3> unfolded;
                    unfolded[vidx0] = v_pos;
                    unfolded[vidx1] = v_pos +
                        rotate(current_angle, normal, ref_edge) *
                        length_of(current_edge);
                    unfolded[vidx2] = v_pos +
                        rotate(next_angle,    normal, ref_edge) *
                        length_of(next_of(next_of(current_edge)));

                    face.neighbors.push_back(fid);
                    face.neighbor_ccw[i].emplace_back(fid, Triangle(unfolded));

                    current_angle = next_angle;
                    current_edge  = opposite_of(next_of(next_of(current_edge)));
                }
            }

            // clock wise
            {
                const Real3 ref_edge = face.triangle.edge_at(i) /
                    length(face.triangle.edge_at(i));

                face.neighbor_cw[i].clear();
                const auto start_edge  = face.edges[i];
                EdgeID   current_edge  = next_of(opposite_of(start_edge));
                Real     current_angle = 0.0;
                while(current_edge != start_edge)
                {
                    const FaceID fid  = face_of(current_edge);

                    const std::size_t vidx0 = this->face_at(fid).index_of(vid);
                    const std::size_t vidx1 = this->face_at(fid).index_of(
                                              target_of(current_edge));
                    const std::size_t vidx2 = this->face_at(fid).index_of(
                                              target_of(next_of(current_edge)));
                    assert(vidx0 < 3);
                    assert(vidx1 < 3);
                    assert(vidx2 < 3);

                    const Real next_angle = current_angle +
                        this->face_at(fid).triangle.angle_at(vidx0);

                    boost::array<Real3, 3> unfolded;
                    unfolded[vidx0] = v_pos;
                    unfolded[vidx1] = v_pos +
                        rotate(-next_angle, normal, ref_edge) *
                        length_of(current_edge);
                    unfolded[vidx2] = v_pos +
                        rotate(-current_angle, normal, ref_edge) *
                        length_of(next_of(next_of(current_edge)));

                    face.neighbors.push_back(fid);
                    face.neighbor_cw[i].emplace_back(fid, Triangle(unfolded));

                    current_angle = next_angle;
                    current_edge  = next_of(opposite_of(current_edge));
                }
            }
        }
        std::sort(face.neighbors.begin(), face.neighbors.end());
        const auto last = std::unique(face.neighbors.begin(), face.neighbors.end());
        face.neighbors.erase(last, face.neighbors.end());
    }
    return;
}

Real Polygon::distance_sq(const std::pair<Real3, FaceID>& pos1,
                          const std::pair<Real3, FaceID>& pos2) const
{
    typedef utils::pair_first_element_unary_predicator<FaceID, Triangle>
            face_finder_type;
    constexpr Real pi = boost::math::constants::pi<Real>();

    // if two particles are on the same face, return just a 3D distance.
    if(pos1.second == pos2.second)
    {
        return length_sq(pos2.first - pos1.first);
    }

    // If positions are on different faces, there can be several cases.
    // 1.)  ______
    //     /\    /\   | The positions are on the faces connected by a vertex.
    //    /  \  /p2\  | using p1-vtx-p2 angle and the low-of-cosines, the
    //   /____\/____\ | minimum distance on the surface can be calculated.
    //        /\    / |
    //       /p1\  /  |
    //      /____\/   |
    //
    // 2.)     ______
    //        /\ p2 / | The positions are on the faces that are not connected
    //       /  \  /  | by any vertex. There can be several options to unfold
    //      /____\/   | the polygon to make the particles on the same plane.
    //     /\    /    | In this case, finding the minimum path is too difficult
    //    /p1\  /     | to use in simulation, so just returns inf. In the SGFRD
    //   /____\/      | simulation, movement from p1 to p2 is inhibited.
    //
    // 3.)  ______
    //     /\    /\   | The positions are on the faces connected by a vertex
    //    /p2\  /  \  | and the apex angle exceeds 360 degree. There can be 2
    //   /____\/____\ | pathes, counter-clockwise and clockwise, and both angle
    //        /\    / | exceeds 180 degree. In this case, the minimum distance
    //       /p1\  /  | pathway goes across the vertex.
    //      /____\/   |
    //
    // 4.)  ......... ______
    //     /\connected\    /\   | The particles are on the faces connected by a
    //    /  \ <=====> \  /  \  | vertex and the apex angle exceeds 360 degree.
    //   /____\.........\/____\ | And the triangles overlaps when unfolded.
    //   \    /\        /\    / | In this case, we must use the minimum angle
    //    \  /p2\      /p1\  /  | because just unfolding triangle makes minimum
    //     \/____\    /____\/   | pathway shorter than the 3D distance.
    //
    // 5.)
    //     \`.p1           | If a polygon is severely deformed, the minimum
    //      \o`.           | angle pathway can protrude the face. To avoid this,
    //       ^  `.         | It is required to check the minimum angle pathway
    //       |\___`.vertex | goes across all the edges.
    //       |/   .'       |
    //       v  .'         |
    //      /o.'           |
    //     /.'p2           |
    //

    const Real3& p1 = pos1.first;
    const Real3& p2 = pos2.first;
    const FaceID f1 = pos1.second;
    const FaceID f2 = pos2.second;

    // for comparison
    const auto min_edge_length = std::min(edge_length_[0],
            std::min(edge_length_[1], edge_length_[2]));
    const auto rel_tol = relative_tolerance * min_edge_length;

    const face_data& face = face_at(f1);
    const Real3&   normal = face.triangle.normal();

    boost::container::static_vector<Real3, 3> connecting_vtxs;

    Real distance_sq = std::numeric_limits<Real>::infinity();
    for(std::size_t i=0; i<3; ++i)
    {
        // counter clockwise
        //
        //          ^ vertices[i]
        // edge[i] /|\
        //        / |~\
        //       /  o  \ next(next(egde[i]))
        //      v______>\
        //    next(edge[i])

        const VertexID    vid = face.vertices[i];
        const Real3&     vpos = position_at(vid);
        const Real3     vtop1 = this->periodic_transpose(p1, vpos) - vpos;
        const Real3     vtop2 = this->periodic_transpose(p2, vpos) - vpos;
        const Real3     p1tov = vtop1 * (-1.0);

        // check p1 or p2 are exactly on the vertex.
        // If they are on, it causes NaN because the length of v->p vector is 0.
        const Real vtop1_len = length(p1tov);
        const Real vtop2_len = length(vtop2);
        if(vtop1_len < relative_tolerance * min_edge_length)
        {
            // pos1 locates exactly on the vtx. distance from pos1 to pos2 is
            // equal to the distance from vtx to pos2.

            // if the face on which pos2 locates has the vertex, the squared
            // distance is just vtop2_len^2.

            const auto& vtxs = this->vertices_of(f2);
            if(std::find(vtxs.begin(), vtxs.end(), vid) != vtxs.end())
            {
                distance_sq = std::min(distance_sq, vtop2_len * vtop2_len);
            }
            continue;
        }
        if(vtop2_len < relative_tolerance * min_edge_length)
        {
            // pos2 locates exactly on the vtx. distance from pos1 to pos2 is
            // equal to the distance from pos1 to vtx.

            const auto& vtxs = this->vertices_of(f2);
            if(std::find(vtxs.begin(), vtxs.end(), vid) != vtxs.end())
            {
                distance_sq = std::min(distance_sq, vtop1_len * vtop1_len);
            }
            continue;
        }

        // calc the initial angle
        const auto init_angle = calc_angle(p1tov, face.triangle.edge_at(i==0?2:i-1));
        Real angle = init_angle;
        const Real apex_angle = apex_angle_at(vid);

        // ------------------------------------------------------------------
        // search `f2`
        bool connected = false;
        for(const auto& neighbor : face.neighbor_ccw[i])
        {
            assert(neighbor.first != f1);

            const auto& nface = face_at(neighbor.first);
            const auto  vidx  = nface.index_of(vid);
            if(neighbor.first == f2)
            {
                angle += calc_angle(vtop2, nface.triangle.edge_at(vidx));
                connected = true;
                break;
            }
            angle += nface.triangle.angle_at(vidx);
        }
        if(!connected)
        {
            continue;
        }
        connecting_vtxs.push_back(vpos);

        // ------------------------------------------------------------------
        // calculate the minimum angle

        const Real angle_ccw = angle;
        const Real angle_cw  = apex_angle - angle;
        assert(angle_cw >= 0.0);
        const Real min_angle = std::min(angle_ccw, angle_cw);

        // skip case 3 (theta > 180 deg).
        if(min_angle > pi) {continue;}

        // if the minimum angle < 180 degree, its case 1.
        // calculate distance using the low of cosine.

        // Before calculating the distance, check whether this is the case 5.
        // If a vertex locates inside of a triangle formed by a vertex, p1, and
        // p2, then it is case 5 and the pathway is not available.

        bool is_case5 = false;
        if(angle_ccw < angle_cw) // the min dist pathway is counter-clockwise
        {
            const auto sin_ccw = std::sin(angle_ccw);
            const auto cos_ccw = std::cos(angle_ccw);

            angle = init_angle;
            for(const auto& neighbor : face.neighbor_ccw[i])
            {
                if(neighbor.first == f2) {break;}

                const auto& nface = face_at(neighbor.first);
                const auto  vidx  = nface.index_of(vid);

                const auto sin_agl = std::sin(angle);
                const auto cos_agl = std::cos(angle);

                // sin of the opposite angle (angle_ccw - angle)
                const auto sin_opp = sin_ccw * cos_agl - cos_ccw * sin_agl;

                const auto threshold = vtop1_len * vtop2_len * sin_ccw /
                            (vtop1_len * sin_agl + vtop2_len * sin_opp);

                if(nface.triangle.length_of_edge_at(vidx) < threshold)
                {
                    is_case5 = true;
                    break;
                }
                angle += nface.triangle.angle_at(vidx);
            }
        }
        else // the minimum distance pathway is clockwise
        {
            const auto sin_cw = std::sin(angle_cw);
            const auto cos_cw = std::cos(angle_cw);

            angle = face.triangle.angle_at(i) - init_angle;
            for(const auto& neighbor : face.neighbor_cw[i])
            {
                if(neighbor.first == f2) {break;}

                const auto& nface = face_at(neighbor.first);
                const auto  vidx  = nface.index_of(vid);

                const auto sin_agl = std::sin(angle);
                const auto cos_agl = std::cos(angle);

                // sin of the opposite angle (angle_cw - angle)
                const auto sin_opp = sin_cw * cos_agl - cos_cw * sin_agl;

                const auto threshold = vtop1_len * vtop2_len * sin_cw /
                            (vtop1_len * sin_agl + vtop2_len * sin_opp);

                const auto edge_idx = (vidx==0) ? 2 : vidx-1;
                if(nface.triangle.length_of_edge_at(edge_idx) < threshold)
                {
                    is_case5 = true;
                    break;
                }
                angle += nface.triangle.angle_at(vidx);
            }
        }
        if(is_case5) // no available path.
        {
            continue;
        }

        Real rotation_angle = min_angle;
        if(angle_cw < angle_ccw)
        {
            rotation_angle *= -1; // rotate it clockwise
        }
        const auto vtop2_unf = rotate(rotation_angle, normal, vtop1) *
                               (vtop2_len / vtop1_len);

        distance_sq = std::min(distance_sq, length_sq(vtop1 - vtop2_unf));
    }
    if(distance_sq != std::numeric_limits<Real>::infinity())
    {
        // distance_sq is updated! It founds the minimum path (case 1).
        return distance_sq;
    }
    if(connecting_vtxs.empty())
    {
        // if `f1` and `f2` are not connected via any vertex, the distance
        // cannot be calculated (case 2).
        return std::numeric_limits<Real>::infinity();
    }

    // Here, the positions are connected via some vertices but the minimum
    // distance is non-trivial. It is case 3. return distance that passes
    // through the vertex that connects `f1` and `f2`.

    // here, distance_sq is still infinity.
    for(const Real3 vpos : connecting_vtxs)
    {
        const Real  l1   = length(this->periodic_transpose(p1, vpos) - vpos);
        const Real  l2   = length(this->periodic_transpose(p2, vpos) - vpos);
        distance_sq = std::min(distance_sq, (l1 + l2) * (l1 + l2));
    }
    return distance_sq;
}

Real Polygon::distance_sq(const std::pair<Real3, VertexID>& pos1,
                          const std::pair<Real3, FaceID>&   pos2) const
{
    for(const auto el : this->vertex_at(pos1.second).outgoing_edges)
    {
        const EdgeID eid = el.first;
        const FaceID fid = this->face_of(eid); // face that connects to the vtx
        if(fid == pos2.second)
        {
            // on the same face.
            return length_sq(
                this->periodic_transpose(pos1.first, pos2.first) - pos2.first);
        }
        //      ______ pos1
        //     ^     /
        //    / \   /
        //   /adj\ /
        //  /_____v

        const FaceID adj = this->face_of(this->opposite_of(this->next_of(eid)));
        if(adj == pos2.second)
        {
            // redirects to distance_sq({Real3, FaceID}, {Real3, FaceID})
            // considering pos1 as a position on a face
            return distance_sq(std::make_pair(pos1.first, fid), pos2);
        }
    }
    return std::numeric_limits<Real>::infinity();
}

Real Polygon::distance_sq(const std::pair<Real3, VertexID>& pos1,
                          const std::pair<Real3, VertexID>& pos2) const
{
    if(pos1.second == pos2.second)
    {
        return 0.0;
    }
    for(const auto el : this->vertex_at(pos1.second).outgoing_edges)
    {
        const EdgeID eid = el.first;
        if(this->target_of(eid) == pos2.second)
        {
            // directly connected.
            const Real l = length_of(eid);
            return l * l;
        }

        const EdgeID inbetween = this->opposite_of(this->next_of(eid));
        if(this->target_of(this->next_of(inbetween)) == pos2.second)
        {
            const FaceID fid = this->face_of(eid);
            const FaceID adj = this->face_of(inbetween);

            // redirects to distance_sq(face, face)
            return distance_sq(std::make_pair(pos1.first, fid),
                               std::make_pair(pos2.first, adj));
        }
    }
    return std::numeric_limits<Real>::infinity();
}

Real Polygon::distance(const std::pair<Real3, VertexID>& pos1,
                       const std::pair<Real3, VertexID>& pos2) const
{
    if(pos1.second == pos2.second)
    {
        return 0.0;
    }
    for(const auto el : this->vertex_at(pos1.second).outgoing_edges)
    {
        const EdgeID eid = el.first;
        if(this->target_of(eid) == pos2.second)
        {
            // directly connected.
            return length_of(eid);
        }

        const EdgeID inbetween = this->opposite_of(this->next_of(eid));
        if(this->target_of(this->next_of(inbetween)) == pos2.second)
        {
            const FaceID fid = this->face_of(eid);
            const FaceID adj = this->face_of(inbetween);

            // redirects to distance(face, face)
            return distance(std::make_pair(pos1.first, fid),
                            std::make_pair(pos2.first, adj));
        }
    }
    return std::numeric_limits<Real>::infinity();
}

// return direction from pos1 to pos2
Real3 Polygon::direction(const std::pair<Real3, FaceID>& pos1,
                         const std::pair<Real3, FaceID>& pos2) const
{
    constexpr Real pi = boost::math::constants::pi<Real>();
    typedef utils::pair_first_element_unary_predicator<FaceID, Triangle>
            face_finder_type;

    if(pos1.second == pos2.second)
    {
        return pos2.first - pos1.first;
    }

    const Real3& p1 = pos1.first;
    const Real3& p2 = pos2.first;
    const FaceID f1 = pos1.second;
    const FaceID f2 = pos2.second;

    // for comparison
    const auto min_edge_length = std::min(edge_length_[0],
            std::min(edge_length_[1], edge_length_[2]));
    const auto rel_tol = relative_tolerance * min_edge_length;

    const face_data& face = face_at(f1);
    const Real3&   normal = face.triangle.normal();

    boost::container::static_vector<Real3, 3> connecting_vtxs;

    Real distance_sq = std::numeric_limits<Real>::infinity();
    Real3 retval(0,0,0);
    for(std::size_t i=0; i<3; ++i)
    {
        // counter clockwise
        //
        //          ^ vertices[i]
        // edge[i] /|\
        //        / |~\
        //       /  o  \ next(next(egde[i]))
        //      v______>\
        //    next(edge[i])

        const VertexID    vid = face.vertices[i];
        const Real3&     vpos = position_at(vid);
        const Real3     vtop1 = this->periodic_transpose(p1, vpos) - vpos;
        const Real3     vtop2 = this->periodic_transpose(p2, vpos) - vpos;
        const Real3     p1tov = vtop1 * (-1.0);

        // check p1 or p2 are exactly on the vertex.
        // If they are on, it causes NaN because the length of v->p vector is 0.
        const Real vtop1_len = length(p1tov);
        const Real vtop2_len = length(vtop2);
        if(vtop1_len < relative_tolerance * min_edge_length)
        {
            // pos1 locates exactly on the vtx. distance from pos1 to pos2 is
            // equal to the distance from vtx to pos2.
            const auto& vtxs = this->vertices_of(f2);
            if(std::find(vtxs.begin(), vtxs.end(), vid) == vtxs.end())
            {
                continue;
            }
            return vtop2;
        }
        if(vtop2_len < relative_tolerance * min_edge_length)
        {
            // pos2 locates exactly on the vtx. distance from pos1 to pos2 is
            // equal to the distance from pos1 to vtx.
            const auto& vtxs = this->vertices_of(f2);
            if(std::find(vtxs.begin(), vtxs.end(), vid) == vtxs.end())
            {
                continue;
            }
            return vtop1;
        }

        // calc the initial angle
        const auto init_angle = calc_angle(p1tov, face.triangle.edge_at(i==0?2:i-1));
        Real angle = init_angle;
        const Real apex_angle = apex_angle_at(vid);

        // ------------------------------------------------------------------
        // search `f2`
        bool connected = false;
        for(const auto& neighbor : face.neighbor_ccw[i])
        {
            assert(neighbor.first != f1);

            const auto& nface = face_at(neighbor.first);
            const auto  vidx  = nface.index_of(vid);
            if(neighbor.first == f2)
            {
                angle += calc_angle(vtop2, nface.triangle.edge_at(vidx));
                connected = true;
                break;
            }
            angle += nface.triangle.angle_at(vidx);
        }
        if(!connected)
        {
            continue;
        }
        connecting_vtxs.push_back(vpos);

        // ------------------------------------------------------------------
        // calculate the minimum angle

        const Real angle_ccw = angle;
        const Real angle_cw  = apex_angle - angle;
        assert(angle_cw >= 0.0);
        const Real min_angle = std::min(angle_ccw, angle_cw);

        // skip case 3 (theta > 180 deg).
        if(min_angle > pi) {continue;}

        // if the minimum angle < 180 degree, its case 1.
        // calculate distance using the low of cosine.

        // Before calculating the distance, check whether this is the case 5.
        // If a vertex locates inside of a triangle formed by a vertex, p1, and
        // p2, then it is case 5 and the pathway is not available.

        bool is_case5 = false;
        if(angle_ccw < angle_cw) // the min dist pathway is counter-clockwise
        {
            const auto sin_ccw = std::sin(angle_ccw);
            const auto cos_ccw = std::cos(angle_ccw);

            angle = init_angle;
            for(const auto& neighbor : face.neighbor_ccw[i])
            {
                if(neighbor.first == f2) {break;}

                const auto& nface = face_at(neighbor.first);
                const auto  vidx  = nface.index_of(vid);

                const auto sin_agl = std::sin(angle);
                const auto cos_agl = std::cos(angle);

                // sin of the opposite angle (angle_ccw - angle)
                const auto sin_opp = sin_ccw * cos_agl - cos_ccw * sin_agl;

                const auto threshold = vtop1_len * vtop2_len * sin_ccw /
                            (vtop1_len * sin_agl + vtop2_len * sin_opp);

                if(nface.triangle.length_of_edge_at(vidx) < threshold)
                {
                    is_case5 = true;
                    break;
                }
                angle += nface.triangle.angle_at(vidx);
            }
        }
        else // the minimum distance pathway is clockwise
        {
            const auto sin_cw = std::sin(angle_cw);
            const auto cos_cw = std::cos(angle_cw);

            angle = face.triangle.angle_at(i) - init_angle;
            for(const auto& neighbor : face.neighbor_cw[i])
            {
                if(neighbor.first == f2) {break;}

                const auto& nface = face_at(neighbor.first);
                const auto  vidx  = nface.index_of(vid);

                const auto sin_agl = std::sin(angle);
                const auto cos_agl = std::cos(angle);

                // sin of the opposite angle (angle_cw - angle)
                const auto sin_opp = sin_cw * cos_agl - cos_cw * sin_agl;

                const auto threshold = vtop1_len * vtop2_len * sin_cw /
                            (vtop1_len * sin_agl + vtop2_len * sin_opp);

                const auto edge_idx = (vidx==0) ? 2 : vidx-1;
                if(nface.triangle.length_of_edge_at(edge_idx) < threshold)
                {
                    is_case5 = true;
                    break;
                }
                angle += nface.triangle.angle_at(vidx);
            }
        }
        if(is_case5) // no available path.
        {
            continue;
        }

        Real rotation_angle = min_angle;
        if(angle_cw < angle_ccw)
        {
            rotation_angle *= -1; // rotate it clockwise
        }
        const auto vtop2_unf = rotate(rotation_angle, normal, vtop1) *
                               (vtop2_len / vtop1_len);
        const auto dir       = vtop2_unf - vtop1; // from p1 to p2
        const auto dist      = length_sq(dir);

        if(dist < distance_sq)
        {
            distance_sq = dist;
            retval      = dir;
        }
    }
    if(distance_sq != std::numeric_limits<Real>::infinity())
    {
        // distance_sq is updated! It founds the minimum path (case 1).
        return retval;
    }
    // otherwise, there is no available direction between the positions.
    throw std::runtime_error((boost::format(
        "polygon::direction: couldn't find the min path between "
        "%1% on %2% <-> %3% on %4%") % pos1.first % pos1.second %
        pos2.first % pos2.second).str());
}

std::pair<Real3, Polygon::FaceID> Polygon::travel(
        const std::pair<Real3, FaceID>& pos, const Real3& disp) const
{
    const Real3& p = pos.first;
    const FaceID f = pos.second;
    const Real3 np = p + disp;

    const face_data& fd = this->face_at(f);
    const Triangle& tri = fd.triangle;
    const Barycentric b2(to_barycentric(np, tri));

    // if pos + disp is inside of the current face, just return the sum.
    if(::ecell4::is_inside(b2))
    {
        // to avoid numerical error that make the particle goes outside of the
        // face that the particle belongs, use `to_absolute`.
        return std::make_pair(to_absolute(b2, tri), f);
    }

    const Barycentric b1(to_barycentric(p, tri));
    const Barycentric db(b2 - b1);
    const std::pair<std::size_t, Real> cs   = first_cross_edge(b1, db);
    const std::pair<FaceID, Triangle>& next = fd.neighbor_cw[cs.first].front();

    // if the position is inside of the adjacent face, return the position
    // reconstructed from unfolded-Barycentric by using folded Triangle.
    const Barycentric unfolded_b(to_barycentric(np, next.second));
    if(::ecell4::is_inside(unfolded_b, 1e-8))
    {
        Real3 nxt = to_absolute(
                force_put_inside(unfolded_b), this->triangle_at(next.first));

        if(!this->is_inside_of_boundary(nxt))
        {
            const Real xlim    = this->edge_length_[0];
            const Real ylim    = this->edge_length_[1];
            const Real zlim    = this->edge_length_[2];
            const Real rel_tol = relative_tolerance;

            // allow numerical error in relative tolerance
            if(nxt[0] <  0.0 && std::abs(nxt[0]) < rel_tol * xlim) {nxt[0] = 0.0;}
            if(nxt[1] <  0.0 && std::abs(nxt[1]) < rel_tol * ylim) {nxt[1] = 0.0;}
            if(nxt[2] <  0.0 && std::abs(nxt[2]) < rel_tol * zlim) {nxt[2] = 0.0;}

            if(nxt[0] > xlim && nxt[0] - xlim < rel_tol * xlim) {nxt[0] = xlim;}
            if(nxt[1] > ylim && nxt[1] - ylim < rel_tol * ylim) {nxt[1] = ylim;}
            if(nxt[2] > zlim && nxt[2] - zlim < rel_tol * zlim) {nxt[2] = zlim;}

            if(!this->is_inside_of_boundary(nxt))
            {

                std::cerr << "travel: initial pos          = " << p   << " on " << f            << std::endl;
                std::cerr << "travel: initial face         = " << f << " -> " << this->triangle_at(f) << std::endl;
                std::cerr << "travel: initial disp         = " << disp                          << std::endl;
                std::cerr << "travel: next face (unfolded) = " << next.second                   << std::endl;
                std::cerr << "travel: next barycentric crd = " << unfolded_b                    << std::endl;
                std::cerr << "travel: next bary (inside)   = " << force_put_inside(unfolded_b)  << std::endl;
                std::cerr << "travel: next face            = " << this->triangle_at(next.first) << std::endl;
                std::cerr << "travel: next pos             = " << nxt << " on " << next.first   << std::endl;
                throw std::runtime_error("out of bound");
            }
        }
        return std::make_pair(nxt, next.first);
        // use folded (normal) Triangle, NOT next.second
    }

    // stride over not only the edge but adjacent face.
    // XXX to make it sure that `on_edge` should be on the edge under the PBC,
    //     to_absolute is used with the next triangle.
    const Barycentric on_edge_b(to_barycentric(p + disp * cs.second, next.second));
    const Real3 edge_over = direction_of(fd.edges[cs.first]);
    Real3 next_pos  = to_absolute(on_edge_b, this->triangle_at(next.first));

    if(!this->is_inside_of_boundary(next_pos))
    {
        const Real xlim    = this->edge_length_[0];
        const Real ylim    = this->edge_length_[1];
        const Real zlim    = this->edge_length_[2];
        const Real rel_tol = relative_tolerance;

        // allow numerical error in relative tolerance
        if(next_pos[0] <  0.0 && std::abs(next_pos[0]) < rel_tol * xlim) {next_pos[0] = 0.0;}
        if(next_pos[1] <  0.0 && std::abs(next_pos[1]) < rel_tol * ylim) {next_pos[1] = 0.0;}
        if(next_pos[2] <  0.0 && std::abs(next_pos[2]) < rel_tol * zlim) {next_pos[2] = 0.0;}

        if(next_pos[0] > xlim && next_pos[0] - xlim < rel_tol * xlim) {next_pos[0] = xlim;}
        if(next_pos[1] > ylim && next_pos[1] - ylim < rel_tol * ylim) {next_pos[1] = ylim;}
        if(next_pos[2] > zlim && next_pos[2] - zlim < rel_tol * zlim) {next_pos[2] = zlim;}

        if(!this->is_inside_of_boundary(next_pos))
        {
            std::cerr << "travel: initial  pos               = " << p        << " on " << f          << std::endl;
            std::cerr << "travel: initial face               = " << f << " -> " << this->triangle_at(f)             << std::endl;
            std::cerr << "travel: initial disp               = " << disp                          << std::endl;
            std::cerr << "travel: next face (unfolded)       = " << next.second                      << std::endl;
            std::cerr << "travel: next barycnetric (on edge) = " << on_edge_b                        << std::endl;
            std::cerr << "travel: next face                  = " << this->triangle_at(next.first)    << std::endl;
            std::cerr << "travel: next  pos                  = " << next_pos << " on " << next.first << std::endl;
            throw std::runtime_error("out of bound");
        }
    }

    return this->travel(std::make_pair(next_pos, next.first),
        rotate(tilt_angle_at(fd.edges[cs.first]),     // rotate disp by tilt_angle
               edge_over * (1.0 / length(edge_over)), // around the edge
               disp * (1 - cs.second)));              // the rest of displacement
}

std::pair<Real3, Polygon::FaceID> Polygon::travel(
        const std::pair<Real3, FaceID>& pos, const Real3& disp,
        const std::size_t restraint) const
{
    if(restraint == 0)
    {
        std::cerr << "movement along surface of a Polygon: "
                     "restraint hits 0. tolerance violated!" << std::endl;
        std::cerr << "the rest of displacement: " << disp    << std::endl;
        return pos;
    }
    const Real3& p = pos.first;
    const FaceID f = pos.second;
    const Real3 np = p + disp;

    const face_data& fd = this->face_at(f);
    const Triangle& tri = fd.triangle;
    const Barycentric b2(to_barycentric(np, tri));

    // if pos + disp is inside of the current face, just return the sum.
    if(::ecell4::is_inside(b2))
    {
        // to avoid numerical error that make the particle goes outside of the
        // face that the particle belongs, use `to_absolute`.
        return std::make_pair(to_absolute(b2, tri), f);
    }

    const Barycentric b1(to_barycentric(p, tri));
    const Barycentric db(b2 - b1);
    const std::pair<std::size_t, Real> cs   = first_cross_edge(b1, db);
    const std::pair<FaceID, Triangle>& next = fd.neighbor_cw[cs.first].front();

    // if the position is inside of the adjacent face, return the position
    // reconstructed from unfolded-Barycentric by using folded Triangle.
    const Barycentric unfolded_b(to_barycentric(np, next.second));
    if(::ecell4::is_inside(unfolded_b, 1e-8))
    {
        Real3 nxt = to_absolute(
                force_put_inside(unfolded_b), this->triangle_at(next.first));

        if(!this->is_inside_of_boundary(nxt))
        {
            const Real xlim    = this->edge_length_[0];
            const Real ylim    = this->edge_length_[1];
            const Real zlim    = this->edge_length_[2];
            const Real rel_tol = relative_tolerance;

            // allow numerical error in relative tolerance
            if(nxt[0] <  0.0 && std::abs(nxt[0]) < rel_tol * xlim) {nxt[0] = 0.0;}
            if(nxt[1] <  0.0 && std::abs(nxt[1]) < rel_tol * ylim) {nxt[1] = 0.0;}
            if(nxt[2] <  0.0 && std::abs(nxt[2]) < rel_tol * zlim) {nxt[2] = 0.0;}

            if(nxt[0] > xlim && nxt[0] - xlim < rel_tol * xlim) {nxt[0] = xlim;}
            if(nxt[1] > ylim && nxt[1] - ylim < rel_tol * ylim) {nxt[1] = ylim;}
            if(nxt[2] > zlim && nxt[2] - zlim < rel_tol * zlim) {nxt[2] = zlim;}

            if(!this->is_inside_of_boundary(nxt))
            {
                std::cerr << "travel: initial  pos         = " << p   << " on " << f            << std::endl;
                std::cerr << "travel: initial face         = " << f << " -> " << this->triangle_at(f)          << std::endl;
                std::cerr << "travel: initial disp         = " << disp                          << std::endl;
                std::cerr << "travel: next face (unfolded) = " << next.second                   << std::endl;
                std::cerr << "travel: next barycentric crd = " << unfolded_b                    << std::endl;
                std::cerr << "travel: next bary (inside)   = " << force_put_inside(unfolded_b)  << std::endl;
                std::cerr << "travel: next face            = " << this->triangle_at(next.first) << std::endl;
                std::cerr << "travel: next  pos            = " << nxt << " on " << next.first   << std::endl;
                throw std::runtime_error("particle goes exceeds the boundary");
            }
        }
        return std::make_pair(nxt, next.first);
        // use folded (normal) Triangle, NOT next.second
    }

    // stride over not only the edge but adjacent face.
    // XXX to make it sure that `on_edge` should be on the edge under the PBC,
    //     to_absolute is used with the next triangle.
    const Barycentric on_edge_b(to_barycentric(p + disp * cs.second, next.second));
    const Real3 edge_over = direction_of(fd.edges[cs.first]);
    Real3 next_pos  = to_absolute(on_edge_b, this->triangle_at(next.first));

    if(!this->is_inside_of_boundary(next_pos))
    {
        const Real xlim    = this->edge_length_[0];
        const Real ylim    = this->edge_length_[1];
        const Real zlim    = this->edge_length_[2];
        const Real rel_tol = relative_tolerance;

        // allow numerical error in relative tolerance
        if(next_pos[0] <  0.0 && std::abs(next_pos[0]) < rel_tol * xlim) {next_pos[0] = 0.0;}
        if(next_pos[1] <  0.0 && std::abs(next_pos[1]) < rel_tol * ylim) {next_pos[1] = 0.0;}
        if(next_pos[2] <  0.0 && std::abs(next_pos[2]) < rel_tol * zlim) {next_pos[2] = 0.0;}

        if(next_pos[0] > xlim && next_pos[0] - xlim < rel_tol * xlim) {next_pos[0] = xlim;}
        if(next_pos[1] > ylim && next_pos[1] - ylim < rel_tol * ylim) {next_pos[1] = ylim;}
        if(next_pos[2] > zlim && next_pos[2] - zlim < rel_tol * zlim) {next_pos[2] = zlim;}

        if(!this->is_inside_of_boundary(next_pos))
        {
            std::cerr << "travel: initial  pos               = " << p        << " on " << f          << std::endl;
            std::cerr << "travel: initial face               = " << f << " -> " << this->triangle_at(f)             << std::endl;
            std::cerr << "travel: initial disp               = " << disp                          << std::endl;
            std::cerr << "travel: next face (unfolded)       = " << next.second                      << std::endl;
            std::cerr << "travel: next barycnetric (on edge) = " << on_edge_b                        << std::endl;
            std::cerr << "travel: next face                  = " << this->triangle_at(next.first)    << std::endl;
            std::cerr << "travel: next  pos                  = " << next_pos << " on " << next.first << std::endl;
            throw std::runtime_error("particle goes exceeds the boundary");
        }
    }
    return this->travel(std::make_pair(next_pos, next.first),
        rotate(tilt_angle_at(fd.edges[cs.first]),     // rotate disp by tilt_angle
               edge_over * (1.0 / length(edge_over)), // around the edge
               disp * (1 - cs.second)),               // the rest of displacement
        restraint - 1);
}

} // ecell4
