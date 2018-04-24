#include <ecell4/core/Polygon.hpp>
#include <boost/container/flat_map.hpp>

namespace ecell4
{

const Real Polygon::absolute_tolerance = 1e-12;
const Real Polygon::relative_tolerance = 1e-8;

void Polygon::assign(const std::vector<Triangle>& ts)
{
    const Real pi = boost::math::constants::pi<Real>();
    const Real tol_abs2 = absolute_tolerance * absolute_tolerance;
    const Real tol_rel2 = relative_tolerance * relative_tolerance;

    vertices_.clear();
       faces_.clear();
       edges_.clear();
    this->total_area_ = 0.0;

    // prepair temporal data storage
    typedef std::pair<FaceID, std::size_t>          fid_vidx_pair;
    typedef std::pair<Real3, std::vector<fid_vidx_pair> > tmp_vtx_type;
    typedef boost::container::flat_map<VertexID, tmp_vtx_type> tmp_vertex_map;
    tmp_vertex_map tmp_vtxs;

    // first, generate (FaceIDs for all triangles) and (EdgeIDs for all Edges).
    // and collect vertices that are at the same position.
    for(std::vector<Triangle>::const_iterator
            t_iter(ts.begin()), t_end(ts.end()); t_iter != t_end; ++t_iter)
    {
        const Triangle& triangle = *t_iter;
        this->total_area_ += triangle.area();

        const FaceID fid(faces_.size());
        face_data fd;
        fd.triangle = triangle;

        for(std::size_t i=0; i<3; ++i)
        {
            const Real3& v1 = triangle.vertices()[i];
            boost::optional<VertexID> found_vtx = boost::none;

            // find near vertex
            for(tmp_vertex_map::iterator
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
                const VertexID new_vid(tmp_vtxs.size());
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
            const EdgeID eid(edges_.size());
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
    // * set outgoing_edges without order
    for(tmp_vertex_map::const_iterator
            vi(tmp_vtxs.begin()), ve(tmp_vtxs.end()); vi != ve; ++vi)
    {
        const VertexID vid = vi->first;
        const Real3          pos = vi->second.first;
        const std::vector<fid_vidx_pair>& face_pos = vi->second.second;

        vertex_data vd;
        vd.position = pos;

        // * set vertex.outgoing_edges, but not sorted.
        for(std::vector<fid_vidx_pair>::const_iterator
                i(face_pos.begin()), e(face_pos.end()); i!=e; ++i)
        {
            const FaceID fid = i->first;
            const std::size_t  idx = i->second;
            face_data& fd = this->faces_.at(fid);

            assert(vid == fd.vertices[idx]);
            vd.outgoing_edges.push_back(std::make_pair(fd.edges[idx], 0.0));
        }
        this->vertices_.push_back(vd);
    }

    // * refine vertex positions
    for(face_container_type::iterator
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
    for(face_container_type::const_iterator
            fi(this->faces_.begin()), fe(this->faces_.end()); fi != fe; ++fi)
    {
        const face_data& fd = *fi;
        for(std::size_t i=0; i<3; ++i)
        {
            const EdgeID eid = fd.edges[i];
            this->edges_.at(eid).length    = fd.triangle.length_of_edge_at(i);
            this->edges_.at(eid).direction = fd.triangle.edge_at(i);
        }
    }

    // search pairs of opposite edges & calculate edge.tilt.
    for(std::size_t i=0; i<edges_.size(); ++i)
    {
        const EdgeID   eid(i);
        const VertexID start  = this->target_of(eid);
        const VertexID target = this->target_of(
                this->next_of(this->next_of(eid)));

        bool opposite_found = false;
        const std::vector<std::pair<EdgeID, Real> >& vd =
            this->vertices_.at(start).outgoing_edges;

        for(std::vector<std::pair<EdgeID, Real> >::const_iterator
                iter(vd.begin()), iend(vd.end()); iter != iend; ++iter)
        {
            const EdgeID outgoing = iter->first;
            if(this->target_of(outgoing) == target)
            {
                // found opposite edge! calculate tilt...
                this->edge_at(eid).opposite_edge = outgoing;

                const FaceID fid1 = face_of(eid);
                const FaceID fid2 = face_of(outgoing);
                const Real3 n1 = this->faces_.at(fid1).triangle.normal();
                const Real3 n2 = this->faces_.at(fid2).triangle.normal();
                const Real3 cr = cross_product(this->edge_at(eid).direction, n1);
                const Real  sg = dot_product(cr, n2);
                const Real ang = calc_angle(n1, n2) * (sg > 0 ? 1 : -1);

                this->edges_.at(i       ).tilt = ang;
                this->edges_.at(outgoing).tilt = ang;

                opposite_found = true;
                break;
            }
        }
//         if(!opposite_found)
//         {
//             std::cerr << "edge " << this->edge_at(eid).direction << std::endl;
//             std::cerr << "towards: "
//                 << this->vertex_at(this->edge_at(eid).target).position << std::endl;
//             std::cerr << "does not have any opposite pair." << std::endl;
//         }
        assert(opposite_found);
    }

    // set vertex_data.angle by traversing edges.
    for(std::size_t i=0; i<vertices_.size(); ++i)
    {
        vertex_data& vtx = this->vertices_[i];

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
                const std::vector<EdgeID>::iterator found = std::find(
                    outgoing_edges_tmp.begin(), outgoing_edges_tmp.end(), current);
                outgoing_edges_tmp.erase(found);
            }

            const face_data& f = this->faces_.at(this->face_of(current));
            Real angle = std::numeric_limits<Real>::infinity();
            for(std::size_t idx=0; idx<3; ++idx)
            {
                if(f.edges[idx] == current)
                {
                    angle = f.triangle.angle_at(idx);
                    break;
                }
            }
            assert(angle != std::numeric_limits<Real>::infinity());

            total_angle += angle;
            vtx.outgoing_edges.push_back(std::make_pair(current, angle));
            current = this->opposite_of(this->next_of(this->next_of(current)));
        }
        while(current != start);

        this->vertices_[i].apex_angle = total_angle;

        assert(outgoing_edges_tmp.empty());
        assert(vtx.outgoing_edges.size() == num_edges);
    }

    // make neighbor list for faces!

    for(std::vector<face_data>::iterator
            fi(this->faces_.begin()), fe(this->faces_.end()); fi != fe; ++fi)
    {
        face_data& face = *fi;
        for(std::size_t i=0; i<3; ++i)
        {
            const VertexID vid = face.vertices[i];
            const Real3 v_pos        = face.triangle.vertex_at(i);
            const Real  offset_angle = face.triangle.angle_at(i);
            const Real3 normal       = face.triangle.normal();

            // counter clock wise
            {
                const Real3 ref_edge = face.triangle.edge_at(i==0?2:i-1) /
                    (-length(face.triangle.edge_at(i==0?2:i-1)));

                EdgeID current_edge  = face.edges[i];
                Real         current_angle = 0.0;
                do
                {
                    current_edge = opposite_of(next_of(next_of(current_edge)));
                    const FaceID fid  = face_of(current_edge);

                    const std::size_t vidx0 = this->face_at(fid).index_of(vid);
                    const std::size_t vidx1 = this->face_at(fid).index_of(
                                              target_of(current_edge));
                    const std::size_t vidx2 = this->face_at(fid).index_of(
                                              target_of(next_of(current_edge)));

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
                    face.neighbor_ccw[i].push_back(
                            std::make_pair(fid, Triangle(unfolded)));

                    current_angle = next_angle;
                }
                while(current_angle + offset_angle <= pi);
            }

            // clock wise
            {
                const Real3 ref_edge = face.triangle.edge_at(i) /
                    length(face.triangle.edge_at(i));

                EdgeID current_edge  = face.edges[i];
                Real         current_angle = 0.0;
                do
                {
                    current_edge  = next_of(opposite_of(current_edge));
                    const FaceID fid  = face_of(current_edge);

                    const std::size_t vidx0 = this->face_at(fid).index_of(vid);
                    const std::size_t vidx1 = this->face_at(fid).index_of(
                                              target_of(current_edge));
                    const std::size_t vidx2 = this->face_at(fid).index_of(
                                              target_of(next_of(current_edge)));

                    const Real prev_angle = current_angle;
                    current_angle += this->face_at(fid).triangle.angle_at(vidx0);

                    boost::array<Real3, 3> unfolded;
                    unfolded[vidx0] = v_pos;
                    unfolded[vidx1] = v_pos +
                        rotate(-current_angle, normal, ref_edge) *
                        length_of(current_edge);
                    unfolded[vidx2] = v_pos +
                        rotate(-prev_angle,    normal, ref_edge) *
                        length_of(next_of(next_of(current_edge)));

                    face.neighbors.push_back(fid);
                    face.neighbor_cw[i].push_back(
                            std::make_pair(fid, Triangle(unfolded)));
                }
                while(current_angle + offset_angle <= pi);
            }
        }
        std::sort(face.neighbors.begin(), face.neighbors.end());
        const std::vector<FaceID>::iterator last =
            std::unique(face.neighbors.begin(), face.neighbors.end());
        face.neighbors.erase(last, face.neighbors.end());
    }
    return;
}

Real Polygon::distance_sq(
        const std::pair<Real3, FaceID>& pos1,
        const std::pair<Real3, FaceID>& pos2) const
{
    typedef utils::pair_first_element_unary_predicator<FaceID, Triangle>
            face_finder_type;

    const Real3& p1 = pos1.first;
    const FaceID f1 = pos1.second;
    const FaceID f2 = pos2.second;
    const Barycentric  b2 = to_barycentric(pos2.first, face_at(f2).triangle);

    const face_data& face = face_at(f1);
    const Real3&   normal = face.triangle.normal();

    Real solution = std::numeric_limits<Real>::infinity();
    for(std::size_t i=0; i<3; ++i)
    {
        const VertexID vid = face.vertices[i];
        const Real3& vpos(position_at(vid));
        const Real3 vtop1(p1 - vpos);

        { // counter clock wise
            const std::vector<std::pair<FaceID, Triangle> >::const_iterator fi =
                std::find_if(face.neighbor_ccw[i].begin(),
                             face.neighbor_ccw[i].end(), face_finder_type(f2));
            if(fi != face.neighbor_ccw[i].end())
            {
                // unfolded place of p2
                const Real3 p2 = to_absolute(b2, fi->second);
                const Real3 vtop2(p2 - vpos);
                // check the angle between p1-v-p2 does not exceeds PI
                if(dot_product(normal, cross_product(vtop1, vtop2)) >= 0)
                {
                    solution = std::min(length_sq(p1 - p2), solution);
                }
            }
        }
        { // clock wise
            const std::vector<std::pair<FaceID, Triangle> >::const_iterator fi =
                std::find_if(face.neighbor_cw[i].begin(),
                             face.neighbor_cw[i].end(), face_finder_type(f2));
            if(fi != face.neighbor_cw[i].end())
            {
                // unfolded place of p2
                const Real3 p2 = to_absolute(b2, fi->second);
                const Real3 vtop2(p2 - vpos);
                // check the angle between p1-v-p2 does not exceeds PI
                if(dot_product(normal, cross_product(vtop1, vtop2)) <= 0)
                {
                    solution = std::min(length_sq(p1 - p2), solution);
                }
                break;
            }
        }
    }
    if(solution != std::numeric_limits<Real>::infinity())
    {
        return solution;
    }

    boost::optional<VertexID> connected = boost::none;
    // search f2 in the connected faces (if the apex angle of the vertex
    // exceeded 2PI, the minimum path can be the path that goes through
    // the vertex).
    for(std::size_t i=0; i<3; ++i)
    {
        const VertexID vid = face.vertices[i];
        const std::vector<std::pair<EdgeID, Real> >&
            oes = this->vertex_at(vid).outgoing_edges;
        for(std::vector<std::pair<EdgeID, Real> >::const_iterator
                iter(oes.begin()), iend(oes.end()); iter!=iend; ++iter)
        {
            if(face_of(iter->first) == f2)
            {
                assert(!connected);
                connected = vid;
            }
        }
    }

    if(connected)
    {
        const Real3& vpos = position_at(*connected);
        const Real   lsq1 = length_sq(p1         - vpos);
        const Real   lsq2 = length_sq(pos2.first - vpos);
        // (x+y)^2 = x^2 + y^2 + 2xy = x^2 + y^2 + 2 * sqrt(x^2y^2)
        return lsq1 + lsq2 + 2 * std::sqrt(lsq1 * lsq2);
    }
    return std::numeric_limits<Real>::infinity();
}

Real3 Polygon::direction(
        const std::pair<Real3, FaceID>& pos1,
        const std::pair<Real3, FaceID>& pos2) const
{
    typedef utils::pair_first_element_unary_predicator<FaceID, Triangle>
            face_finder_type;

    const Real3& p1 = pos1.first;
    const FaceID f1 = pos1.second;
    const FaceID f2 = pos2.second;
    const Barycentric b2 = to_barycentric(pos2.first, face_at(f2).triangle);

    const face_data& face = face_at(f1);
    const Real3&   normal = face.triangle.normal();

    Real mindist2 = std::numeric_limits<Real>::infinity();
    Real3 direction(0,0,0);
    for(std::size_t i=0; i<3; ++i)
    {
        const VertexID vid = face.vertices[i];
        const Real3& vpos(position_at(vid));
        const Real3 vtop1(p1 - vpos);

        { // counter clock wise
            const std::vector<std::pair<FaceID, Triangle> >::const_iterator fi =
                std::find_if(face.neighbor_ccw[i].begin(),
                             face.neighbor_ccw[i].end(), face_finder_type(f2));
            if(fi != face.neighbor_ccw[i].end())
            {
                // unfolded place of p2
                const Real3 p2 = to_absolute(b2, fi->second);
                const Real3 vtop2(p2 - vpos);
                // check the angle between p1-v-p2 does not exceeds PI
                if(dot_product(normal, cross_product(vtop1, vtop2)) >= 0)
                {
                    const Real3 dr = p2 - p1;
                    const Real  d2 = length_sq(dr);
                    if(d2 < mindist2)
                    {
                        mindist2  = d2;
                        direction = dr;
                    }
                }
                break;
            }
        }
        { // clock wise
            const std::vector<std::pair<FaceID, Triangle> >::const_iterator fi =
                std::find_if(face.neighbor_cw[i].begin(),
                             face.neighbor_cw[i].end(), face_finder_type(f2));
            if(fi != face.neighbor_cw[i].end())
            {
                // unfolded place of p2
                const Real3 p2 = to_absolute(b2, fi->second);
                const Real3 vtop2(p2 - vpos);
                // check the angle between p1-v-p2 does not exceeds PI
                if(dot_product(normal, cross_product(vtop1, vtop2)) <= 0)
                {
                    const Real3 dr = p2 - p1;
                    const Real  d2 = length_sq(dr);
                    if(d2 < mindist2)
                    {
                        mindist2  = d2;
                        direction = dr;
                    }
                }
                break;
            }
        }
    }
    if(mindist2 == std::numeric_limits<Real>::infinity())
    {
        throw std::runtime_error("polygon::direction: couldn't find the min path");
    }
    return direction;
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
        return std::make_pair(np, f);
    }

    const Barycentric b1(to_barycentric(p, tri));
    const Barycentric db(b2 - b1);
    const std::pair<std::size_t, Real> cs   = first_cross_edge(b1, db);
    const std::pair<FaceID, Triangle>& next = fd.neighbor_cw[cs.first].front();

    // if the position is inside of the adjacent face, return the position
    // reconstructed from unfolded-Barycentric by using folded Triangle.
    const Barycentric unfolded_b(to_barycentric(np, next.second));
    if(::ecell4::is_inside(unfolded_b))
    {
        return std::make_pair(// use folded (normal) Triangle, NOT next.second
            to_absolute(unfolded_b, this->triangle_at(next.first)), next.first);
    }

    // stride over not only the edge but adjacent face.
    // XXX to make it sure that `on_edge` should be on the edge under the PBC,
    //     to_absolute is used with the next triangle.
    const Barycentric on_edge_b(to_barycentric(p + disp * cs.second, next.second));
    return this->travel(std::make_pair(
        to_absolute(on_edge_b, this->triangle_at(next.first)), next.first),
        rotate(tilt_angle_at(fd.edges[cs.first]), // rotate disp by tilt_angle
               direction_of(fd.edges[cs.first]),  // around the edge
               disp * (1 - cs.second)));          // the rest of displacement
}

} // ecell4
