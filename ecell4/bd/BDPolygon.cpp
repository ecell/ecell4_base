#include "BDPolygon.hpp"
#include "rotate_vector.hpp"
#include <set>
#include <limits>

namespace ecell4
{

namespace bd
{

void BDPolygon::detect_connectivity()
{
    detect_shared_edges();
    detect_shared_vertices();
    return;
}

std::pair<bool, edge_index_type>
BDPolygon::is_connected(const face_id_type& lhs, const face_id_type& rhs) const
{
    if(edge_pairs_.find(traits::make_pair(lhs, 0))->second.first == rhs)
        return std::make_pair(true, 0);
    if(edge_pairs_.find(traits::make_pair(lhs, 1))->second.first == rhs)
        return std::make_pair(true, 1);
    if(edge_pairs_.find(traits::make_pair(lhs, 2))->second.first == rhs)
        return std::make_pair(true, 2);
    return std::make_pair(false, 3);
}

std::pair<bool, vertex_index_type>
BDPolygon::is_share_vertex(const face_id_type& lhs, const face_id_type& rhs) const
{
    for(vertex_index_type i=0; i<3; ++i)
    {
        vertex_id_list const& list =
            vertex_groups_.find(std::make_pair(lhs, i));

        const face_finder cmp(rhs);
        if(std::find_if(list.begin(), list.end(), cmp) != list.end())
            return std::make_pair(true, i);
    }
    return std::make_pair(false, 3);
}

Real BDPolygon::distance(const std::pair<Real3, face_id_type>& lhs,
                         const std::pair<Real3, face_id_type>& rhs) const
{
    // CASE 0: on the same face
    if(lhs.second == rhs.second)
        return length(lhs.first - rhs.first);

    // CASE 1: lhs and rhs share an edge
    { 
    const std::pair<bool, edge_index_type> edg =
        this->is_connected(lhs.second, rhs.second);

    if(edg.first) // connected
    {
        // make two faces being same plane and then calculate distance
        const Triangle& lhs_t = faces_[lhs.second];
        const Triangle& rhs_t = faces_[rhs.second];
        const Real ang = angle(lhs_t.normal(), rhs_t.normal());

        const Real3 developped = lhs_t.vertex_at(edg.second) + 
            rotate(-ang,
                   rhs.first - lhs_t.vertex_at(edg.second),
                   lhs_t.edge_at(edg.second));
        return length(lhs.first - developped);
    }
    }// end CASE 1

    // CASE 2: lhs and rhs share a vertex
    {
    const std::pair<bool, vertex_index_type> vtx =
        this->is_share_vertex(lhs.second, rhs.second);

    if(vtx.first) // XXX!
    {
        Real whole_angle = 0.; // whole angle around the vertex
        Real inter_angle = 0.; // angle lhs--vtx--rhs

        Real r0 = 0.; // distance from the shared vertex to lhs.position
        Real r1 = 0.; // ...                                rhs.position

        vertex_id_list const& list = const_at(vertex_groups_,
                traits::make_vertex_id(lhs.second, vtx.second));

        bool inter = false;
        for(vertex_id_list::const_iterator
                iter = list.begin(); iter != list.end(); ++iter) // ccw
        {
            const Triangle& f = this->faces_.at(traits::get_face_id(*iter));
            whole_angle += f.angle_at(traits::get_vertex_index(*iter));

            if(!inter && traits::get_face_id(*iter) == lhs.second)
            {
                inter = true;
                angle_sum += angle(f.vertex_at(iter->second) - lhs.first,
                                   f.edge_at(iter->second==0?2:iter->second-1));
            }
            else if(inter && traits::get_face_id(*iter) == rhs.second)
            {
                inter = false;
                angle_sum += angle(rhs.first - f.vertex_at(iter->second),
                                   f.edge_at(iter->second));
            }
            else if(inter)
            {
                angle_sum += this->faces_.at(iter->first).angle_at(iter->second);
            }
        }
        if(inter)
        {
            for(vertex_id_list::const_iterator iter = list.begin();
                    iter != list.end(); ++iter)
            {
                if(iter->first == rhs.second)
                {
                    const Triangle& f = this->faces_.at(iter->first);
                    angle_sum += angle(rhs.first - f.vertex_at(iter->second),
                                       f.edge_at(iter->second));
                    break;
                }
                angle_sum += this->faces_.at(iter->first).angle_at(iter->second);
            }
        }
        const Real min_angle = std::min(angle_sum, whole_angle - angle_sum);
        return r0 * r0 + r1 * r1 - 2. * r0 * r1 * std::cos(min_angle);
    }
    }// end CASE 2

    // lhs and rhs don't share edge nor vertex
    return std::numeric_limits<Real>::infinity();
}

std::pair<std::pair<Real3, BDPolygon::face_id_type>, Real3>
BDPolygon::move_next_face(const std::pair<Real3, face_id_type>& pos,
                          const Real3& disp) const
{
    const face_type& f = faces_[pos.second];

    const boost::array<Real, 3> newpos =
        this->to_barycentric(pos.first + disp, f);

    if(is_inside(newpos))
    {
        return std::make_pair(
                std::make_pair(this->to_absolute(newpos, f), pos.second),
                Real3(0.,0.,0.));
    }

#ifdef ECELL4_TWOD_BD_DUMP
    const boost::array<Real, 3> curpos =
        this->to_barycentric(pos.first, f);
    std::cerr << "curpos       = " << curpos[0] << ", " << curpos[1] << ", " << curpos[2] << std::endl;
    std::cerr << "newpos       = " << newpos[0] << ", " << newpos[1] << ", " << newpos[2] << std::endl;
    std::cerr << "length(disp) = " << length(disp) << std::endl;
    std::cerr << "dot(norm, nd)= " << dot_product(disp, f.normal()) << std::endl;
#endif

    const boost::array<Real, 3> bpos = to_barycentric(pos.first, f);
    boost::array<Real, 3> bdis;
    bdis[0] = newpos[0] - bpos[0];
    bdis[1] = newpos[1] - bpos[1];
    bdis[2] = newpos[2] - bpos[2];

#ifdef ECELL4_TWOD_BD_DUMP
    std::cerr << "pos(bary)    = " << bpos[0] << ", " << bpos[1] << ", " << bpos[2] << std::endl;
    std::cerr << "disp(bary)   = " << bdis[0] << ", " << bdis[1] << ", " << bdis[2] << std::endl;
#endif

    const std::pair<uint32_t, Real> cross = crossed_edge(bpos, bdis);
    const edge_finder cmp(std::make_pair(pos.second, cross.first));
    const face_id_type next_face_id =
        std::find_if(edge_pairs_.begin(), edge_pairs_.end(), cmp)->second.first;
    const face_type& next_face = faces_[next_face_id];

#ifdef ECELL4_TWOD_BD_DUMP
    std::cerr << "cross ratio  = " << cross.second << std::endl;
    std::cerr << "next_face id = " << next_face_id << std::endl;
#endif

    // turn displacement
    const Real ang = angle(f.normal(), next_face.normal());
    const Real3 axis = f.edge_at(cross.first) / length(f.edge_at(cross.first));
    const Real3 next_disp = rotate(ang, axis, disp) * (1. - cross.second);

#ifdef ECELL4_TWOD_BD_DUMP
    std::cerr << "next disp    = " << next_disp << std::endl;
    std::cerr << "next disp len= " << length(next_disp) << std::endl;
    std::cerr << "dot(norm, nd)= " << dot_product(next_disp, next_face.normal()) << std::endl;
    std::cerr << std::endl;
#endif

    return std::make_pair(std::make_pair(pos.first + disp * cross.second, next_face_id),
            next_disp);
}

//!
//  find first-crossing edge and
//  return its index and displacement ratio to the cross section
std::pair<edge_index_type, Real>
BDPolygon::crossed_edge(const barycentric_type& pos, const barycentric_type& disp) const
{
    barycentric_type npos = pos + disp;

    if(npos[0] * npos[1] * npos[2] < 0.) // (+, +, -) or one of its permutations
    {
             if(npos[0] < 0.) return std::make_pair(1, cross_section(pos, disp, 1));
        else if(npos[1] < 0.) return std::make_pair(2, cross_section(pos, disp, 2));
        else if(npos[2] < 0.) return std::make_pair(0, cross_section(pos, disp, 0));
        else throw std::logic_error("Polygon::crossed_edge: never reach here");
    }
    else // (+, -, -) or one of its permutations
    {
        assert(!is_inside(npos)); // not (+, +, +) case
//         assert(on_plane(npos));   // summation is 1
        if(npos[0] > 0.)
        {
            const Real ab = cross_section(pos, disp, 0);
            const Real ca = cross_section(pos, disp, 2);
            return (ab > ca) ? std::make_pair(2, ca) : std::make_pair(0, ab);
        }
        else if(npos[1] > 0.)
        {
            const Real ab = cross_section(pos, disp, 0);
            const Real bc = cross_section(pos, disp, 1);
            return (bc > ab) ? std::make_pair(0, ab) : std::make_pair(1, bc);
        }
        else if(npos[2] > 0.)
        {
            const Real bc = cross_section(pos, disp, 1);
            const Real ca = cross_section(pos, disp, 2);
            return (ca > bc) ? std::make_pair(1, bc) : std::make_pair(2, ca);
        }
    }
}

Real BDPolygon::cross_section(const barycentric_type& pos,
        const barycentric_type& disp, const edge_index_type e) const
{
    switch(e)
    {
        case 0:  return -pos[2] / disp[2];
        case 1:  return -pos[0] / disp[0];
        case 2:  return -pos[1] / disp[1];
        default: return std::numeric_limits<Real>::infinity();
    }
}


void BDPolygon::detect_shared_vertices()
{
    std::set<vertex_id_type> is_detected;
    const Real same_position_tolerance = 1e-6;
    for(std::size_t fidx = 0; fidx < faces_.size(); ++fidx)
    {
        for(std::size_t vidx = 0; vidx < 3; ++vidx)
        {
            const vertex_id_type current_vtx = std::make_pair(fidx, vidx);

            if(is_detected.count(current_vtx) == 1) continue;
            vertex_id_list sharing_list;
            sharing_list.push_back(current_vtx);
            is_detected.insert(current_vtx);

            edge_id_type lookup = current_vtx;
            std::size_t i=0;
            for(; i<100; ++i)
            {
                const std::size_t f = lookup.first;
                const std::size_t e = (lookup.second == 0) ? 2: lookup.second - 1;
                const edge_id_type eid  = std::make_pair(f, e);
                const edge_id_type next = this->edge_pairs_[eid];
                lookup = next;
                if(lookup.first == current_vtx.first) break;
                sharing_list.push_back(lookup);
            }
            if(i==100) throw std::logic_error("too many sharing vertex");

            for(std::size_t i=0; i<sharing_list.size(); ++i)
                vertex_groups_[sharing_list.at(i)] = sharing_list;
        }
    }

    return;
}


void BDPolygon::detect_shared_edges()
{
    std::set<edge_id_type> is_detected;
    const Real same_position_tolerance = 1e-6;

    for(std::size_t fidx = 0; fidx < faces_.size(); ++fidx)
    {
        const std::size_t currentf = fidx;
        for(std::size_t eidx = 0; eidx < 3; ++eidx)
        {
            const edge_id_type current_edge = std::make_pair(currentf, eidx);
            if(is_detected.count(current_edge) == 1) continue;

            const std::size_t start_v = eidx;
            const std::size_t end_v   = (eidx==2) ? 0 : eidx+1;
            const Real3 start_pos = faces_.at(currentf).vertex_at(start_v);
            const Real3 end_pos   = faces_.at(currentf).vertex_at(end_v);

            const Real start_dist = length(start_pos);
            const Real start_length_scale =
                (start_dist >= same_position_tolerance) ?
                (start_dist) : (same_position_tolerance);

            const Real end_dist = length(end_pos);
            const Real end_length_scale =
                (end_dist >= same_position_tolerance) ?
                (end_dist) : (same_position_tolerance);

            bool found = false;
            for(std::size_t f = currentf + 1; f < faces_.size(); ++f)
            {
                for(std::size_t e = 0; e < 3; ++e)
                {
                    const Real start_pos_dist =
                        length(start_pos - faces_.at(f).vertex_at(e));
                    if(start_pos_dist >
                            start_length_scale * same_position_tolerance)
                        continue;

                    const Real end_pos_dist =
                        length(end_pos - faces_.at(f).vertex_at((e==0)?2:e-1));
                    if(end_pos_dist <= 
                            end_length_scale * same_position_tolerance)
                    {
                        const edge_id_type detected = std::make_pair(f, (e==0)?2:e-1);
                        this->edge_pairs_[current_edge] = detected;
                        this->edge_pairs_[detected] = current_edge;
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
 

}// bd
}// ecell
