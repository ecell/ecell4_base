#include "BDPolygon.hpp"
#include "rotate_vector.hpp"
#include <set>
#include <limits>

namespace ecell4
{

namespace bd
{

void BDPolygon::detect_connectivity(const Real tolerance)
{
    detect_shared_edges(tolerance);
    detect_shared_vertices(tolerance);
    return;
}

std::pair<bool, BDPolygon::edge_index_type>
BDPolygon::is_connected(const face_id_type& lhs, const face_id_type& rhs) const
{
    if(traits::get_face_id(const_at(edge_pairs_, traits::make_edge_id(lhs, 0)))
        == rhs)
        return std::make_pair(true, 0);
    if(traits::get_face_id(const_at(edge_pairs_, traits::make_edge_id(lhs, 1)))
        == rhs)
        return std::make_pair(true, 1);
    if(traits::get_face_id(const_at(edge_pairs_, traits::make_edge_id(lhs, 2)))
        == rhs)
        return std::make_pair(true, 2);

    return std::make_pair(false, 3);
}

std::pair<bool, BDPolygon::vertex_index_type>
BDPolygon::is_share_vertex(const face_id_type& lhs, const face_id_type& rhs) const
{
    for(vertex_index_type i=0; i<3; ++i)
    {
        vertex_id_list const& list =
            const_at(vertex_groups_, traits::make_vertex_id(lhs, i)).first;

        typedef vertex_id_list::const_iterator iter_type;
        for(iter_type iter = list.begin(); iter != list.end(); ++iter)
            if(traits::get_face_id(*iter) == rhs)
                return std::make_pair(true, i);
    }
    return std::make_pair(false, 3);
}

Real BDPolygon::distance_sq(const std::pair<Real3, face_id_type>& lhs,
                            const std::pair<Real3, face_id_type>& rhs) const
{
    // CASE 0: on the same face
    if(lhs.second == rhs.second)
        return length_sq(lhs.first - rhs.first);

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
        return length_sq(lhs.first - developped);
    }
    }// end CASE 1

    // CASE 2: lhs and rhs share a vertex
    {
    const std::pair<bool, vertex_index_type> vtx =
        this->is_share_vertex(lhs.second, rhs.second);

    if(vtx.first)
    {
        const Real3 vtx_position = faces_.at(lhs.second).vertex_at(vtx.second);
        const Real3 lhs_to_vtx = vtx_position - lhs.first;
        const Real3 vtx_to_rhs = rhs.first - vtx_position;
        const Real lhs_to_vtx_lensq = length_sq(lhs_to_vtx);
        const Real rhs_to_vtx_lensq = length_sq(vtx_to_rhs);

        vertex_list_angle_pair const& vlap = const_at(vertex_groups_,
                traits::make_vertex_id(lhs.second, vtx.second));
        const Real whole_angle = vlap.second;
        Real inter_angle = angle(lhs_to_vtx, faces_.at(lhs.second).edge_at(
                    (vtx.second == 0) ? 2 : vtx.second - 1));

        vertex_id_list const& vlist = vlap.first;
        vertex_id_list::const_iterator iter = vlist.begin();
        ++iter;// the first element is same as lhs.second + vtx.second
        while(iter != vlist.end())
        {
            const face_id_type fid = traits::get_face_id(*iter);
            const Triangle&      f = faces_.at(fid);
            if(fid == rhs.second)
            {
                inter_angle += angle(vtx_to_rhs,
                        f.edge_at(traits::get_edge_index(traits::vtoe(*iter))));
                break;
            }
            const vertex_index_type vid = traits::get_vertex_index(*iter);
            inter_angle += f.angle_at(vid);
            ++iter;
        }
//         assert(inter_angle <= whole_angle);
        if(inter_angle > whole_angle)
        {
            std::cerr << "inter particle angle exceeds the whole angle" << std::endl;
            std::cerr << "p1 is on " << lhs.second << "-th face" << std::endl;
            std::cerr << "p2 is on " << rhs.second << "-th face" << std::endl;
            std::cerr << "inter angle " << inter_angle << std::endl;
            std::cerr << "whole angle " << whole_angle << std::endl;
            std::cerr << "faces sharing the vertex ";
            for(vertex_id_list::const_iterator iter = vlist.begin(); iter != vlist.end(); ++iter)
            {
                std::cerr << traits::get_face_id(*iter) << ", ";
            }
            std::cerr << std::endl;
        }

        const Real min_angle = std::min(inter_angle, whole_angle - inter_angle);
        return lhs_to_vtx_lensq + rhs_to_vtx_lensq - 2. *
               std::sqrt(lhs_to_vtx_lensq * rhs_to_vtx_lensq) * std::cos(min_angle);
    }
    }// end CASE 2

    // lhs and rhs don't share edge nor vertex
    return std::numeric_limits<Real>::infinity();
}

Real3 BDPolygon::inter_position_vector(
        const std::pair<Real3, face_id_type>& lhs,
        const std::pair<Real3, face_id_type>& rhs) const
{
    // CASE 0: on the same face
    if(lhs.second == rhs.second)
        return rhs.first - lhs.first;

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
        return lhs.first - developped;
    }
    }// end CASE 1

    // CASE 2: lhs and rhs share a vertex
    {
    const std::pair<bool, vertex_index_type> vtx =
        this->is_share_vertex(lhs.second, rhs.second);

    if(vtx.first)
    {
        const Real3 vtx_position = faces_.at(lhs.second).vertex_at(vtx.second);
        const Real3 lhs_to_vtx(vtx_position - lhs.first);
        const Real3 vtx_to_rhs(rhs.first - vtx_position);
        const Real lhs_to_vtx_lensq = length_sq(lhs_to_vtx);
        const Real rhs_to_vtx_lensq = length_sq(vtx_to_rhs);
        const Real3 normal(faces_.at(lhs.second).normal());

        vertex_list_angle_pair const& vlap = const_at(vertex_groups_,
                traits::make_vertex_id(lhs.second, vtx.second));
        const Real whole_angle = vlap.second;
        Real inter_angle = angle(lhs_to_vtx, faces_.at(lhs.second).edge_at(
                    (vtx.second == 0) ? 2 : vtx.second - 1));

        vertex_id_list const& vlist = vlap.first;
        vertex_id_list::const_iterator iter = vlist.begin();
        ++iter;
        while(iter != vlist.end())
        {
            const face_id_type fid = traits::get_face_id(*iter);
            const Triangle&      f = faces_.at(fid);
            if(fid == rhs.second)
            {
                inter_angle += angle(vtx_to_rhs,
                        f.edge_at(traits::get_edge_index(traits::vtoe(*iter))));
                break;
            }
            const vertex_index_type vid = traits::get_vertex_index(*iter);
            inter_angle += f.angle_at(vid);
            ++iter;
        }
//         assert(inter_angle <= whole_angle);
        if(inter_angle > whole_angle)
        {
            std::cerr << "inter particle angle exceeds the whole angle" << std::endl;
            std::cerr << "p1 is on " << lhs.second << "-th face" << std::endl;
            std::cerr << "p2 is on " << rhs.second << "-th face" << std::endl;
            std::cerr << "inter angle " << inter_angle << std::endl;
            std::cerr << "whole angle " << whole_angle << std::endl;
            std::cerr << "faces sharing the vertex ";
            for(vertex_id_list::const_iterator iter = vlist.begin(); iter != vlist.end(); ++iter)
            {
                std::cerr << traits::get_face_id(*iter) << ", ";
            }
            std::cerr << std::endl;
        }
        // XXX:NOTE
        // in the case of inter_angle == whole_angle - inter_angle, it is good
        // to choose inter-position-vector randomly. but here it is unidirectional.
        const Real min_angle = std::min(inter_angle, whole_angle - inter_angle);

        return lhs_to_vtx + rotate(min_angle, normal, lhs_to_vtx * (-1.0)) *
                         std::sqrt(rhs_to_vtx_lensq / lhs_to_vtx_lensq);
    }
    }// end CASE 2

    throw NotSupported("Polygon::inter_position_vector couldn't determine the path");
}

std::pair<std::pair<Real3, BDPolygon::face_id_type>, Real3>
BDPolygon::move_next_face(const std::pair<Real3, face_id_type>& pos,
                          const Real3& disp) const
{
    const face_type& f = faces_.at(pos.second);
    const barycentric_type newpos = to_barycentric(pos.first + disp, f);

    if(is_inside(newpos))
        return std::make_pair(
            std::make_pair(to_absolute(newpos, f), pos.second), Real3(0.,0.,0.));

    // calculate partial displacement to the edge
    const barycentric_type bpos = to_barycentric(pos.first, f);
    const barycentric_type bdis = newpos - bpos;

    const std::pair<edge_index_type, Real> cross = crossed_edge(bpos, bdis);

    const face_id_type next_face_id = traits::get_face_id(
        const_at(edge_pairs_, traits::make_edge_id(pos.second, cross.first)));
    const face_type& next_face = faces_.at(next_face_id);
    const Real3 next_pos  = pos.first + disp * cross.second;

    // turn displacement
    const Real  ang  = angle(f.normal(), next_face.normal());
    const Real3 axis = f.edge_at(cross.first) / f.length_of_edge_at(cross.first);
    const Real3 next_disp = rotate(ang, axis, disp) * (1. - cross.second);

    return std::make_pair(std::make_pair(next_pos, next_face_id), next_disp);
}

//!
//  find first-crossing edge and
//  return its index and displacement ratio to the cross section
std::pair<BDPolygon::edge_index_type, Real>
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
        if(npos[0] > 0. && npos[1] > 0. && npos[2] > 0) // not (+, +, +) case
            throw std::invalid_argument("BDPolygon::crossed_edge");

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
        else // if(npos[2] > 0.)
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

void BDPolygon::detect_shared_vertices(const Real tolerance)
{
    for(std::size_t fidx = 0; fidx < faces_.size(); ++fidx)
    {
        for(vertex_index_type vidx = 0; vidx < 3; ++vidx)
        {
            const vertex_id_type current_vtx = traits::make_vertex_id(fidx, vidx);

            Real whole_angle = faces_.at(fidx).angle_at(vidx);
            vertex_id_list sharing_list;
            sharing_list.push_back(current_vtx);

            edge_id_type lookup = traits::vtoe(current_vtx);
            std::size_t i=0;
            for(; i<100; ++i)
            {
                const std::size_t     f = traits::get_face_id(lookup);
                const edge_index_type e = traits::get_edge_index(lookup);
                const edge_id_type  eid = traits::make_edge_id(f, (e==0?2:e-1));
                lookup = this->edge_pairs_[eid]; // update lookup face

                if(traits::get_face_id(lookup) == traits::get_face_id(current_vtx))
                    break;
                const vertex_id_type vid = traits::etov(lookup);
                whole_angle += faces_.at(traits::get_face_id(vid)).angle_at(
                        traits::get_vertex_index(vid));
                sharing_list.push_back(vid);
            }
            if(i == 100)
                throw std::runtime_error("too many faces share a certain vertex");
            vertex_groups_[current_vtx] = std::make_pair(sharing_list, whole_angle);
        }
    }

    // XXX: this is assertion.
    typedef vertex_group_type::const_iterator iter_type;
    for(iter_type iter = vertex_groups_.begin();
            iter != vertex_groups_.end(); ++iter)
    {
        const Real wangle = iter->second.second;
        for(vertex_id_list::const_iterator i = iter->second.first.begin();
                i != iter->second.first.end(); ++i)
        {
            assert(std::abs(vertex_groups_[*i].second - wangle) < 1e-12);
        }
    }

    return;
}

void BDPolygon::detect_shared_edges(const Real tolerance)
{
    std::set<edge_id_type> is_detected;

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
                (start_dist >= tolerance) ? (start_dist) : (tolerance);

            const Real end_dist = length(end_pos);
            const Real end_length_scale =
                (end_dist >= tolerance) ? (end_dist) : (tolerance);

            bool found = false;
            for(std::size_t f = currentf + 1; f < faces_.size(); ++f)
            {
                for(std::size_t e = 0; e < 3; ++e)
                {
                    const Real start_pos_dist =
                        length(start_pos - faces_.at(f).vertex_at(e));
                    if(start_pos_dist > start_length_scale * tolerance)
                        continue;

                    const Real end_pos_dist =
                        length(end_pos - faces_.at(f).vertex_at((e==0)?2:e-1));
                    if(end_pos_dist <= end_length_scale * tolerance)
                    {
                        const edge_id_type detected =
                            std::make_pair(f, (e==0)?2:e-1);
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
