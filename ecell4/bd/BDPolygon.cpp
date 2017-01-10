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

std::pair<bool, uint32_t>
BDPolygon::is_connected(const face_id_type& lhs, const face_id_type& rhs) const
{
    if(edge_pairs_.find(std::make_pair(lhs, 0))->second.first == rhs)
        return std::make_pair(true, 0);
    if(edge_pairs_.find(std::make_pair(lhs, 1))->second.first == rhs)
        return std::make_pair(true, 1);
    if(edge_pairs_.find(std::make_pair(lhs, 2))->second.first == rhs)
        return std::make_pair(true, 2);
    return std::make_pair(false, 3);
}

std::pair<bool, uint32_t>
BDPolygon::is_share_vertex(const face_id_type& lhs, const face_id_type& rhs) const
{
    for(uint32_t i=0; i<3; ++i)
    {
        vertex_id_list list = vertex_groups_.find(std::make_pair(lhs, 0))->second;
        const face_finder cmp(rhs);
        if(std::find_if(list.begin(), list.end(), cmp) != list.end())
            return std::make_pair(true, i);
    }
    return std::make_pair(false, 3);
}

Real BDPolygon::distance(const std::pair<Real3, face_id_type>& lhs,
                         const std::pair<Real3, face_id_type>& rhs) const
{
    if(lhs.second == rhs.second)
        return length(lhs.first - rhs.first);

    const std::pair<bool, uint32_t> edg =
        this->is_connected(lhs.second, rhs.second);
    if(edg.first)
    {
        const Triangle& lhs_t = faces_[lhs.second];
        const Triangle& rhs_t = faces_[rhs.second];
        const Real ang = angle(lhs_t.normal(), rhs_t.normal());

//         const edge_id_type redg = edge_pairs_.find(edg)->second;
        const Real3 developped = lhs_t.vertex_at(edg.second) + 
            rotate(-ang,
                   rhs.first - lhs_t.vertex_at(edg.second),
                   lhs_t.edge_at(edg.second));
        return length(lhs.first - developped);
    }

    const std::pair<bool, uint32_t> vtx =
        this->is_connected(lhs.second, rhs.second);
    if(vtx.first)
    {
        Real whole_angle = 0.;
        Real angle_sum = 0.;

        bool inter = false;
        Real r0 = 0.;
        Real r1 = 0.;

        vertex_id_list const& list =
            vertex_groups_.find(std::make_pair(lhs.second, vtx.second))->second;
        for(vertex_id_list::const_iterator iter = list.begin();
                iter != list.end(); ++iter) // ccw
        {
            const Triangle& f = this->faces_.at(iter->first);
            whole_angle += f.angle_at(iter->second);
            if(!inter && iter->first == lhs.second)
            {
                inter = true;
                angle_sum += angle(f.vertex_at(iter->second) - lhs.first,
                                   f.edge_at(iter->second==0?2:iter->second-1));
            }
            else if(inter && iter->first == rhs.second)
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
        return r0 * r0 + r1 * r1 - 2. * r0 * r1 * cos(min_angle);
    }

//     throw std::logic_error("cannot determine the distance");
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

    const boost::array<Real, 3> bpos = to_barycentric(pos.first, f);
    boost::array<Real, 3> bdis;
    bdis[0] = newpos[0] - bpos[0];
    bdis[1] = newpos[1] - bpos[1];
    bdis[2] = newpos[2] - bpos[2];

    const std::pair<uint32_t, Real> cross = crossed_edge(bpos, bdis);
    const edge_finder cmp(std::make_pair(pos.second, cross.first));
    const face_id_type next_face =
        std::find_if(edge_pairs_.begin(), edge_pairs_.end(), cmp)->first.first;

    return std::make_pair(std::make_pair(pos.first + disp * cross.second, next_face),
            disp * (1. - cross.second));
}

std::pair<uint32_t, Real>
BDPolygon::crossed_edge(const boost::array<Real, 3>& pos, const boost::array<Real, 3>& disp) const
{
    boost::array<Real, 3> npos;
    npos[0] = pos[0] + disp[0];
    npos[1] = pos[1] + disp[1];
    npos[2] = pos[2] + disp[2];

    if(npos[0] * npos[1] * npos[2] < 0)
    {
             if(npos[0] < 0) return std::make_pair(1, cross_section(pos, disp, 1));
        else if(npos[1] < 0) return std::make_pair(2, cross_section(pos, disp, 2));
        else if(npos[2] < 0) return std::make_pair(0, cross_section(pos, disp, 0));
    }
    else
    {
        assert(!is_inside(npos));
        if(npos[0] > 0)
        {
            const Real ab = cross_section(pos, disp, 0);
            const Real ca = cross_section(pos, disp, 2);
            return (ca > ab) ? std::make_pair(0, ab) : std::make_pair(2, ca);
        }
        else if(npos[1] > 0)
        {
            const Real ab = cross_section(pos, disp, 0);
            const Real bc = cross_section(pos, disp, 1);
            return (ab > bc) ? std::make_pair(1, bc) : std::make_pair(0, ab);
        }
        else if(npos[2] > 0)
        {
            const Real bc = cross_section(pos, disp, 1);
            const Real ca = cross_section(pos, disp, 2);
            return (bc > ca) ? std::make_pair(2, ca) : std::make_pair(1, bc);
        }
    }
}

Real BDPolygon::cross_section(const boost::array<Real, 3>& pos,
        const boost::array<Real, 3>& disp, const uint32_t e) const
{
    switch(e)
    {
        case 0:  return -pos[1] / disp[1];
        case 1:  return -pos[0] / disp[0];
        case 2:  return -pos[2] / disp[2];
        default: return std::numeric_limits<Real>::infinity();
    }
}

boost::array<Real, 3>
BDPolygon::to_barycentric(const Real3& pos, const face_type& face) const
{
    const Real3& a = face.vertex_at(0);
    const Real3& b = face.vertex_at(1);
    const Real3& c = face.vertex_at(2);
    const Real3  m = cross_product(face.edge_at(0), face.edge_at(2)) * (-1.);
    const Real x = std::abs(m[0]);
    const Real y = std::abs(m[1]);
    const Real z = std::abs(m[2]);

    Real nu, nv, ood;
    if (x >= y && x >= z)
    {
        nu = triangle_area_2D(pos[1], pos[2], b[1], b[2], c[1], c[2]);
        nv = triangle_area_2D(pos[1], pos[2], c[1], c[2], a[1], a[2]);
        ood = 1.0 / m[0];
    }
    else if (y >= x && y >= z)
    {
        nu = triangle_area_2D(pos[0], pos[2], b[0], b[2], c[0], c[2]);
        nv = triangle_area_2D(pos[0], pos[2], c[0], c[2], a[0], a[2]);
        ood = 1.0 / -m[1];
    }
    else
    {
        nu = triangle_area_2D(pos[0], pos[1], b[0], b[1], c[0], c[1]);
        nv = triangle_area_2D(pos[0], pos[1], c[0], c[1], a[0], a[1]);
        ood = 1.0 / m[2];
    }
    boost::array<Real, 3> bary;
    bary[0] = nu * ood;
    bary[1] = nv * ood;
    bary[2] = 1.0 - bary[0] - bary[1];
    return bary;
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
