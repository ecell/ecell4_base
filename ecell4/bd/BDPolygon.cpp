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
    {
    std::set<edge_id_type> is_detected;
    for(std::size_t fidx = 0; fidx < faces_.size(); ++fidx)
    for(std::size_t eidx = 0; eidx < 3; ++eidx)
    {
        const edge_id_type current = std::make_pair(fidx, eidx);
        if(is_detected.count(current) == 1) continue;
        const Real3 pos1 = faces_.at(fidx).edge_at(eidx); // start pos
        const Real3 pos2 = faces_.at(fidx).edge_at(eidx==2?0:eidx+1); // end pos

        bool found = false;
        for(std::size_t f=fidx+1; f<faces_.size(); ++f)
        {
            for(std::size_t e=0; e<3; ++e)
            {
                const edge_id_type lookup = std::make_pair(f,e);
                if(is_detected.count(lookup) == 1)
                    continue;
                else if(// each vertex positions are same
                   length(faces_.at(f).vertex_at(e)          - pos2) < 1e-10 &&
                   length(faces_.at(f).vertex_at(e==2?0:e+1) - pos1) < 1e-10)
                {
                    found = true;
                    is_detected.insert(lookup);
                    this->edge_pairs_[current] = lookup;
                    this->edge_pairs_[lookup] = current;
                    break;
                }
            }
            if(found) break;
        }
        if(!found) throw std::logic_error("the polygon is not closed");
    }
    }// edge

    {
    std::set<vertex_id_type> is_detected;
    for(std::size_t fidx = 0; fidx < faces_.size(); ++fidx)
    for(std::size_t vidx = 0; vidx < 3; ++vidx)
    {
        const vertex_id_type current = std::make_pair(fidx, vidx);
        if(is_detected.count(current) == 1) continue;
        this->vertex_groups_[current] = vertex_id_list();
        vertex_id_type lookup = current;
        while(true)
        {
            const edge_id_type prev_edge =
                std::make_pair(lookup.first, lookup.second==0?2:lookup.second-1);
            const typename edge_pair_type::const_iterator pairing =
                std::find(edge_pairs_.begin(), edge_pairs_.end(), prev_edge);
            if(pairing == edge_pairs_.end())
                throw std::logic_error("the polygon is not closed");

            lookup = edge_pairs_[prev_edge]; // XXX dangerous cast!!!
            if(lookup.first == current.first) break;
            this->vertex_groups_[current].push_back(lookup);
            is_detected.insert(lookup);
        }
    }
    }// vertex

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
        if(std::find(list.begin(), list.end(), cmp) != list.end())
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
        return std::make_pair(
                std::make_pair(this->to_absolute(newpos, f), pos.second),
                Real3(0.,0.,0.));

    const boost::array<Real, 3> bpos = to_barycentric(pos.first, f);
    boost::array<Real, 3> bdis;
    bdis[0] = newpos[0] - bpos[0];
    bdis[1] = newpos[1] - bpos[1];
    bdis[2] = newpos[2] - bpos[2];

    const std::pair<uint32_t, Real> cross = crossed_edge(bpos, bdis);
    const face_id_type next_face =
        std::find(edge_pairs_.begin(), edge_pairs_.end(), std::make_pair(f, cross.first)
                )->first.first;

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
 

}// bd
}// ecell
