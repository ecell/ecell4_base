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

        const edge_id_type redg = edge_pairs_.find(edg)->second;
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


}// bd

}// ecell
