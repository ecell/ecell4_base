#include "BDPolygon.hpp"
#include <set>

namespace ecell4
{

namespace bd
{

void Polygon::detect_connectivity()
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
        }
    }
    }// vertex

    return;
}



}// bd

}// ecell
