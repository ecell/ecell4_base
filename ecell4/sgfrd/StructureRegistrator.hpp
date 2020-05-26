#ifndef ECELL4_SGFRD_STRUCTURE_REGISTRATOR
#define ECELL4_SGFRD_STRUCTURE_REGISTRATOR
#include <ecell4/core/Polygon.hpp>
#include <unordered_map>

namespace ecell4
{
namespace sgfrd
{

template<typename T_element_id, typename T_structure_id>
struct StructureRegistrator
{
public:

    typedef T_element_id     element_id_type;
    typedef T_structure_id   structure_id_type;
    typedef std::vector<element_id_type> element_id_array_type;
    typedef std::unordered_map<structure_id_type, element_id_array_type> container_type;
    typedef std::unordered_map<element_id_type, structure_id_type> elemid_to_strid_map_type;

    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;
    typedef ecell4::Polygon polygon_type;
public:

    StructureRegistrator()  = default;
    ~StructureRegistrator() = default;

    void emplace(const element_id_type&, const structure_id_type&); // add new relation
    void update (const element_id_type&, const structure_id_type&); // remove old relation and add new one
    void remove (const element_id_type&, const structure_id_type&); // use hint
    void remove (const element_id_type& eid)
    {
        this->remove(eid, elemid_to_strid_map_.at(eid));
        return ;
    }

    bool have(const element_id_type& eid) const
    {
        return elemid_to_strid_map_.count(eid) != 0;
    }

    element_id_array_type&       elements_over(const structure_id_type& sid)
    {
        return container_.at(sid);
    }
    element_id_array_type const& elements_over(const structure_id_type& sid) const
    {
        return container_.at(sid);
    }
    structure_id_type&           structure_on(const element_id_type& eid)
    {
        return elemid_to_strid_map_.at(eid);
    }
    structure_id_type const&     structure_on(const element_id_type& eid) const
    {
        return elemid_to_strid_map_.at(eid);
    }

    void reset(){elemid_to_strid_map_.clear(); container_.clear();}
    bool empty() const throw() {return container_.empty();}
    std::size_t size() const throw() {return container_.size();}
    void resize(std::size_t i){return container_.resize(i);}

    iterator begin() throw() {return container_.begin();}
    iterator end()   throw() {return container_.begin();}
    const_iterator begin()  const throw() {return container_.begin();}
    const_iterator end()    const throw() {return container_.end();}
    const_iterator cbegin() const throw() {return container_.begin();}
    const_iterator cend()   const throw() {return container_.end();}

    void dump(std::ostream& os) const;

protected:

    elemid_to_strid_map_type elemid_to_strid_map_;   //ex {pID -> fID}
    container_type container_;  // {fid -> {pid,...}, ...}
};


template<typename Te, typename Ts>
void StructureRegistrator<Te, Ts>::emplace(
        const element_id_type& eid, const structure_id_type& sid)
{
    if(this->have(eid))
    {
        throw std::logic_error("already have");
    }
    elemid_to_strid_map_[eid] = sid;

    if(container_.count(sid) == 0)
    {
        container_[sid] = element_id_array_type{};
    }
    container_[sid].push_back(eid);
    return;
}

template<typename Te, typename Ts>
void StructureRegistrator<Te, Ts>::update(
        const element_id_type& eid, const structure_id_type& sid)
{
    // remove older eid-sid relationship
    const structure_id_type old_sid   = elemid_to_strid_map_[eid];
    element_id_array_type&  old_value = container_[old_sid];

    const auto found = std::find(old_value.begin(), old_value.end(), eid);
    assert(found != old_value.end()); // should be found
    old_value.erase(found);

    // add new relationship
    elemid_to_strid_map_[eid] = sid;
    if(container_.count(sid) == 0)
    {
        container_[sid] = element_id_array_type{};
    }
    container_[sid].push_back(eid);
    return;
}

template<typename Te, typename Ts>
void StructureRegistrator<Te, Ts>::remove(
        const element_id_type& eid, const structure_id_type& sid)
{
    element_id_array_type& old_value = container_[sid];

    const auto found = std::find(old_value.begin(), old_value.end(), eid);
    assert(found != old_value.end());
    old_value.erase(found);

    elemid_to_strid_map_.erase(eid);
    return;
}

template<typename Te, typename Ts>
void StructureRegistrator<Te, Ts>::dump(std::ostream& os) const
{
//  elemid_to_strid_map_type elemid_to_strid_map_;   //ex {pID -> fID}
//  container_type           container_;  //ex {<fid, {pid,...}>, ...}
    os << "StructureRegistrator::dump\n";
    os << "{element ID -> structure ID}\n";
    for(const auto& eid_sid : this->elemid_to_strid_map_)
    {
        os << "{ " << eid_sid.first << " -> " << eid_sid.second << " }\n";
    }
    os << std::endl;

    os << "{structure ID -> {list of elements...}}\n";
    for(const auto& sid_es : this->container_)
    {
        os << "{ " << sid_es.first << " -> { ";
        for(const auto& eid : sid_es.second) {os << eid << ' ';}
        os << "}}\n";
    }
    os << std::endl;
    return ;
}

} // sgfrd
} // ecell4
#endif // ECELL4_SGFRD_STRUCTURE_REGISTRATOR
