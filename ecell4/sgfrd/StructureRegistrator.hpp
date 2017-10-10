#ifndef ECELL4_SGFRD_STRUCTURE_REGISTRATOR
#define ECELL4_SGFRD_STRUCTURE_REGISTRATOR
#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/get_mapper_mf.hpp>

namespace ecell4
{
namespace sgfrd
{

template<typename T_element_id, typename T_structure_id,
         typename T_polygon_traits>
struct StructureRegistrator
{
public:

    typedef T_element_id     element_id_type;
    typedef T_structure_id   structure_id_type;
    typedef T_polygon_traits polygon_traits_type;
    typedef std::vector<element_id_type> element_id_array_type;
    typedef std::pair<structure_id_type, element_id_array_type> value_type;
    typedef std::vector<value_type>                 container_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

    typedef ecell4::Polygon<polygon_traits_type> polygon_type;
    typedef typename polygon_traits_type::converter_type
        structure_id_converter_type;
    typedef typename utils::get_mapper_mf<element_id_type, structure_id_type>::type
        elemid_to_strid_map_type;

public:

    StructureRegistrator(const polygon_type& poly)
        : structure_id_converter_(poly.converter()),
          container_(poly.num_faces())
    {}
    ~StructureRegistrator(){}

    void emplace(const element_id_type&, const structure_id_type&);
    void update(const element_id_type&, const structure_id_type&);
    void remove(const element_id_type&, const structure_id_type&); // use hint
    void remove(const element_id_type&);

    bool have(const element_id_type&) const;

    element_id_array_type&       elements_over(const structure_id_type&);
    element_id_array_type const& elements_over(const structure_id_type&) const;
    structure_id_type&           structure_on(const element_id_type&);
    structure_id_type const&     structure_on(const element_id_type&) const;

    void reset(){elemid_to_strid_map_.clear(); container_.clear();}
    bool empty() const throw() {return container_.empty();}
    std::size_t size() const throw() {return container_.size();}
    void resize(std::size_t i){return container_.resize(i);}

    value_type&       operator[](std::size_t i)       throw() {return container_[i];}
    value_type const& operator[](std::size_t i) const throw() {return container_[i];}
    value_type&       at(std::size_t i)       {return container_.at(i);}
    value_type const& at(std::size_t i) const {return container_.at(i);}

    element_id_array_type&       element_ids_at(std::size_t i);
    element_id_array_type const& element_ids_at(std::size_t i) const;
    structure_id_type&           structure_id_at(std::size_t i);
    structure_id_type const&     structure_id_at(std::size_t i) const;

    iterator begin() throw() {return container_.begin();}
    iterator end()   throw() {return container_.begin();}
    const_iterator begin()  const throw() {return container_.begin();}
    const_iterator end()    const throw() {return container_.end();}
    const_iterator cbegin() const throw() {return container_.begin();}
    const_iterator cend()   const throw() {return container_.end();}

    void dump(std::ostream& os) const;

protected:

    std::size_t to_index(const structure_id_type& sid) const
    {
        return structure_id_converter_.to_index(sid);
    }

    std::size_t to_index(const element_id_type& eid) const
    {
        typename elemid_to_strid_map_type::const_iterator iter(
            elemid_to_strid_map_.find(eid));
        if(iter == elemid_to_strid_map_.end())
            throw std::out_of_range("no structure id");
        return to_index(iter->second);
    }

protected:

    structure_id_converter_type const& structure_id_converter_;//ex {fID -> idx}
    elemid_to_strid_map_type           elemid_to_strid_map_;   //ex {pID -> fID}
    container_type                     container_;  //ex {<fid, {pid,...}>, ...}
};


template<typename Te, typename Ts, typename Tp>
void StructureRegistrator<Te, Ts, Tp>::emplace(
        const element_id_type& eid, const structure_id_type& sid)
{
    if(this->have(eid))
        throw std::logic_error("already have");
    const std::size_t idx = this->to_index(sid);
    if(container_.size() <= idx) container_.resize(idx+1);
    elemid_to_strid_map_[eid] = sid;

    value_type& contained = container_[idx];
    contained.first = sid;
    contained.second.push_back(eid);

    return;
}

template<typename Te, typename Ts, typename Tp>
void StructureRegistrator<Te, Ts, Tp>::update(
        const element_id_type& eid, const structure_id_type& sid)
{
    const std::size_t idx = this->to_index(sid);
    if(container_.size() <= idx) container_.resize(idx+1);

    // cleanup eid->sid map
    const structure_id_type old_sid = elemid_to_strid_map_[eid];
    value_type& old_value = container_[this->to_index(old_sid)];

    // cleanup container.second
    const typename element_id_array_type::iterator found =
        std::find(old_value.second.begin(), old_value.second.end(), eid);
    assert(found != old_value.second.end());
    old_value.second.erase(found);

    // update
    container_[idx].first = sid;
    container_[idx].second.push_back(eid);
    elemid_to_strid_map_[eid] = sid;
    return;
}

template<typename Te, typename Ts, typename Tp>
void StructureRegistrator<Te, Ts, Tp>::remove(
        const element_id_type& eid, const structure_id_type& sid)
{
    const std::size_t idx = this->to_index(sid);
    value_type& old_value = container_[idx];

    const typename element_id_array_type::iterator found =
        std::find(old_value.second.begin(), old_value.second.end(), eid);
    assert(found != old_value.second.end());
    old_value.second.erase(found);

    elemid_to_strid_map_.erase(eid);
    return;
}

template<typename Te, typename Ts, typename Tp>
void StructureRegistrator<Te, Ts, Tp>::remove(const element_id_type& eid)
{
    remove(eid, structure_id_at(to_index(eid)));
    return;
}

template<typename Te, typename Ts, typename Tp>
inline bool StructureRegistrator<Te, Ts, Tp>::have(const element_id_type& eid) const
{
    return elemid_to_strid_map_.count(eid) == 1;
}

template<typename Te, typename Ts, typename Tp>
inline typename StructureRegistrator<Te, Ts, Tp>::element_id_array_type&
StructureRegistrator<Te, Ts, Tp>::elements_over(const structure_id_type& sid)
{
    return element_ids_at(this->to_index(sid));
}

template<typename Te, typename Ts, typename Tp>
inline typename StructureRegistrator<Te, Ts, Tp>::element_id_array_type const&
StructureRegistrator<Te, Ts, Tp>::elements_over(const structure_id_type& sid) const
{
    return element_ids_at(this->to_index(sid));
}

template<typename Te, typename Ts, typename Tp>
inline typename StructureRegistrator<Te, Ts, Tp>::structure_id_type&
StructureRegistrator<Te, Ts, Tp>::structure_on(const element_id_type& eid)
{
    return structure_id_at(this->to_index(eid));
}

template<typename Te, typename Ts, typename Tp>
inline typename StructureRegistrator<Te, Ts, Tp>::structure_id_type const&
StructureRegistrator<Te, Ts, Tp>::structure_on(const element_id_type& eid) const
{
    return structure_id_at(this->to_index(eid));
}

template<typename Te, typename Ts, typename Tp>
inline typename StructureRegistrator<Te, Ts, Tp>::element_id_array_type&
StructureRegistrator<Te, Ts, Tp>::element_ids_at(std::size_t i)
{
    return container_.at(i).second;
}

template<typename Te, typename Ts, typename Tp>
inline typename StructureRegistrator<Te, Ts, Tp>::element_id_array_type const&
StructureRegistrator<Te, Ts, Tp>::element_ids_at(std::size_t i) const
{
    return container_.at(i).second;
}

template<typename Te, typename Ts, typename Tp>
inline typename StructureRegistrator<Te, Ts, Tp>::structure_id_type&
StructureRegistrator<Te, Ts, Tp>::structure_id_at(std::size_t i)
{
    return container_.at(i).first;
}

template<typename Te, typename Ts, typename Tp>
inline typename StructureRegistrator<Te, Ts, Tp>::structure_id_type const&
StructureRegistrator<Te, Ts, Tp>::structure_id_at(std::size_t i) const
{
    return container_.at(i).first;
}

template<typename Te, typename Ts, typename Tp>
void StructureRegistrator<Te, Ts, Tp>::dump(std::ostream& os) const
{
//     elemid_to_strid_map_type           elemid_to_strid_map_;   //ex {pID -> fID}
//     container_type                     container_;  //ex {<fid, {pid,...}>, ...}
    os << "StructureRegistrator::dump\n";
    os << "{element ID -> structure ID}\n";
    for(typename elemid_to_strid_map_type::const_iterator
        i(elemid_to_strid_map_.begin()), e(elemid_to_strid_map_.end()); i!=e; ++i)
    {
        os << "{ " << i->first << " -> " << i->second " }\n";
    }
    os << std::endl;

    os << "{structure ID -> {list of elements...}}\n";
    for(typename container_type::const_iterator
        i(container_.begin()), e(container_.end()); i!=e; ++i)
    {
        os << "{ " << i->first << " -> { ";
        for(typename element_id_array_type::const_iterator
            ei(i->second.begin()), ee(i->second.end()); ei != ee; ++ei)
        {
            os << *ei << ' ';
        }
        os << " }\n"
    }
    os << std::endl;
    return ;
}

} // sgfrd
} // ecell4
#endif // ECELL4_SGFRD_STRUCTURE_REGISTRATOR
