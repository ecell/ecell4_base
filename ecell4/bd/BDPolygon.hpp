#ifndef EGFRD_POLYGON
#define EGFRD_POLYGON
#include <ecell4/core/Triangle.hpp>
#include <ecell4/core/get_mapper_mf.hpp>
#include <utility>
#include <vector>

namespace ecell4
{
namespace bd
{

struct Polygon
{
  public:
    typedef Triangle face_type;
    typedef std::vector<face_type> face_container_type;
    typedef std::size_t face_id_type;
    typedef std::pair<face_id_type, uint32_t> vertex_id_type;
    typedef std::pair<face_id_type, uint32_t> edge_id_type;
    typedef std::vector<vertex_id_type> vertex_id_list;
    typedef typename utils::get_mapper_mf<edge_id_type, edge_id_type>::type
        edge_pair_type;
    typedef typename utils::get_mapper_mf<vertex_id_type, vertex_id_list>::type
        vertex_group_type;

  public:
    Polygon(){}
    ~Polygon(){}

    void add_face(const face_type& face)
    {
        faces_.push_back(face);
    }

    void detect_connectivity();

    std::size_t num_faces() const
    {
        return faces_.size();
    }

    edge_id_type const& connecting_edge(const edge_id_type& eid)
    {
        return edge_pairs_[eid];
    }

    vertex_id_list const& connecting_vertices(const vertex_id_type& vid)
    {
        return vertex_groups_[vid];
    }

    bool empty() const {return faces_.empty();}
    void clear() {return faces_.clear();}
    face_type&       operator[](const std::size_t i)       {return faces_[i];}
    face_type const& operator[](const std::size_t i) const {return faces_[i];}
    face_type&       at(const std::size_t i)       {return faces_.at(i);}
    face_type const& at(const std::size_t i) const {return faces_.at(i);}

  private:

    face_container_type faces_;
    edge_pair_type      edge_pairs_;
    vertex_group_type   vertex_groups_;
};

}// bd

}// ecell
#endif //EGFRD_POLYGON
