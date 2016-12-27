#ifndef BD_POLYGON
#define BD_POLYGON
#include <ecell4/core/Triangle.hpp>
#include <ecell4/core/get_mapper_mf.hpp>
#include <utility>
#include <functional>
#include <vector>

namespace ecell4
{
namespace bd
{

struct BDPolygon
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
    BDPolygon(){}
    ~BDPolygon(){}

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

    std::pair<bool, uint32_t>
    is_connected(const face_id_type& lhs, const face_id_type& rhs) const;
    std::pair<bool, uint32_t>
    is_share_vertex(const face_id_type& lhs, const face_id_type& rhs) const;

    Real distance(const std::pair<Real3, face_id_type>& lhs,
                  const std::pair<Real3, face_id_type>& rhs) const;


    std::pair<face_id_type, Real3> bend_displacement(
            const std::pair<face_id_type, Real3>& pos, const Real3& disp) const;


    bool empty() const {return faces_.empty();}
    void clear() {return faces_.clear();}
    face_type&       operator[](const std::size_t i)       {return faces_[i];}
    face_type const& operator[](const std::size_t i) const {return faces_[i];}
    face_type&       at(const std::size_t i)       {return faces_.at(i);}
    face_type const& at(const std::size_t i) const {return faces_.at(i);}

    struct face_finder
        : public std::unary_function<vertex_id_type, bool>
    {
        face_finder(const face_id_type& fid): fid_(fid){}
        bool operator()(const vertex_id_type& vid) const
        {
            return vid.first == fid_;
        }
      protected:
        face_id_type fid_;
    };

  private:

    face_container_type faces_;
    edge_pair_type      edge_pairs_;
    vertex_group_type   vertex_groups_;
};

}// bd

}// ecell
#endif //EGFRD_POLYGON
