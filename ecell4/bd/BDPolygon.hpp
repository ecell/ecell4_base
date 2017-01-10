#ifndef BD_POLYGON
#define BD_POLYGON
#include <ecell4/core/Triangle.hpp>
#include <ecell4/core/get_mapper_mf.hpp>
#include <ecell4/core/hash.hpp>
#include <utility>
#include <functional>
#include <vector>

ECELL4_DEFINE_HASH_BEGIN()

template<>
struct hash<std::pair<std::size_t, uint32_t> >
{
    typedef std::pair<std::size_t, uint32_t> argument_type;

    std::size_t operator()(argument_type const& val) const
    {
        return hash<std::size_t>()(val.first) ^
               hash<uint32_t>()(val.second);
    }
};

ECELL4_DEFINE_HASH_END()

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


    std::pair<std::pair<Real3, face_id_type>, Real3>
    move_next_face(
        const std::pair<Real3, face_id_type>& pos, const Real3& disp) const;

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


    struct edge_finder
        : public std::unary_function<edge_id_type, bool>
    {
        edge_finder(const edge_id_type& eid): eid_(eid){}
        bool operator()(const std::pair<edge_id_type, edge_id_type>& eid) const
        {
            return eid.first == eid_;
        }
      protected:
        edge_id_type eid_;
    };

  private:

    void detect_shared_vertices();
    void detect_shared_edges();

    boost::array<Real, 3> to_barycentric(const Real3& pos, const face_type& face) const;
    Real3 to_absolute(const boost::array<Real, 3>& pos, const face_type& face) const;
    bool is_inside(const boost::array<Real, 3>& pos) const;
    std::pair<uint32_t, Real>
        crossed_edge(const boost::array<Real, 3>& pos,
                     const boost::array<Real, 3>& disp) const;
    Real cross_section(const boost::array<Real, 3>& pos,
                       const boost::array<Real, 3>& disp, const uint32_t e) const;
    Real triangle_area_2D(const Real& x1, const Real& y1,
        const Real& x2, const Real& y2, const Real& x3, const Real& y3) const;

  private:

    face_container_type faces_;
    edge_pair_type      edge_pairs_;
    vertex_group_type   vertex_groups_;
};


inline Real3 BDPolygon::to_absolute(
        const boost::array<Real, 3>& bary, const face_type& f) const
{
    return f.vertex_at(0) * bary[0] + f.vertex_at(1) * bary[1] +
           f.vertex_at(2) * bary[2];
}

inline bool BDPolygon::is_inside(const boost::array<Real, 3>& barycentric) const
{
    return (0. <= barycentric[0] && barycentric[0] <= 1.0) &&
           (0. <= barycentric[1] && barycentric[1] <= 1.0) &&
           (0. <= barycentric[2] && barycentric[2] <= 1.0);
}

inline Real BDPolygon::triangle_area_2D(const Real& x1, const Real& y1,
    const Real& x2, const Real& y2, const Real& x3, const Real& y3) const
{
    return (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2);
}

}// bd

}// ecell
#endif //EGFRD_POLYGON
