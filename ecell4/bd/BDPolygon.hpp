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

  private:

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
