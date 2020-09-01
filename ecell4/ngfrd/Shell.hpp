#ifndef ECELL4_NGFRD_SHELL_HPP
#define ECELL4_NGFRD_SHELL_HPP
#include <ecell4/core/Sphere.hpp>
#include <ecell4/core/Cylinder.hpp>
#include <ecell4/core/Circle.hpp>
#include <ecell4/core/Cone.hpp>
#include <ecell4/core/Polygon.hpp>
#include <ecell4/ngfrd/DomainID.hpp>
#include <boost/variant.hpp>
#include <boost/optional.hpp>

namespace ecell4
{
namespace ngfrd
{

struct SphericalShell
{
    using shape_type = Sphere;

    explicit SphericalShell(shape_type sph): shape_(sph) {}

    Real3& position()       noexcept {return shape_.position();}
    Real3  position() const noexcept {return shape_.position();}

    shape_type&       shape()       noexcept {return shape_;}
    shape_type const& shape() const noexcept {return shape_;}

  private:
    shape_type shape_;
};

struct CylindricalShell
{
    using shape_type = Cylinder;

    explicit CylindricalShell(shape_type cyl): shape_(cyl) {}

    Real3& position()       noexcept {return shape_.position();}
    Real3  position() const noexcept {return shape_.position();}

    shape_type&       shape()       noexcept {return shape_;}
    shape_type const& shape() const noexcept {return shape_;}

  private:

    shape_type shape_;
};

// 2D Shell on a face
struct CircularShell
{
    using shape_type = Circle;

    CircularShell(shape_type crc, const FaceID& fid): fid_(fid), shape_(crc) {}

    Real3& position()       noexcept {return shape_.position();}
    Real3  position() const noexcept {return shape_.position();}

    shape_type&       shape()       noexcept {return shape_;}
    shape_type const& shape() const noexcept {return shape_;}

    FaceID&       fid()       noexcept {return fid_;}
    FaceID const& fid() const noexcept {return fid_;}

  private:

    FaceID     fid_;
    shape_type shape_;
};

// 2D shell on a vertex
struct ConicalShell
{
    using shape_type = ConicalSurface;

    ConicalShell(shape_type con, const VertexID& vid): vid_(vid), shape_(con) {}

    Real3& position()       noexcept {return shape_.position();}
    Real3  position() const noexcept {return shape_.position();}

    shape_type&       shape()       noexcept {return shape_;}
    shape_type const& shape() const noexcept {return shape_;}

    VertexID&       vid()       noexcept {return vid_;}
    VertexID const& vid() const noexcept {return vid_;}

  private:

    VertexID   vid_;
    shape_type shape_;
};

enum class ShellKind : int // boost::variant::which returns an int.
{
    Shperical   = 0,
    Cylindrical = 1,
    Circular    = 2,
    Conical     = 3,
    Unknown     = 4
};

struct Shell
{
public:

    using storage_type = boost::variant<
        SphericalShell, CylindricalShell, CircularShell, ConicalShell>;

public:

    template<typename S>
    explicit Shell(S&& s): did_(boost::none), storage_(std::forward<S>(s)) {}

    template<typename S>
    Shell(S&& s, const DomainID& did): did_(did), storage_(std::forward<S>(s)) {}

    ShellKind kind() const noexcept {return ShellKind(storage_.which());}

    bool is_spherical  () const noexcept {return storage_.which() == 0;}
    bool is_cylindrical() const noexcept {return storage_.which() == 1;}
    bool is_circular   () const noexcept {return storage_.which() == 2;}
    bool is_conical    () const noexcept {return storage_.which() == 3;}

    SphericalShell const&   as_spherical()   const {return boost::get<SphericalShell  >(storage_);}
    SphericalShell&         as_spherical()         {return boost::get<SphericalShell  >(storage_);}
    CylindricalShell const& as_cylindrical() const {return boost::get<CylindricalShell>(storage_);}
    CylindricalShell&       as_cylindrical()       {return boost::get<CylindricalShell>(storage_);}
    CircularShell const&    as_circular()    const {return boost::get<CircularShell   >(storage_);}
    CircularShell&          as_circular()          {return boost::get<CircularShell   >(storage_);}
    ConicalShell const&     as_conical()     const {return boost::get<ConicalShell    >(storage_);}
    ConicalShell&           as_conical()           {return boost::get<ConicalShell    >(storage_);}

    storage_type const& as_variant() const noexcept {return storage_;}
    storage_type&       as_variant()       noexcept {return storage_;}

    Real3& position()
    {
        switch(storage_.which())
        {
            case 0: {return boost::get<SphericalShell  >(storage_).position();}
            case 1: {return boost::get<CylindricalShell>(storage_).position();}
            case 2: {return boost::get<CircularShell   >(storage_).position();}
            case 3: {return boost::get<ConicalShell    >(storage_).position();}
            default:{throw std::runtime_error("Shell::position: bad_visit");}
        }
    }
    Real3 position() const
    {
        switch(storage_.which())
        {
            case 0: {return boost::get<SphericalShell  >(storage_).position();}
            case 1: {return boost::get<CylindricalShell>(storage_).position();}
            case 2: {return boost::get<CircularShell   >(storage_).position();}
            case 3: {return boost::get<ConicalShell    >(storage_).position();}
            default:{throw std::runtime_error("Shell::position: bad_visit");}
        }
    }

    boost::optional<DomainID> const& domain_id() const noexcept {return did_;}
    boost::optional<DomainID>&       domain_id()       noexcept {return did_;}

private:

    boost::optional<DomainID> did_;
    storage_type storage_;
};

template<typename Visitor, typename ... Shells>
typename std::remove_const<typename std::remove_reference<Visitor>::type>::type::result_type
visit(Visitor&& visitor, Shells&& ... sh)
{
    return boost::apply_visitor(std::forward<Visitor>(visitor),
            std::forward<Shells>(sh).as_variant() ...);
}

struct ShellDistanceCalculator
    : boost::static_visitor<Real>
{
    Real3           pos;
    Boundary const* boundary;

    ShellDistanceCalculator(const Real3& p, const Boundary& b)
        : pos(p), boundary(&b)
    {}

    // if the position is inside of the sphere, returns a negative value.
    Real operator()(const SphericalShell& sh) const noexcept
    {
        return length(pos - boundary->periodic_transpose(sh.position(), pos)) -
               sh.shape().radius();
    }
    // if the position is inside of the sphere, returns a negative value.
    Real operator()(const CylindricalShell& sh) const noexcept
    {
        // a candidate of nearest point
        const auto p1 = boundary->periodic_transpose(pos, sh.position());
        const auto d  = distance(sh.shape(), p1);

        const auto& cyl_r = sh.shape().radius();
        const auto& cyl_h = sh.shape().half_height();
        const Real  R = std::sqrt(cyl_r * cyl_r + cyl_h * cyl_h) + d;

        const auto& edge = boundary->edge_lengths();

        // likely
        if(2 * R < edge[0] && 2 * R < edge[1] && 2 * R < edge[2])
        {
            return d;
        }

        // an AABB that includes minmaxdist region. If a point exceeds this,
        // the point will never be the mindist point.
        const Real3 lower(sh.position() - Real3(R,R,R));
        const Real3 upper(sh.position() + Real3(R,R,R));

        // check all the mirror image. using the AABB we constructed, we can
        // skip most of the images.
        Real dist = d;
        for(int ix=-1; ix<=1; ++ix)
        {
            const Real px = p1[0] + ix * edge[0];
            if(px < lower[0] || upper[0] < px) {continue;}
            for(int iy=-1; iy<=1; ++iy)
            {
                const Real py = p1[1] + iy * edge[1];
                if(py < lower[1] || upper[1] < py) {continue;}

                for(int iz=-1; iz<=1; ++iz)
                {
                    const Real pz = p1[2] + iz * edge[2];
                    if(pz < lower[2] || upper[2] < pz) {continue;}
                    dist = std::min(dist,
                            this->distance(sh.shape(), Real3(px,py,pz)));
                }
            }
        }
        return dist;
    }
    Real operator()(const CircularShell& sh) const noexcept
    {
        // sometimes circle wraps two triangles. it is difficult to calculate
        // distance between position and such a shell. But, in any case, shell
        // does not go beyond the bounding sphere. Here, it returns a distance
        // to the bounding sphere. So it UNDER-ESTIMATES the distance.

        return length(pos - boundary->periodic_transpose(sh.position(), pos)) -
               sh.shape().radius();
    }
    Real operator()(const ConicalShell& sh) const noexcept
    {
        // sometimes circle wraps two triangles. it is difficult to calculate
        // distance between position and such a shell. But, in any case, shell
        // does not go beyond the bounding sphere. Here, it returns a distance
        // to the bounding sphere. So it UNDER-ESTIMATES the distance.

        return length(pos - boundary->periodic_transpose(sh.position(), pos)) -
               sh.shape().slant_height();
    }

private:

    // a helper function to calculate cylinder-point distance
    static Real distance(const Cylinder& cyl, const Real3& pos) noexcept
    {
        const auto dr = pos - cyl.center();
        const auto ax = cyl.axis() / length(cyl.axis());
        const auto z  = ax * dot_product(dr, ax);
        const auto r  = dr - z;
        const auto lz = length(z);
        const auto lr = length(r);

        //  | A  |  B
        //   .--.  ___  A) distance to the center - half height
        //  |'--'|  C   B) distance to the edge
        //  |    | ___  C) distance to the axis - radius
        //   '--'

        if(lz < cyl.half_height()) // C
        {
            return lr - cyl.radius();
        }
        else if(lr < cyl.radius()) // A
        {
            return lz - cyl.half_height();
        }
        else // B
        {
            const auto a = lz - cyl.half_height();
            const auto b = lr - cyl.radius();
            return std::sqrt(a * a + b * b);
        }
    }
};

} // ngfrd
} // ecell4
#endif//ECELL4_NGFRD_SPHERICAL_SHELL_HPP
