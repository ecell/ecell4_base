#ifndef ECELL4_NGFRD_SHELL_HPP
#define ECELL4_NGFRD_SHELL_HPP
#include <ecell4/core/Sphere.hpp>
#include <ecell4/core/Cylinder.hpp>
#include <ecell4/core/Circle.hpp>
#include <ecell4/core/Cone.hpp>
#include <boost/variant.hpp>

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

    explicit CircularShell(shape_type crc): shape_(crc) {}

    Real3& position()       noexcept {return shape_.position();}
    Real3  position() const noexcept {return shape_.position();}

    shape_type&       shape()       noexcept {return shape_;}
    shape_type const& shape() const noexcept {return shape_;}

  private:

    shape_type shape_;
};

// 2D shell on a vertex
struct ConicalShell
{
    using shape_type = ConicalSurface;

    explicit ConicalShell(shape_type con): shape_(con) {}

    Real3& position()       noexcept {return shape_.position();}
    Real3  position() const noexcept {return shape_.position();}

    shape_type&       shape()       noexcept {return shape_;}
    shape_type const& shape() const noexcept {return shape_;}

  private:

    shape_type shape_;
};

enum class ShellKind : int // boost::variant::which returns an int.
{
    Shperical     = 0,
    Cylindrical   = 1,
    CircularShell = 2,
    ConicalShell  = 3,
    Unknown       = 4
};

struct Shell
{
public:

    using storage_type = boost::variant<
        SphericalShell, CylindricalShell, CircularShell, ConicalShell>;

public:

    template<typename S>
    explicit Shell(S&& s): storage_(std::forward<S>(s)) {}

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

private:

    storage_type storage_;
};

template<typename Visitor, typename ... Shells>
typename Visitor::result_type visit(Visitor&& visitor, Shells&& ... sh)
{
    return boost::apply_visitor(std::forward<Visitor>(visitor),
            std::forward<Shells>(sh).as_variant() ...);
}

} // ngfrd
} // ecell4
#endif//ECELL4_NGFRD_SPHERICAL_SHELL_HPP
