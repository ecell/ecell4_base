#ifndef ECELL4_PYTHON_API_SHAPE_HPP
#define ECELL4_PYTHON_API_SHAPE_HPP

#include <pybind11/pybind11.h>
#include <ecell4/core/Shape.hpp>
#include <ecell4/core/shape_operators.hpp>
#include <ecell4/core/Sphere.hpp>
#include <ecell4/core/Cylinder.hpp>
#include <ecell4/core/PlanarSurface.hpp>
#include <ecell4/core/Rod.hpp>
#include <ecell4/core/AABB.hpp>
#include <ecell4/core/Mesh.hpp>

namespace py = pybind11;

namespace ecell4
{

namespace python_api
{

    template<class Base = Shape>
    class PyShape: public Base
    {
    public:
        using Base::Base;
        using dimension_kind = Shape::dimension_kind;

        dimension_kind dimension() const override
        {
            PYBIND11_OVERLOAD_PURE(dimension_kind, Base, dimension,);
        }

        Real is_inside(const Real3& coord) const override
        {
            PYBIND11_OVERLOAD_PURE(Real, Base, is_inside, coord);
        }

        Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const override
        {
            PYBIND11_OVERLOAD_PURE(Real3, Base, draw_position, rng);
        }

        bool test_AABB(const Real3& l, const Real3& u) const override
        {
            PYBIND11_OVERLOAD_PURE(bool, Base, test_AABB, l, u);
        }

        void bounding_box(const Real3& edge_lengths, Real3& lower, Real3& upper) const override
        {
            PYBIND11_OVERLOAD(void, Base, bounding_box, edge_lengths, lower, upper);
        }
    };

    void define_shape(py::module& m)
    {
        py::class_<Shape, PyShape<>, boost::shared_ptr<Shape>>(m, "Shape")
            .def("dimension", [](const Shape& self) { return static_cast<Integer>(self.dimension()); })
            .def("is_inside", &Shape::is_inside)
            ;

        py::class_<Surface, PyShape<Surface>, boost::shared_ptr<Surface>>(m, "Surface")
            .def("root", &Surface::root)
            .def(py::pickle(
                [](const Surface& self)
                {
                    return py::make_tuple(self.root());
                },
                [](py::tuple t)
                {
                    if (t.size() != 1)
                        throw std::runtime_error("Invalid state");
                    return Surface(t[0].cast<const boost::shared_ptr<Shape>&>());
                }
            ));

        py::class_<Union, PyShape<Union>, boost::shared_ptr<Union>>(m, "Union")
            .def(py::init<const boost::shared_ptr<Shape>&, const boost::shared_ptr<Shape>&>())
            .def("surface", &Union::surface)
            .def("one", &Union::one)
            .def("another", &Union::another)
            .def(py::pickle(
                [](const Union& self)
                {
                    return py::make_tuple(self.one(), self.another());
                },
                [](py::tuple t)
                {
                    if (t.size() != 2)
                        throw std::runtime_error("Invalid state");
                    return Union(
                        t[0].cast<const boost::shared_ptr<Shape>&>(),
                        t[1].cast<const boost::shared_ptr<Shape>&>()
                    );
                }
            ));

        py::class_<Complement, PyShape<Complement>, boost::shared_ptr<Complement>>(m, "Complement")
            .def(py::init<const boost::shared_ptr<Shape>&, const boost::shared_ptr<Shape>&>())
            .def("surface", &Complement::surface)
            .def("one", &Complement::one)
            .def("another", &Complement::another)
            .def(py::pickle(
                [](const Complement& self)
                {
                    return py::make_tuple(self.one(), self.another());
                },
                [](py::tuple t)
                {
                    if (t.size() != 2)
                        throw std::runtime_error("Invalid state");
                    return Complement(
                        t[0].cast<const boost::shared_ptr<Shape>&>(),
                        t[1].cast<const boost::shared_ptr<Shape>&>()
                    );
                }
            ));

        py::class_<AffineTransformation, PyShape<AffineTransformation>,
            boost::shared_ptr<AffineTransformation>>(m, "AffineTransformation")
            .def(py::init<>())
            .def(py::init<const boost::shared_ptr<Shape>&>())
            .def(py::init<const boost::shared_ptr<Shape>&, const Real3&, const Real3&, const Real3&, const Real3&>())
            .def("translate", &AffineTransformation::translate)
            .def("rescale", &AffineTransformation::rescale)
            .def("xroll", &AffineTransformation::xroll)
            .def("yroll", &AffineTransformation::yroll)
            .def("zroll", &AffineTransformation::zroll)
            .def("surface", &AffineTransformation::surface)
            .def("root", &AffineTransformation::root)
            .def("first", &AffineTransformation::first)
            .def("second", &AffineTransformation::second)
            .def("third", &AffineTransformation::third)
            .def("shift", &AffineTransformation::shift)
            .def(py::pickle(
                [](const AffineTransformation& self)
                {
                    return py::make_tuple(self.root(), self.first(), self.second(), self.third(), self.shift());
                },
                [](py::tuple t)
                {
                    if (t.size() != 5)
                        throw std::runtime_error("Invalid state");
                    return AffineTransformation(
                        t[0].cast<const boost::shared_ptr<Shape>&>(),
                        t[1].cast<const Real3&>(),
                        t[2].cast<const Real3&>(),
                        t[3].cast<const Real3&>(),
                        t[4].cast<const Real3&>()
                    );
                }
            ));

        py::class_<Sphere, PyShape<Sphere>, boost::shared_ptr<Sphere>>(m, "Sphere")
            .def(py::init<const Real3&, const Real>())
            .def("distance", &Sphere::distance)
            .def("surface", &Sphere::surface)
            .def("center", &Sphere::center)
            .def("radius", &Sphere::radius)
            .def(py::pickle(
                [](const Sphere& self)
                {
                    return py::make_tuple(self.center(), self.radius());
                },
                [](py::tuple t)
                {
                    if (t.size() != 2)
                        throw std::runtime_error("Invalid state");
                    return Sphere(t[0].cast<const Real3&>(), t[1].cast<const Real>());
                }
            ));

        py::class_<SphericalSurface, PyShape<SphericalSurface>,
            boost::shared_ptr<SphericalSurface>>(m, "SphericalSurface")
            .def(py::init<const Real3&, const Real>())
            .def("distance", &SphericalSurface::distance)
            .def("inside", &SphericalSurface::inside)
            .def("center", &SphericalSurface::center)
            .def("radius", &SphericalSurface::radius)
            .def(py::pickle(
                [](const SphericalSurface& self)
                {
                    return py::make_tuple(self.center(), self.radius());
                },
                [](py::tuple t)
                {
                    if (t.size() != 2)
                        throw std::runtime_error("Invalid state");
                    return SphericalSurface(t[0].cast<const Real3&>(), t[1].cast<const Real>());
                }
            ));

        py::class_<Cylinder, PyShape<Cylinder>, boost::shared_ptr<Cylinder>>(m, "Cylinder")
            .def(py::init<const Real3&, const Real, const Real3&, const Real>())
            .def("distance", &Cylinder::distance)
            .def("surface", &Cylinder::surface)
            .def("center", &Cylinder::center)
            .def("axis", &Cylinder::axis)
            .def("half_height", &Cylinder::half_height)
            .def(py::pickle(
                [](const Cylinder& self)
                {
                    return py::make_tuple(self.center(), self.radius(), self.axis(), self.half_height());
                },
                [](py::tuple t)
                {
                    if (t.size() != 4)
                        throw std::runtime_error("Invalid state");
                    return Cylinder(
                        t[0].cast<const Real3&>(),
                        t[1].cast<const Real>(),
                        t[2].cast<const Real3&>(),
                        t[3].cast<const Real>()
                    );
                }
            ));

        py::class_<CylindricalSurface, PyShape<CylindricalSurface>,
            boost::shared_ptr<CylindricalSurface>>(m, "CylindricalSurface")
            .def(py::init<const Real3&, const Real, const Real3&, const Real>())
            .def("distance", &CylindricalSurface::distance)
            .def("inside", &CylindricalSurface::inside)
            .def("center", &CylindricalSurface::center)
            .def("radius", &CylindricalSurface::radius)
            .def("axis", &CylindricalSurface::axis)
            .def("half_height", &CylindricalSurface::half_height)
            .def(py::pickle(
                [](const CylindricalSurface& self)
                {
                    return py::make_tuple(self.center(), self.radius(), self.axis(), self.half_height());
                },
                [](py::tuple t)
                {
                    if (t.size() != 4)
                        throw std::runtime_error("Invalid state");
                    return CylindricalSurface(
                        t[0].cast<const Real3&>(),
                        t[1].cast<const Real>(),
                        t[2].cast<const Real3&>(),
                        t[3].cast<const Real>()
                    );
                }
            ));

        py::class_<PlanarSurface, PyShape<PlanarSurface>,
            boost::shared_ptr<PlanarSurface>>(m, "PlanarSurface")
            .def(py::init<const Real3&, const Real3&, const Real3&>())
            .def("origin", &PlanarSurface::origin)
            .def("e0", &PlanarSurface::e0)
            .def("e1", &PlanarSurface::e1)
            .def("normal", &PlanarSurface::normal)
            .def(py::pickle(
                [](const PlanarSurface& self)
                {
                    return py::make_tuple(self.origin(), self.e0(), self.e1());
                },
                [](py::tuple t)
                {
                    if (t.size() != 3)
                        throw std::runtime_error("Invalid state");
                    return PlanarSurface(
                        t[0].cast<const Real3&>(),
                        t[1].cast<const Real3&>(),
                        t[2].cast<const Real3&>()
                    );
                }
            ));

        py::class_<Rod, PyShape<Rod>, boost::shared_ptr<Rod>>(m, "Rod")
            .def(py::init<const Real&, const Real&>())
            .def(py::init<const Real&, const Real&, const Real3&>())
            .def("distance", &Rod::distance)
            .def("origin", &Rod::origin)
            .def("length", &Rod::length)
            .def("radius", &Rod::radius)
            .def("shift", &Rod::shift)
            .def("surface", &Rod::surface)
            .def(py::pickle(
                [](const Rod& self)
                {
                    return py::make_tuple(self.length(), self.radius(), self.origin());
                },
                [](py::tuple t)
                {
                    if (t.size() != 3)
                        throw std::runtime_error("Invalid state");
                    return Rod(
                        t[0].cast<const Real>(),
                        t[1].cast<const Real>(),
                        t[2].cast<const Real3&>()
                    );
                }
            ));

        py::class_<RodSurface, PyShape<RodSurface>, boost::shared_ptr<RodSurface>>(m, "RodSurface")
            .def(py::init<const Real&, const Real&>())
            .def(py::init<const Real&, const Real&, const Real3&>())
            .def("distance", &RodSurface::distance)
            .def("origin", &RodSurface::origin)
            .def("length", &RodSurface::length)
            .def("radius", &RodSurface::radius)
            .def("shift", &RodSurface::shift)
            .def(py::pickle(
                [](const RodSurface& self)
                {
                    return py::make_tuple(self.length(), self.radius(), self.origin());
                },
                [](py::tuple t)
                {
                    if (t.size() != 3)
                        throw std::runtime_error("Invalid state");
                    return RodSurface(
                        t[0].cast<const Real>(),
                        t[1].cast<const Real>(),
                        t[2].cast<const Real3&>()
                    );
                }
            ));

        py::class_<AABB, PyShape<AABB>, boost::shared_ptr<AABB>>(m, "AABB")
            .def(py::init<const Real3&, const Real3&>())
            .def("distance", &AABB::distance)
            .def("upper", &AABB::upper)
            .def("lower", &AABB::lower)
            .def("surface", &AABB::surface)
            .def(py::pickle(
                [](const AABB& self)
                {
                    return py::make_tuple(self.lower(), self.upper());
                },
                [](py::tuple t)
                {
                    if (t.size() != 2)
                        throw std::runtime_error("Invalid state");
                    return AABB(
                        t[0].cast<const Real3&>(),
                        t[1].cast<const Real3&>()
                    );
                }
            ));

        py::class_<MeshSurface, PyShape<MeshSurface>, boost::shared_ptr<MeshSurface>>(m, "MeshSurface")
            .def(py::init<const std::string, const Real3&>())
            .def("filename", &MeshSurface::filename)
            .def("edge_lengths", &MeshSurface::edge_lengths)
            .def(py::pickle(
                [](const MeshSurface& self)
                {
                    return py::make_tuple(self.filename(), self.edge_lengths());
                },
                [](py::tuple t)
                {
                    if (t.size() != 2)
                        throw std::runtime_error("Invalid state");
                    return MeshSurface(
                        t[0].cast<const std::string>(),
                        t[1].cast<const Real3&>()
                    );
                }
            ));

        m.def("create_x_plane",
            [](Real x)
            {
                return PlanarSurface(Real3(x, 0, 0), Real3(0, 1, 0), Real3(0, 0, 1));
            });

        m.def("create_y_plane",
            [](Real y)
            {
                return PlanarSurface(Real3(0, y, 0), Real3(1, 0, 0), Real3(0, 0, 1));
            });

        m.def("create_z_plane",
            [](Real z)
            {
                return PlanarSurface(Real3(0, 0, z), Real3(1, 0, 0), Real3(0, 1, 0));
            });
    }
}

}

#endif /* ECELL4_PYTHON_API_SHAPE_HPP */
