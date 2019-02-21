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

    template<class Base>
    class PyShapeImpl: public PyShape<Base>
    {
    public:
        using PyShape<Base>::PyShape;
        using dimension_kind = Shape::dimension_kind;

        dimension_kind dimension() const override
        {
            PYBIND11_OVERLOAD(dimension_kind, Base, dimension,);
        }

        Real is_inside(const Real3& coord) const override
        {
            PYBIND11_OVERLOAD(Real, Base, is_inside, coord);
        }

        Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const override
        {
            PYBIND11_OVERLOAD(Real3, Base, draw_position, rng);
        }

        bool test_AABB(const Real3& l, const Real3& u) const override
        {
            PYBIND11_OVERLOAD(bool, Base, test_AABB, l, u);
        }
    };

}

}

#endif /* ECELL4_PYTHON_API_SHAPE_HPP */
