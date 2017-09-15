#ifndef ECELL4_MESH_HPP
#define ECELL4_MESH_HPP

#include <ecell4/core/config.h>
#include "Shape.hpp"

#ifdef HAVE_VTK
#include <vtkSmartPointer.h>
#include <vtkCellLocator.h>
#include <vtkSTLReader.h>
#include <vtkOBBTree.h>
#include <vtkPolyData.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkPoints.h>
#include <vtkTriangle.h>
#endif

namespace ecell4
{

struct MeshSurface
    : public Shape
{
public:

    MeshSurface(const std::string filename, const Real3& edge_lengths);
    MeshSurface(const MeshSurface& rhs);

    std::string filename() const
    {
        return filename_;
    }

    Real3 edge_lengths() const
    {
        return edge_lengths_;
    }

    virtual dimension_kind dimension() const
    {
        return TWO;
    }

    virtual Real is_inside(const Real3& pos) const;
    virtual Real3 draw_position(
        boost::shared_ptr<RandomNumberGenerator>& rng) const;
    virtual bool test_AABB(const Real3& l, const Real3& u) const;

#ifdef HAVE_VTK
    virtual void bounding_box(
        const Real3& edge_lengths, Real3& lower, Real3& upper) const
    {
        double bounds[6];
        reader_->GetOutput()->GetBounds(bounds);

        const Real xlim(ratio_ * (bounds[1] - bounds[0]));
        const Real ylim(ratio_ * (bounds[3] - bounds[2]));
        const Real zlim(ratio_ * (bounds[5] - bounds[4]));

        lower = Real3(0.0, 0.0, 0.0);
        upper = Real3(xlim, ylim, zlim);
    }
#endif

protected:

    std::string filename_;
    Real3 edge_lengths_;
    Real ratio_;
    Real3 shift_;

#ifdef HAVE_VTK
    vtkSmartPointer<vtkSTLReader> reader_;
    // vtkSmartPointer<vtkOBBTree> tree_;
#endif
};

} // ecell4

#endif /* ECELL4_MESH_HPP */
