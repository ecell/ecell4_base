#include "Mesh.hpp"
#include "exceptions.hpp"

#include <vtkPoints.h>

namespace ecell4
{

MeshSurface::MeshSurface(const std::string filename, const Real3& edge_lengths)
    : filename_(filename), edge_lengths_(edge_lengths)
{
    {
        reader_ = vtkSmartPointer<vtkSTLReader>::New();
        reader_->SetFileName(filename_.c_str());
        reader_->Update();
        // tree_ = vtkSmartPointer<vtkOBBTree>::New();
        // tree_->SetDataSet(reader_->GetOutput());
        // tree_->BuildLocator();
    }

    {
        Real bounds[6];
        reader_->GetOutput()->GetBounds(bounds);

        const Real xratio(edge_lengths_[0] / (bounds[1] - bounds[0]));
        const Real yratio(edge_lengths_[1] / (bounds[3] - bounds[2]));
        const Real zratio(edge_lengths_[2] / (bounds[5] - bounds[4]));
        ratio_ = std::min(std::min(xratio, yratio), zratio);
        shift_ = Real3(-bounds[0], -bounds[2], -bounds[4]); //XXX: align origin
    }
}

MeshSurface::MeshSurface(const MeshSurface& rhs)
    : filename_(rhs.filename()), edge_lengths_(rhs.edge_lengths())
{
    {
        reader_ = vtkSmartPointer<vtkSTLReader>::New();
        reader_->SetFileName(filename_.c_str());
        reader_->Update();
        // tree_ = vtkSmartPointer<vtkOBBTree>::New();
        // tree_->SetDataSet(reader_->GetOutput());
        // tree_->BuildLocator();
    }

    {
        Real bounds[6];
        reader_->GetOutput()->GetBounds(bounds);

        const Real xratio(edge_lengths_[0] / (bounds[1] - bounds[0]));
        const Real yratio(edge_lengths_[1] / (bounds[3] - bounds[2]));
        const Real zratio(edge_lengths_[2] / (bounds[5] - bounds[4]));
        ratio_ = std::min(std::min(xratio, yratio), zratio);
        shift_ = Real3(-bounds[0], -bounds[2], -bounds[4]); //XXX: align origin
    }
}

Real MeshSurface::is_inside(const Real3& pos) const
{
    double lineP0[3];
    lineP0[0] = pos[0] / ratio_ - shift_[0];
    lineP0[1] = pos[1] / ratio_ - shift_[1];
    lineP0[2] = pos[2] / ratio_ - shift_[2];
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(lineP0);

    vtkSmartPointer<vtkPolyData> pointsPolydata = vtkSmartPointer<vtkPolyData>::New();
    pointsPolydata->SetPoints(points);
    vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints
        = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
    selectEnclosedPoints->SetInput(pointsPolydata);
    selectEnclosedPoints->SetSurface(reader_->GetOutput());
    selectEnclosedPoints->Update();
    return (selectEnclosedPoints->IsInside(0) ? 0.0 : inf);

    // double lineP0[3];
    // lineP0[0] = pos[0] / ratio_ - shift_[0];
    // lineP0[1] = pos[1] / ratio_ - shift_[1];
    // lineP0[2] = pos[2] / ratio_ - shift_[2];
    // double lineP1[3];
    // lineP1[0] = 0.0 / ratio_ - shift_[0];
    // lineP1[1] = pos[1] / ratio_ - shift_[1];
    // lineP1[2] = pos[2] / ratio_ - shift_[2];
    // vtkSmartPointer<vtkPoints> intersectPoints = vtkSmartPointer<vtkPoints>::New();
    // tree_->IntersectWithLine(lineP0, lineP1, intersectPoints, NULL);
    // return (intersectPoints->GetNumberOfPoints() % 2 == 1 ? 0.0 : inf);
}

Real3 MeshSurface::draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    throw NotImplemented("not implemented yet.");
}

bool MeshSurface::test_AABB(const Real3& l, const Real3& u) const
{
    throw NotImplemented("not implemented yet.");
}

} // ecell4
