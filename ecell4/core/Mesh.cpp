#include <numeric>
#include "Mesh.hpp"
#include "exceptions.hpp"

namespace ecell4
{

MeshSurface::MeshSurface(const std::string filename, const Real3& edge_lengths)
    : filename_(filename), edge_lengths_(edge_lengths)
{
#ifdef HAVE_VTK
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
#endif
}

MeshSurface::MeshSurface(const MeshSurface& rhs)
    : filename_(rhs.filename()), edge_lengths_(rhs.edge_lengths())
{
#ifdef HAVE_VTK
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
#endif
}

Real MeshSurface::is_inside(const Real3& pos) const
{
#ifdef HAVE_VTK
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
#else
    throw NotImplemented("not implemented yet.");
#endif
}

Real3 MeshSurface::draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const
{
#ifdef HAVE_VTK
    vtkPolyData* polydata = reader_->GetOutput();
    std::vector<double> areas(polydata->GetNumberOfCells());
    for (vtkIdType i(0); i < polydata->GetNumberOfCells(); i++)
    {
        vtkCell* cell = polydata->GetCell(i);
        vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(cell);
        double p0[3];
        double p1[3];
        double p2[3];
        triangle->GetPoints()->GetPoint(0, p0);
        triangle->GetPoints()->GetPoint(1, p1);
        triangle->GetPoints()->GetPoint(2, p2);
        const double area = vtkTriangle::TriangleArea(p0, p1, p2);
        // std::cout << "p0: " << p0[0] << " " << p0[1] << " " << p0[2] << std::endl;
        // std::cout << "p1: " << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
        // std::cout << "p2: " << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
        // std::cout << "area of triangle " << i << ": " << area << std::endl;
        areas[i] = area;
    }
    const double rnd = rng->uniform(0.0, std::accumulate(areas.begin(), areas.end(), 0.0));
    double totarea = 0.0;
    for (vtkIdType i(0); i < polydata->GetNumberOfCells(); i++)
    {
        totarea += areas[i];
        if (rnd < totarea)
        {
            vtkCell* cell = polydata->GetCell(i);
            vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(cell);
            double p0[3];
            double p1[3];
            double p2[3];
            triangle->GetPoints()->GetPoint(0, p0);
            triangle->GetPoints()->GetPoint(1, p1);
            triangle->GetPoints()->GetPoint(2, p2);
            const Real3 P0(p0[0], p0[1], p0[2]);
            const Real3 P1(p1[0], p1[1], p1[2]);
            const Real3 P2(p2[0], p2[1], p2[2]);
            const Real p(rng->uniform(0.0, 1.0)), q(rng->uniform(0.0, 1.0 - p));
            return (((P1 - P0) * p + (P2 - P0) * q + P0) + shift_) * ratio_;
        }
    }
    throw IllegalState("Never reach here.");
#else
    throw NotImplemented("not implemented yet.");
#endif
}

bool MeshSurface::test_AABB(const Real3& l, const Real3& u) const
{
    throw NotImplemented("not implemented yet.");
}

} // ecell4
