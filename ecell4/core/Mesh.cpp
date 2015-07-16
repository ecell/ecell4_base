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
        tree_ = vtkSmartPointer<vtkOBBTree>::New();
        tree_->SetDataSet(reader_->GetOutput());
        tree_->BuildLocator();
    }

    {
        Real bounds[6];
        reader_->GetOutput()->GetBounds(bounds);
        ratio_ = std::min(std::min(edge_lengths_[0] / (bounds[1] - bounds[0]), edge_lengths_[1] / (bounds[3] - bounds[2])), edge_lengths_[2] / (bounds[5] - bounds[4]));
        shift_ = Real3(-bounds[0], -bounds[2], -bounds[4]); //XXX: align origin

        // std::cout  << "xmin: " << bounds[0] << " "
        //     << "xmax: " << bounds[1] << std::endl
        //     << "ymin: " << bounds[2] << " "
        //     << "ymax: " << bounds[3] << std::endl
        //     << "zmin: " << bounds[4] << " "
        //     << "zmax: " << bounds[5] << std::endl;
    }
}

MeshSurface::MeshSurface(const MeshSurface& rhs)
    : filename_(rhs.filename()), edge_lengths_(rhs.edge_lengths())
{
    {
        reader_ = vtkSmartPointer<vtkSTLReader>::New();
        reader_->SetFileName(filename_.c_str());
        reader_->Update();
        tree_ = vtkSmartPointer<vtkOBBTree>::New();
        tree_->SetDataSet(reader_->GetOutput());
        tree_->BuildLocator();
    }

    {
        Real bounds[6];
        reader_->GetOutput()->GetBounds(bounds);
        ratio_ = std::min(std::min(edge_lengths_[0] / (bounds[1] - bounds[0]), edge_lengths_[1] / (bounds[3] - bounds[2])), edge_lengths_[2] / (bounds[5] - bounds[4]));
        shift_ = Real3(-bounds[0], -bounds[2], -bounds[4]); //XXX: align origin
    }
}

Real MeshSurface::is_inside(const Real3& pos) const
{
    double lineP0[3];
    lineP0[0] = pos[0] / ratio_ - shift_[0];
    lineP0[1] = pos[1] / ratio_ - shift_[1];
    lineP0[2] = pos[2] / ratio_ - shift_[2];
    double lineP1[3];
    lineP1[0] = 0.0 / ratio_ - shift_[0];
    lineP1[1] = pos[1] / ratio_ - shift_[1];
    lineP1[2] = pos[2] / ratio_ - shift_[2];
    // std::cout << "(" << lineP0[0] << ", " << lineP0[1] << ", " << lineP0[2] << ")" << std::endl;
    vtkSmartPointer<vtkPoints> intersectPoints = vtkSmartPointer<vtkPoints>::New();
    tree_->IntersectWithLine(lineP0, lineP1, intersectPoints, NULL);
    // std::cout << "=> " << intersectPoints->GetNumberOfPoints() << std::endl;
    return (intersectPoints->GetNumberOfPoints() % 2 == 1 ? 0.0 : inf);
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
