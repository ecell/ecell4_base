#include "Triangle.hpp"

namespace ecell4
{

Triangle::Triangle()
{
    // do nothing
}

Triangle::Triangle(const boost::array<Real3, 3>& vertices)
{
    edges_[0] = vertices[1] - vertices[0];
    edges_[1] = vertices[2] - vertices[1];
    edges_[2] = vertices[0] - vertices[2];
    lengths_[0] = length(edges_[0]);
    lengths_[1] = length(edges_[1]);
    lengths_[2] = length(edges_[2]);
    angles_[0] = angle(edges_[0], edges_[2] * -1.0);
    angles_[1] = angle(edges_[1], edges_[0] * -1.0);
    angles_[2] = angle(edges_[2], edges_[1] * -1.0);
    normal_ = cross_product(edges_[0], edges_[2] * (-1));
    normal_ /= length(normal_);
}

Triangle::Triangle(const Real3& a, const Real3& b, const Real3& c)
{
    vertices_[0] = a;
    vertices_[1] = b;
    vertices_[2] = c;
    edges_[0] = vertices_[1] - vertices_[0];
    edges_[1] = vertices_[2] - vertices_[1];
    edges_[2] = vertices_[0] - vertices_[2];
    lengths_[0] = length(edges_[0]);
    lengths_[1] = length(edges_[1]);
    lengths_[2] = length(edges_[2]);
    angles_[0] = angle(edges_[0], edges_[2] * -1.0);
    angles_[1] = angle(edges_[1], edges_[0] * -1.0);
    angles_[2] = angle(edges_[2], edges_[1] * -1.0);
    normal_ = cross_product(edges_[0], edges_[2] * (-1));
    normal_ /= length(normal_);
}

}
