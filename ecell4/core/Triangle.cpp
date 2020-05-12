#include "Triangle.hpp"

namespace ecell4
{

Triangle::Triangle()
{
    // do nothing
}

Triangle::Triangle(const std::array<Real3, 3>& vertices)
{
    vertices_[0] = vertices[0];
    vertices_[1] = vertices[1];
    vertices_[2] = vertices[2];
    edges_[0] = this->vertices_[1] - this->vertices_[0];
    edges_[1] = this->vertices_[2] - this->vertices_[1];
    edges_[2] = this->vertices_[0] - this->vertices_[2];
    lengths_[0] = length(edges_[0]);
    lengths_[1] = length(edges_[1]);
    lengths_[2] = length(edges_[2]);
    angles_[0] = calc_angle(edges_[0], edges_[2] * -1.0);
    angles_[1] = calc_angle(edges_[1], edges_[0] * -1.0);
    angles_[2] = calc_angle(edges_[2], edges_[1] * -1.0);
    normal_ = cross_product(edges_[0], edges_[2] * (-1));
    normal_ /= length(normal_);
}

Triangle::Triangle(const TriangleView& tv)
{
    vertices_[0] = tv.vertices(0);
    vertices_[1] = tv.vertices(1);
    vertices_[2] = tv.vertices(2);
    edges_[0] = this->vertices_[1] - this->vertices_[0];
    edges_[1] = this->vertices_[2] - this->vertices_[1];
    edges_[2] = this->vertices_[0] - this->vertices_[2];
    lengths_[0] = length(edges_[0]);
    lengths_[1] = length(edges_[1]);
    lengths_[2] = length(edges_[2]);
    angles_[0] = calc_angle(edges_[0], edges_[2] * -1.0);
    angles_[1] = calc_angle(edges_[1], edges_[0] * -1.0);
    angles_[2] = calc_angle(edges_[2], edges_[1] * -1.0);
    normal_ = cross_product(edges_[0], edges_[2] * (-1));
    normal_ /= length(normal_);
}

Triangle::Triangle(const TriangleConstView& tcv)
{
    vertices_[0] = tcv.vertices(0);
    vertices_[1] = tcv.vertices(1);
    vertices_[2] = tcv.vertices(2);
    edges_[0] = this->vertices_[1] - this->vertices_[0];
    edges_[1] = this->vertices_[2] - this->vertices_[1];
    edges_[2] = this->vertices_[0] - this->vertices_[2];
    lengths_[0] = length(edges_[0]);
    lengths_[1] = length(edges_[1]);
    lengths_[2] = length(edges_[2]);
    angles_[0] = calc_angle(edges_[0], edges_[2] * -1.0);
    angles_[1] = calc_angle(edges_[1], edges_[0] * -1.0);
    angles_[2] = calc_angle(edges_[2], edges_[1] * -1.0);
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
    angles_[0] = calc_angle(edges_[0], edges_[2] * -1.0);
    angles_[1] = calc_angle(edges_[1], edges_[0] * -1.0);
    angles_[2] = calc_angle(edges_[2], edges_[1] * -1.0);
    normal_ = cross_product(edges_[0], edges_[2] * (-1));
    normal_ /= length(normal_);
}

Triangle& Triangle::operator=(const Triangle& rhs)
{
    vertices_ = rhs.vertices_;
    edges_    = rhs.edges_;
    lengths_  = rhs.lengths_;
    angles_   = rhs.angles_;
    normal_   = rhs.normal_;
    return *this;
}

Triangle& Triangle::operator=(const TriangleView& tv)
{
    vertices_[0] = tv.vertices(0);
    vertices_[1] = tv.vertices(1);
    vertices_[2] = tv.vertices(2);
    edges_[0] = this->vertices_[1] - this->vertices_[0];
    edges_[1] = this->vertices_[2] - this->vertices_[1];
    edges_[2] = this->vertices_[0] - this->vertices_[2];
    lengths_[0] = length(edges_[0]);
    lengths_[1] = length(edges_[1]);
    lengths_[2] = length(edges_[2]);
    angles_[0] = calc_angle(edges_[0], edges_[2] * -1.0);
    angles_[1] = calc_angle(edges_[1], edges_[0] * -1.0);
    angles_[2] = calc_angle(edges_[2], edges_[1] * -1.0);
    normal_ = cross_product(edges_[0], edges_[2] * (-1));
    normal_ /= length(normal_);
    return *this;
}

Triangle& Triangle::operator=(const TriangleConstView& tcv)
{
    vertices_[0] = tcv.vertices(0);
    vertices_[1] = tcv.vertices(1);
    vertices_[2] = tcv.vertices(2);
    edges_[0] = this->vertices_[1] - this->vertices_[0];
    edges_[1] = this->vertices_[2] - this->vertices_[1];
    edges_[2] = this->vertices_[0] - this->vertices_[2];
    lengths_[0] = length(edges_[0]);
    lengths_[1] = length(edges_[1]);
    lengths_[2] = length(edges_[2]);
    angles_[0] = calc_angle(edges_[0], edges_[2] * -1.0);
    angles_[1] = calc_angle(edges_[1], edges_[0] * -1.0);
    angles_[2] = calc_angle(edges_[2], edges_[1] * -1.0);
    normal_ = cross_product(edges_[0], edges_[2] * (-1));
    normal_ /= length(normal_);
    return *this;
}

}// ecell4
