#include "Triangle.hpp"

namespace ecell4
{

Triangle::Triangle()
{
    // do nothing
}

Triangle::Triangle(const boost::array<Real3, 3>& vertices)
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

Real3 closest_point_on_Triangle(const Real3& pos, const boost::array<Real3, 3>& vertices)
{
    // this implementation is from Real-Time Collision Detection by Christer Ericson,
    // published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.
    // pp.141-142

    const Real3 a = vertices[0];
    const Real3 b = vertices[1];
    const Real3 c = vertices[2];

    const Real3 ab = b - a;
    const Real3 ac = c - a;
    const Real3 ap = pos - a;
    const Real  d1 = dot_product(ab, ap);
    const Real  d2 = dot_product(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0)
    {
        return a;
    }

    const Real3 bp = pos - b;
    const Real  d3 = dot_product(ab, bp);
    const Real  d4 = dot_product(ac, bp);
    if (d3 >= 0.0 && d4 <= d3)
    {
        return b;
    }

    const Real vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        Real v = d1 / (d1 - d3);
        return a + ab * v;
    }

    const Real3 cp = pos - c;
    const Real  d5 = dot_product(ab, cp);
    const Real  d6 = dot_product(ac, cp);
    if (d6 >= 0.0 && d5 <= d6)
    {
        return c;
    }

    const Real vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const Real w = d2 / (d2 - d6);
        return a + ac * w;
    }

    const Real va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const Real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + (c - b) * w;
    }

    const Real denom = 1.0 / (va + vb + vc);
    const Real v = vb * denom;
    const Real w = vc * denom;
    return a + ab * v + ac * w;
}

Real distance_sq_point_Triangle(const Real3& pos, const Triangle& tri)
{
    return length_sq(closest_point_on_Triangle(pos, tri.vertices()) - pos);
}

}// ecell4
