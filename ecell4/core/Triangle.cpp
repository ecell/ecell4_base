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

Real3 closest_point_on_Triangle(const Real3& pos, const std::array<Real3, 3>& vertices)
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

// ----------------------------------------------------------------------------
// considering boundary condition

Real minmaxdist_sq(const Real3& lw, const Real3& up, const Real3& p) noexcept
{
    const auto sq = [](const Real x) noexcept -> Real {return x * x;};

    Real3 rm_sq(sq(lw[0] - p[0]), sq(lw[1] - p[1]), sq(lw[2] - p[2]));
    Real3 rM_sq(sq(up[0] - p[0]), sq(up[1] - p[1]), sq(up[2] - p[2]));

    using std::swap;
    if((up[0] + lw[0]) * 0.5 < p[0])
    {
        swap(rm_sq[0], rM_sq[0]);
    }
    if((up[1] + lw[1]) * 0.5 < p[1])
    {
        swap(rm_sq[1], rM_sq[1]);
    }
    if((up[2] + lw[2]) * 0.5 < p[2])
    {
        swap(rm_sq[2], rM_sq[2]);
    }

    const Real dx = rm_sq[0] + rM_sq[1] + rM_sq[2];
    const Real dy = rM_sq[0] + rm_sq[1] + rM_sq[2];
    const Real dz = rM_sq[0] + rM_sq[1] + rm_sq[2];
    return std::min(dx, std::min(dy, dz));
}

Real distance_sq_point_Triangle_impl(const Real3& pos, const Triangle& tri, const Boundary* b)
{
    const auto& vtxs = tri.vertices();

    // first, calculate the AABB of triangle
    constexpr Real inf = std::numeric_limits<Real>::infinity();
    Real3 lower, upper;
    lower[0] = std::min(vtxs[0][0], std::min(vtxs[1][0], vtxs[2][0]));
    lower[1] = std::min(vtxs[0][1], std::min(vtxs[1][1], vtxs[2][1]));
    lower[2] = std::min(vtxs[0][2], std::min(vtxs[1][2], vtxs[2][2]));
    upper[0] = std::max(vtxs[0][0], std::max(vtxs[1][0], vtxs[2][0]));
    upper[1] = std::max(vtxs[0][1], std::max(vtxs[1][1], vtxs[2][1]));
    upper[2] = std::max(vtxs[0][2], std::max(vtxs[1][2], vtxs[2][2]));

    const Real3 center = (lower + upper) * 0.5;
    const Real3 width  = (upper - lower) * 0.5;
    const Real3 edge   = b->edge_lengths();

    assert(0.0 <= width[0] && 0.0 <= width[1] && 0.0 <= width[2]);

    // transpose `pos` according to the center of the AABB
    const Real3 p1 = b->periodic_transpose(pos, center);
    const Real  d1 = length(closest_point_on_Triangle(p1, vtxs) - p1);
    const Real  D  = 2 * (length(width) + d1);

    // Here, D is the diameter of sphere that represents the region in which
    // point can be closer to the triangle than the original position, p1.
    // It means that if a periodic image of p1 exceeds this range, we don't
    // need to check the image.
    if(D < edge[0] && D < edge[1] && D < edge[2]) // likely
    {
        // we don't need to check any of periodic images. The minimum distance
        // between p1 and its periodic image is larger than the diameter of the
        // mindist-bounding sphere.
        return d1 * d1; // return square distance
    }

    // expand the AABB of Triangle by default mindist.
    // If periodic image exceeds this range along any axis,
    // we don't need to check it.
    lower[0] -= d1;
    lower[1] -= d1;
    lower[2] -= d1;

    upper[0] += d1;
    upper[1] += d1;
    upper[2] += d1;

    Real dist_sq = d1 * d1;

    // check all the possible transpose and find the minimum distance
    for(std::int32_t i_x=-1; i_x<=1; ++i_x)
    {
        const Real p_x = p1[0] + i_x * edge[0];
        if(p_x < lower[0] || upper[0] < p_x) {continue;}

        for(std::int32_t i_y=-1; i_y<=1; ++i_y)
        {
            const Real p_y = p1[1] + i_y * edge[1];
            if(p_y < lower[1] || upper[1] < p_y) {continue;}

            for(std::int32_t i_z=-1; i_z<=1; ++i_z)
            {
                const Real p_z = p1[2] + i_z * edge[2];
                if(p_z < lower[2] || upper[2] < p_z) {continue;}
                if(i_x == 0 && i_y == 0 && i_z == 0) {continue;}

                const Real3 p(p_x, p_y, p_z);
                dist_sq = std::min(dist_sq,
                    length_sq(closest_point_on_Triangle(p, vtxs) - p));
            }
        }
    }
    return dist_sq;
}

Real distance_sq_point_Triangle(const Real3& pos, const Triangle& tri,
                                const Boundary& b)
{
    return distance_sq_point_Triangle_impl(pos, tri, std::addressof(b));
}

Real distance_sq_point_Triangle(const Real3& pos, const Triangle& tri,
                                const std::unique_ptr<Boundary>& b)
{
    return distance_sq_point_Triangle_impl(pos, tri, b.get());
}
Real distance_sq_point_Triangle(const Real3& pos, const Triangle& tri,
                                const std::shared_ptr<Boundary>& b)
{
    return distance_sq_point_Triangle_impl(pos, tri, b.get());
}

}// ecell4
