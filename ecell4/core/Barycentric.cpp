#include "Barycentric.hpp"

namespace ecell4
{

Barycentric<Real>
to_barycentric(const Real3& pos, const Triangle& face)
{
    const Real3& a = face.vertex_at(0);
    const Real3& b = face.vertex_at(1);
    const Real3& c = face.vertex_at(2);
    const Real3  m = cross_product(face.edge_at(0), face.edge_at(2)) * (-1.);
    const Real   x = std::abs(m[0]);
    const Real   y = std::abs(m[1]);
    const Real   z = std::abs(m[2]);

    Real nu, nv, ood;
    if(x >= y && x >= z)
    {
        nu = detail::triangle_area_2D(pos[1], pos[2], b[1], b[2], c[1], c[2]);
        nv = detail::triangle_area_2D(pos[1], pos[2], c[1], c[2], a[1], a[2]);
        ood = 1.0 / m[0];
    }
    else if(y >= x && y >= z)
    {
        nu = detail::triangle_area_2D(pos[0], pos[2], b[0], b[2], c[0], c[2]);
        nv = detail::triangle_area_2D(pos[0], pos[2], c[0], c[2], a[0], a[2]);
        ood = 1.0 / -m[1];
    }
    else
    {
        nu = detail::triangle_area_2D(pos[0], pos[1], b[0], b[1], c[0], c[1]);
        nv = detail::triangle_area_2D(pos[0], pos[1], c[0], c[1], a[0], a[1]);
        ood = 1.0 / m[2];
    }
    Barycentric<Real> bary;
    bary[0] = nu * ood;
    bary[1] = nv * ood;
    bary[2] = 1.0 - bary[0] - bary[1];
    return bary;
}


}// ecell4
