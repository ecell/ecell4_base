
#ifndef CYLINDER_HPP
#define CYLINDER_HPP

#include <ostream>
#include "Vector3.hpp"
#include <cmath>

// Todo. Make sure cylinder is never larger than 1 cellsize or something.  
template<typename T_>
class Cylinder
{
public:
    typedef T_ value_type;
    typedef Vector3<T_> position_type;
    typedef T_ length_type;

public:
    Cylinder()
        : position_(), radius_(0), orientationZ_(), size_(0) {}

    Cylinder(const position_type& position, const length_type& radius, const position_type& orientationZ, const length_type& size )
        : position_(position), radius_(radius), orientationZ_(orientationZ), size_(size)
    { }

private:
    length_type calculateDistanceToSelfWithNewOrigin( position_type pos, position_type newOrigin )
    {
        /* First compute the (z,r) components of pos in a coordinate system 
         * defined by the vectors unitR and orientationZ, where unitR is
         * choosen such that unitR and orientationZ define a plane in which
         * pos lies. */
        position_type posVector = pos - newOrigin;

        length_type z = dot_product(posVector, orientationZ_); // Can be <0.
        position_type posVectorZ = multiply(orientationZ_, z);

        position_type posVectorR = posVector - posVectorZ;
        length_type r = length(posVectorR);


        /* Then compute distance to cylinder. */
        length_type dz = fabs(z) - size_;
        length_type dr = r - radius_;
        length_type distance;
        if(dz > 0){
            // pos is (either) to the right or to the left of the cylinder.
            if(r > radius_){
                // Compute distance to edge.
                distance = std::sqrt( dz*dz + dr*dr );
            }
            else{
                distance = dz;
            }
        }
        else{
            if(dr > radius_){
                // pos is somewhere 'parellel' to the cylinder.
                distance = dr;
            }
            else{
                // Inside cylinder. 
                distance = std::max(dr, dz);
            }
        }
        return distance;
    }

public:
    length_type calculateDistanceToSelf( position_type pos )
    {
        return calculateDistanceToSelfWithNewOrigin( pos, position_ );
    }

    length_type calculateDistanceToSelfWithOffset( position_type pos, position_type offset )
    {
        // Because this cylinder is on the other side of one of the periodic 
        // boundaries compared to pos, add offset (for example (-L, 0, 0) to 
        // position before calculating distance between this cylinder and pos.  
        return calculateDistanceToSelfWithNewOrigin( pos, position_ + offset );
    }

    bool operator==(const Cylinder& rhs) const
    {
        return position_ == rhs.position_ && radius_ == rhs.radius_ && orientationZ_ == rhs.orientationZ_ && size_ == rhs.size_;
    }

    bool operator!=(const Cylinder& rhs) const
    {
        return !operator==(rhs);
    }

    position_type const& position() const
    {
        return position_;
    }

    position_type& position()
    {
        return position_;
    }

    length_type const& radius() const
    {
        return radius_;
    }

    length_type& radius()
    {
        return radius_;
    }

    position_type const& orientationZ() const
    {
        return orientationZ_;
    }

    position_type& orientationZ()
    {
        return orientationZ_;
    }

    length_type const& size() const
    {
        return size_;
    }

    length_type& size()
    {
        return size_;
    }

private:
    // Todo. These need to be public for comparison operator==?
    position_type position_; // centre.
    length_type radius_;
    position_type orientationZ_; // should be normalized.
    length_type size_; // half length.
};

template<typename Tstrm_, typename T_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const Cylinder<T_>& v)
{
    strm << "{" << v.position() <<  ", " << v.radius() << ", " << v.orientationZ() << ", " << v.size() << "}";
    return strm;
}

#endif /* CYLINDER_HPP */
