#ifndef SHAPE_HPP
#define SHAPE_HPP

template<typename Tobj_>
inline typename Tobj_::shape_type const& shape(Tobj_ const& obj)
{
    return obj.shape();
}

template<typename Tobj_>
inline typename Tobj_::shape_type& shape(Tobj_& obj)
{
    return obj.shape();
}

template<typename Tshape_>
inline Tshape_ offset(Tshape_ const& shape, typename Tshape_::position_type off)
{
    Tshape_ retval(shape);
    retval.position() += off;
    return retval;
}

#endif /* SHAPE_HPP */
