#ifndef BINDING_FCPAIR_REFERENCE_ACCESSOR_WRAPPER_HPP
#define BINDING_FCPAIR_REFERENCE_ACCESSOR_WRAPPER_HPP

namespace binding {

template<typename T_, typename Tval_,
        Tval_ const&(T_::*Vgetter_)() const,
        Tval_ &(T_::*Vsetter_)()>
struct fcpair_reference_accessor_wrapper
{
    static Tval_ const& get(T_ const& impl)
    {
        return (impl.*Vgetter_)();
    }

    static void set(T_& impl, Tval_ const& v)
    {
        (impl.*Vsetter_)().second = v.second;
    }
};

} // namespace binding

#endif /* BINDING_FCPAIR_REFERENCE_ACCESSOR_WRAPPER_HPP */
