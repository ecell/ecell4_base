#ifndef PEER_UTIL_REFERENCE_ACCESSOR_WRAPPER
#define PEER_UTIL_REFERENCE_ACCESSOR_WRAPPER

namespace peer { namespace util {

template<typename T_, typename Tval_,
        Tval_ const&(T_::*Vgetter_)() const,
        Tval_ &(T_::*Vsetter_)()>
struct reference_accessor_wrapper
{
    static Tval_ const& get(T_ const& impl)
    {
        return (impl.*Vgetter_)();
    }

    static void set(T_& impl, Tval_ const& v)
    {
        (impl.*Vsetter_)() = v;
    }
};

} } //namespace peer::util

#endif /* PEER_UTIL_REFERENCE_ACCESSOR_WRAPPER */
