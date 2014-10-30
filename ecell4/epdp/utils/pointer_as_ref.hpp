#ifndef UTILS_POINTER_AS_REF_HPP
#define UTILS_POINTER_AS_REF_HPP

#include <boost/get_pointer.hpp>
#include <boost/call_traits.hpp>

template<typename T_, typename Tholder_ = T_*>
struct pointer_as_ref
{
    typedef T_ element_type;
    typedef Tholder_ holder_type;

    operator T_&() const { return *ptr_; }

    T_& operator=(element_type const& val) const
    {
        *ptr_ = val;
        return *ptr_;
    }

    T_* get() const
    {
        return ptr_;
    }

    void set(T_* ptr)
    {
        ptr_ = ptr;
    }

    operator bool() const
    {
        return !!ptr_;
    }

    explicit pointer_as_ref(typename boost::call_traits<holder_type>::param_type ptr)
        : ptr_(ptr) {}

    explicit pointer_as_ref(): ptr_(0) {}

private:
    holder_type ptr_;
};

namespace boost {

template<typename T, typename Tholder>
inline T* get_pointer(::pointer_as_ref<T, Tholder> const& p)
{
    return p.get();
}

} // namespace boost

#endif /* UTILS_POINTER_AS_REF_HPP */
