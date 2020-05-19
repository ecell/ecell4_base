#ifndef ECELL4_EGFRD_UTILS_POINTER_AS_REF_HPP
#define ECELL4_EGFRD_UTILS_POINTER_AS_REF_HPP

#include <boost/get_pointer.hpp>

namespace ecell4
{
namespace egfrd
{

template<typename T_>
struct pointer_as_ref
{
    typedef T_ element_type;
    typedef T_* holder_type;

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

    explicit pointer_as_ref(holder_type ptr)
        : ptr_(ptr) {}

    explicit pointer_as_ref(): ptr_(0) {}

private:
    holder_type ptr_;
};

} // egfrd
} // ecell4

namespace boost {

template<typename T>
inline T* get_pointer(ecell4::egfrd::pointer_as_ref<T> const& p)
{
    return p.get();
}

} // namespace boost

#endif /* UTILS_POINTER_AS_REF_HPP */
