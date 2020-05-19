#ifndef ECELL4_EGFRD_UTILS_FUN_COMPOSITION_HPP
#define ECELL4_EGFRD_UTILS_FUN_COMPOSITION_HPP
#include <utility>

namespace ecell4
{
namespace egfrd
{

namespace detail
{
template<typename F1, typename F2>
struct fun_compose_impl
{
    using result_type = typename F1::result_type;

    fun_compose_impl(const F1& f1, const F2& f2): f1_(f1), f2_(f2) {}

    template<typename ... Args>
    result_type operator()(Args&& ... args) const
    {
        return f1_( f2_(std::forward<Args>(args)...) );
    }

    template<typename ... Args>
    result_type operator()(Args&& ... args)
    {
        return f1_( f2_(std::forward<Args>(args)...) );
    }

private:
    F1 f1_;
    F2 f2_;
};
} // namespace detail

template<typename F1, typename F2>
inline detail::fun_compose_impl<F1, F2> fun_composition(const F1& f1, const F2& f2)
{
    return detail::fun_compose_impl<F1, F2>(f1, f2);
}

} // egfrd
} // ecell4
#endif /* FUN_COMPOSITION_HPP */
