#ifndef ECELL4_PYTHON_API_TYPE_CASTER
#define ECELL4_PYTHON_API_TYPE_CASTER

#include <memory>
#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <pybind11/stl.h>
// #include <pybind11/stl_bind.h>

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace pybind11
{
    namespace detail
    {
        template<typename... Ts>
        struct type_caster<boost::variant<Ts...>> : variant_caster<boost::variant<Ts...>>
        {
        };

        template<typename T>
        struct type_caster<boost::optional<T>> : optional_caster<boost::optional<T>>
        {
        };

        template <>
        struct visit_helper<boost::variant> {
            template <typename... Args>
            static auto call(Args &&... args) -> decltype(boost::apply_visitor(args...))
            {
                return boost::apply_visitor(args...);
            }
        };
    }
}


#endif /* ECELL4_PYTHON_API_TYPE_CASTER */
