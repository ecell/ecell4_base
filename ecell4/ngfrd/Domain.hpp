#ifndef ECELL4_NGFRD_DOMAIN_HPP
#define ECELL4_NGFRD_DOMAIN_HPP
#include <ecell4/MultiDomain.hpp>

namespace ecell4
{
namespace ngfrd
{

struct Domain
{
public:
    enum class DomainKind : int // boost::variant::which returns an int.
    {
        Multi = 0,
    };

    using storage_type = boost::variant<MultiDomain>;

public:
    template<typename D>
    explicit Domain(D&& d): storage_(std::forward<D>(d)) {}

    DomainKind kind() const noexcept {return DomainKind(storage_.which());}

    bool is_multi() const noexcept {return this->kind() == DomainKind::Multi;}

    MultiDomain const& as_multi() const noexcept {return boost::get<MultiDomain>(storage_);}
    MultiDomain&       as_multi()       noexcept {return boost::get<MultiDomain>(storage_);}

    storage_type const& as_variant() const noexcept {return storage_;}
    storage_type&       as_variant()       noexcept {return storage_;}

private:

    storage_type storage_;
};

} //ngfrd
} //ecell4
#endif// ECELL4_NGFRD_DOMAIN_HPP
