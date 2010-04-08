#ifndef DOMAIN_UTILS_HPP
#define DOMAIN_UTILS_HPP

#include "exceptions.hpp"
#include "Single.hpp"
#include "Pair.hpp"
#include "FirstPassageGreensFunction.hpp"

template<typename Ttraits_>
struct DomainUtils
{
    typedef Ttraits_ traits_type;
    typedef typename traits_type::spherical_shell_type spherical_shell_type;
    typedef typename traits_type::cylindrical_shell_type cylindrical_shell_type;
    typedef typename traits_type::domain_type domain_type;
    typedef Single<Ttraits_, spherical_shell_type> spherical_single_type;
    typedef Single<Ttraits_, cylindrical_shell_type> cylindrical_single_type;
    typedef Pair<Ttraits_, spherical_shell_type> spherical_pair_type;
    typedef Pair<Ttraits_, cylindrical_shell_type> cylindrical_pair_type;

    static length_type calculate_mobility_radius(spherical_single_type const& dom)
    {
        return dom.shell().second.shape().radius() - dom.particle().second.shape().radius();
    }

    static length_type calculate_mobility_radius(cylindrical_single_type const& dom)
    {
        return dom.shell().second.shape().size() - dom.particle().second.shape().radius();
    }

    static length_type calculate_mobility_radius(domain_type const& dom)
    {
        {
            spherical_single_type const* x(
                dynamic_cast<spherical_single_type const*>(&dom));
            if (x)
            {
                return calculate_mobility_radius(*x);
            }
        }
        {
            cylindrical_single_type const* x(
                dynamic_cast<cylindrical_single_type const*>(&dom));
            if (x)
            {
                return calculate_mobility_radius(*x);
            }
        }
        throw unsupported("unsupported domain type");
    }

    static length_type get_shell_size(spherical_single_type const& dom)
    {
        return dom.shell().second.radius();
    }

    static length_type get_shell_size(cylindrical_single_type const& dom)
    {
        return dom.shell().second.radius();
    }

    static length_type get_shell_size(spherical_pair_type const& dom)
    {
        return dom.shell().second.radius();
    }

    static length_type get_shell_size(cylindrical_pair_type const& dom)
    {
        return dom.shell().second.radius();
    }

    static length_type get_shell_size(domain_type const& dom)
    {
        {
            spherical_single_type const* x(
                dynamic_cast<spherical_single_type const*>(&dom));
            if (x)
            {
                return get_shell_size(*x);
            }
        }
        {
            cylindrical_single_type const* x(
                dynamic_cast<cylindrical_single_type const*>(&dom));
            if (x)
            {
                return get_shell_size(*x);
            }
        }
        {
            spherical_pair_type const* x(
                dynamic_cast<spherical_pair_type const*>(&dom));
            if (x)
            {
                return get_shell_size(*x);
            }
        }
        {
            cylindrical_pair_type const* x(
                dynamic_cast<cylindrical_pair_type const*>(&dom));
            if (x)
            {
                return get_shell_size(*x);
            }
        }
        throw unsupported("unsupported domain type");
    }

    template<typename Tdom_>
    static FirstPassageGreensFunction get_com_greens_function(Tdom_ const& dom)
    {
        return FirstPassageGreensFunction(dom)
    }
};

#endif /* DOMAIN_UTILS_HPP */
