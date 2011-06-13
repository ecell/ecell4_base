#ifndef SERIAL_ID_GENERATOR_HPP
#define SERIAL_ID_GENERATOR_HPP

#include <functional>
#include <boost/type_traits/is_integral.hpp>

namespace detail {

template<bool Vis_integral, typename Tid_>
struct identifier_lot_helper
{
    typedef typename Tid_::lot_type type;
};

template<typename Tid_>
struct identifier_lot_helper<true, Tid_>
{
    typedef Tid_ type;
};

template<bool Vis_integral, typename Tid_>
struct identifier_lot_adder_helper: public std::binary_function<
        Tid_, typename Tid_::lot_type, Tid_>
{
    Tid_ operator()(Tid_ const& lhs, typename Tid_::lot_type rhs)
    {
        return lhs.lot_add(rhs);
    }
};

template<typename Tid_>
struct identifier_lot_adder_helper<true, Tid_>: public std::binary_function<
        Tid_, Tid_, Tid_>
{
    Tid_ operator()(Tid_ const& lhs, Tid_ const& rhs)
    {
        return lhs + rhs;
    }
};

template<bool Vis_integral, typename Tid_>
struct identifier_lot_advancer_helper: public std::binary_function<
        Tid_&, typename Tid_::lot_type, Tid_&>
{
    Tid_& operator()(Tid_& lhs, typename Tid_::lot_type const& rhs)
    {
        lhs.lot_advance(rhs);
        return lhs;
    }
};

template<typename Tid_>
struct identifier_lot_advancer_helper<true, Tid_>: public std::binary_function<
        Tid_&, Tid_, Tid_&>
{
    Tid_& operator()(Tid_& lhs, Tid_ const& rhs)
    {
        lhs += rhs;
        return lhs;
    }
};

template<bool Vis_integral, typename Tid_>
struct identifier_lot_retracer_helper: public std::binary_function<
        Tid_&, typename Tid_::lot_type, Tid_&>
{
    Tid_& operator()(Tid_& lhs, typename Tid_::lot_type const& rhs)
    {
        lhs.lot_retrace(rhs);
        return lhs;
    }
};

template<typename Tid_>
struct identifier_lot_retracer_helper<true, Tid_>: public std::binary_function<
        Tid_, Tid_, Tid_>
{
    Tid_& operator()(Tid_& lhs, Tid_ const& rhs)
    {
        lhs -= rhs;
        return lhs;
    }
};

template<bool Vis_integral, typename Tid_>
struct identifier_lot_retriever_helper: public std::binary_function<
        Tid_&, typename Tid_::lot_type, Tid_&>
{
    typename identifier_lot_helper<Vis_integral, Tid_>::type const& operator()(Tid_ const& lhs)
    {
        return lhs.lot();
    }
};

template<typename Tid_>
struct identifier_lot_retriever_helper<true, Tid_>: public std::binary_function<
        Tid_, Tid_, Tid_>
{
    typename identifier_lot_helper<true, Tid_>::type& operator()(Tid_& lhs)
    {
        return lhs;
    }
};

template<bool Vis_integral, typename Tid_>
struct identifier_serial_helper
{
    typedef typename Tid_::serial_type type;
};

template<typename Tid_>
struct identifier_serial_helper<true, Tid_>
{
    typedef Tid_ type;
};
template<bool Vis_integral, typename Tid_>
struct identifier_serial_advancer_helper: public std::binary_function<
        Tid_&, typename Tid_::serial_type, Tid_&>
{
    Tid_& operator()(Tid_& lhs, typename Tid_::serial_type const& rhs)
    {
        lhs.serial_advance(rhs);
        return lhs;
    }
};

template<typename Tid_>
struct identifier_serial_advancer_helper<true, Tid_>: public std::binary_function<
        Tid_&, Tid_, Tid_&>
{
    Tid_& operator()(Tid_& lhs, Tid_ const& rhs)
    {
        lhs += rhs;
        return lhs;
    }
};

template<bool Vis_integral, typename Tid_>
struct identifier_serial_retracer_helper: public std::binary_function<
        Tid_&, typename Tid_::serial_type, Tid_&>
{
    Tid_& operator()(Tid_& lhs, typename Tid_::serial_type const& rhs)
    {
        lhs.serial_retrace(rhs);
        return lhs;
    }
};

template<typename Tid_>
struct identifier_serial_retracer_helper<true, Tid_>: public std::binary_function<
        Tid_, Tid_, Tid_>
{
    Tid_& operator()(Tid_& lhs, Tid_ const& rhs)
    {
        lhs -= rhs;
        return lhs;
    }
};

template<bool Vis_integral, typename Tid_>
struct identifier_serial_retriever_helper: public std::binary_function<
        Tid_&, typename Tid_::serial_type, Tid_&>
{
    typename identifier_serial_helper<Vis_integral, Tid_>::type const& operator()(Tid_ const& lhs)
    {
        return lhs.serial();
    }
};

template<typename Tid_>
struct identifier_serial_retriever_helper<true, Tid_>: public std::binary_function<
        Tid_, Tid_, Tid_>
{
    typename identifier_serial_helper<true, Tid_>::type& operator()(Tid_& lhs)
    {
        return lhs;
    }
};

} // namespace detail

template<typename Tid_>
struct identifier_lot: public detail::identifier_lot_helper<
        boost::is_integral<Tid_>::value, Tid_>
{
};

template<typename Tid_>
struct identifier_lot_adder
    : public detail::identifier_lot_adder_helper<
    boost::is_integral<Tid_>::value, Tid_>
{
};

template<typename Tid_>
struct identifier_lot_advancer
        : public detail::identifier_lot_advancer_helper<
            boost::is_integral<Tid_>::value, Tid_>
{
};

template<typename Tid_>
struct identifier_lot_retracer
        : public detail::identifier_lot_retracer_helper<
            boost::is_integral<Tid_>::value, Tid_>
{
};

template<typename Tid_>
struct identifier_lot_retriever
        : public detail::identifier_lot_retriever_helper<
            boost::is_integral<Tid_>::value, Tid_>
{
};

template<typename Tid_>
struct identifier_serial: public detail::identifier_serial_helper<
        boost::is_integral<Tid_>::value, Tid_>
{
};

template<typename Tid_>
struct identifier_serial_advancer
        : public detail::identifier_serial_advancer_helper<
            boost::is_integral<Tid_>::value, Tid_>
{
};

template<typename Tid_>
struct identifier_serial_retracer
        : public detail::identifier_serial_retracer_helper<
            boost::is_integral<Tid_>::value, Tid_>
{
};

template<typename Tid_>
struct identifier_serial_retriever
        : public detail::identifier_serial_retriever_helper<
            boost::is_integral<Tid_>::value, Tid_>
{
};

template<typename Tid_>
Tid_ lot_add(Tid_ const& lhs, typename identifier_lot<Tid_>::type const& rhs)
{
    return identifier_lot_adder<Tid_>()(lhs, rhs);
}

template<typename Tid_>
Tid_& lot_advance(Tid_& lhs, typename identifier_lot<Tid_>::type const& rhs)

{
    return identifier_lot_advancer<Tid_>()(lhs, rhs);
}

template<typename Tid_>
Tid_& lot_retrace(Tid_& lhs, typename identifier_lot<Tid_>::type const& rhs)
{
    return identifier_lot_retracer<Tid_>()(lhs, rhs);
}

template<typename Tid_>
typename identifier_lot<Tid_>::type lot(Tid_& lhs)
{
    return identifier_lot_retriever<Tid_>()(lhs);
}

template<typename Tid_>
Tid_& serial_advance(Tid_& lhs, typename identifier_serial<Tid_>::type const& rhs)
{
    return identifier_serial_advancer<Tid_>()(lhs, rhs);
}

template<typename Tid_>
Tid_& serial_retrace(Tid_& lhs, typename identifier_serial<Tid_>::type const& rhs)
{
    return identifier_serial_retracer<Tid_>()(lhs, rhs);
}

template<typename Tid_>
typename identifier_serial<Tid_>::type serial(Tid_& lhs)
{
    return identifier_serial_retriever<Tid_>()(lhs);
}

template<typename Tid_>
struct SerialIDGenerator
{
    typedef Tid_ identifier_type;
    typedef typename identifier_lot<identifier_type>::type lot_type;

    SerialIDGenerator(lot_type const& lot = lot_type())
        : next_(lot_add(identifier_type(), lot))
    {
    }

    identifier_type operator()()
    {
        return serial_advance(next_, 1);
    }

private:
    identifier_type next_;
};

#endif /* SERIAL_ID_GENERATOR_HPP */
