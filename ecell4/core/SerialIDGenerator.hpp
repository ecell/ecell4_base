#ifndef ECELL4_SERIAL_ID_GENERATOR_HPP
#define ECELL4_SERIAL_ID_GENERATOR_HPP

#include <functional>
#include <boost/type_traits/is_integral.hpp>
#include <boost/scoped_ptr.hpp>

#include <ecell4/core/config.h>

#ifdef WITH_HDF5
#include <hdf5.h>
#include <H5Cpp.h>
#endif


namespace ecell4
{

namespace detail
{

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
struct identifier_lot_adder_helper
    : public std::binary_function<Tid_, typename Tid_::lot_type, Tid_>
{
    Tid_ operator()(const Tid_& lhs, typename Tid_::lot_type rhs)
    {
        return lhs.lot_add(rhs);
    }
};

template<typename Tid_>
struct identifier_lot_adder_helper<true, Tid_>
    : public std::binary_function<Tid_, Tid_, Tid_>
{
    Tid_ operator()(const Tid_& lhs, const Tid_& rhs)
    {
        return lhs + rhs;
    }
};

template<bool Vis_integral, typename Tid_>
struct identifier_lot_advancer_helper
    : public std::binary_function<Tid_&, typename Tid_::lot_type, Tid_&>
{
    Tid_& operator()(Tid_& lhs, const typename Tid_::lot_type& rhs)
    {
        lhs.lot_advance(rhs);
        return lhs;
    }
};

template<typename Tid_>
struct identifier_lot_advancer_helper<true, Tid_>
    : public std::binary_function<Tid_&, Tid_, Tid_&>
{
    Tid_& operator()(Tid_& lhs, const Tid_& rhs)
    {
        lhs += rhs;
        return lhs;
    }
};

template<bool Vis_integral, typename Tid_>
struct identifier_lot_retracer_helper
    : public std::binary_function<Tid_&, typename Tid_::lot_type, Tid_&>
{
    Tid_& operator()(Tid_& lhs, const typename Tid_::lot_type& rhs)
    {
        lhs.lot_retrace(rhs);
        return lhs;
    }
};

template<typename Tid_>
struct identifier_lot_retracer_helper<true, Tid_>
    : public std::binary_function<Tid_, Tid_, Tid_>
{
    Tid_& operator()(Tid_& lhs, const Tid_& rhs)
    {
        lhs -= rhs;
        return lhs;
    }
};

template<bool Vis_integral, typename Tid_>
struct identifier_lot_retriever_helper
    : public std::binary_function<Tid_&, typename Tid_::lot_type, Tid_&>
{
    const typename identifier_lot_helper<Vis_integral, Tid_>::type&
    operator()(const Tid_& lhs)
    {
        return lhs.lot();
    }
};

template<typename Tid_>
struct identifier_lot_retriever_helper<true, Tid_>
    : public std::binary_function<Tid_, Tid_, Tid_>
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
struct identifier_serial_advancer_helper
    : public std::binary_function<Tid_&, typename Tid_::serial_type, Tid_&>
{
    Tid_& operator()(Tid_& lhs, const typename Tid_::serial_type& rhs)
    {
        lhs.serial_advance(rhs);
        return lhs;
    }
};

template<typename Tid_>
struct identifier_serial_advancer_helper<true, Tid_>
    : public std::binary_function<Tid_&, Tid_, Tid_&>
{
    Tid_& operator()(Tid_& lhs, const Tid_& rhs)
    {
        lhs += rhs;
        return lhs;
    }
};

template<bool Vis_integral, typename Tid_>
struct identifier_serial_retracer_helper
    : public std::binary_function<Tid_&, typename Tid_::serial_type, Tid_&>
{
    Tid_& operator()(Tid_& lhs, const typename Tid_::serial_type& rhs)
    {
        lhs.serial_retrace(rhs);
        return lhs;
    }
};

template<typename Tid_>
struct identifier_serial_retracer_helper<true, Tid_>
    : public std::binary_function<Tid_, Tid_, Tid_>
{
    Tid_& operator()(Tid_& lhs, const Tid_& rhs)
    {
        lhs -= rhs;
        return lhs;
    }
};

template<bool Vis_integral, typename Tid_>
struct identifier_serial_retriever_helper
    : public std::binary_function<Tid_&, typename Tid_::serial_type, Tid_&>
{
    const typename identifier_serial_helper<Vis_integral, Tid_>::type&
    operator()(const Tid_& lhs)
    {
        return lhs.serial();
    }
};

template<typename Tid_>
struct identifier_serial_retriever_helper<true, Tid_>
    : public std::binary_function<Tid_, Tid_, Tid_>
{
    typename identifier_serial_helper<true, Tid_>::type& operator()(Tid_& lhs)
    {
        return lhs;
    }
};

} // namespace detail

template<typename Tid_>
struct identifier_lot
    : public detail::identifier_lot_helper<boost::is_integral<Tid_>::value, Tid_>
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
struct identifier_serial
    : public detail::identifier_serial_helper<
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
Tid_ lot_add(const Tid_& lhs, const typename identifier_lot<Tid_>::type& rhs)
{
    return identifier_lot_adder<Tid_>()(lhs, rhs);
}

template<typename Tid_>
Tid_& lot_advance(Tid_& lhs, const typename identifier_lot<Tid_>::type& rhs)

{
    return identifier_lot_advancer<Tid_>()(lhs, rhs);
}

template<typename Tid_>
Tid_& lot_retrace(Tid_& lhs, const typename identifier_lot<Tid_>::type& rhs)
{
    return identifier_lot_retracer<Tid_>()(lhs, rhs);
}

template<typename Tid_>
typename identifier_lot<Tid_>::type lot(Tid_& lhs)
{
    return identifier_lot_retriever<Tid_>()(lhs);
}

template<typename Tid_>
Tid_& serial_advance(
    Tid_& lhs, const typename identifier_serial<Tid_>::type& rhs)
{
    return identifier_serial_advancer<Tid_>()(lhs, rhs);
}

template<typename Tid_>
Tid_& serial_retrace(
    Tid_& lhs, const typename identifier_serial<Tid_>::type& rhs)
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
public:

    typedef Tid_ identifier_type;
    typedef typename identifier_lot<identifier_type>::type lot_type;

public:

    SerialIDGenerator(const lot_type& lot = lot_type())
        : next_(lot_add(identifier_type(), lot))
    {
        ;
    }

    identifier_type operator()()
    {
        return serial_advance(next_, 1);
    }

#ifdef WITH_HDF5
    void save(H5::CommonFG* root) const
    {
        using namespace H5;

        boost::scoped_ptr<DataType> optype(new DataType(H5T_OPAQUE, 1));
        hsize_t bufsize(sizeof(identifier_type));
        DataSpace dataspace(1, &bufsize);
        optype->setTag("SerialIDGenerator state type");
        boost::scoped_ptr<DataSet> dataset(
            new DataSet(root->createDataSet("idgen", *optype, dataspace)));
        dataset->write((unsigned char*)(&next_), *optype);
    }

    void load(const H5::CommonFG& root)
    {
        using namespace H5;

        const DataSet dataset(DataSet(root.openDataSet("idgen")));
        boost::scoped_ptr<DataType> optype(new DataType(H5T_OPAQUE, 1));
        optype->setTag("SerialIDGenerator state type");
        identifier_type state;
        dataset.read((unsigned char*)(&state), *optype);
        next_ = state;
    }
#endif

private:

    identifier_type next_;
};

} // ecell4

#endif /* ECELL4_SERIAL_ID_GENERATOR_HPP */
