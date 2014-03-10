#ifndef OBJECTMATRIX_PEER_SET_INDEXING_SUITE_HPP
#define OBJECTMATRIX_PEER_SET_INDEXING_SUITE_HPP

#include <cstddef>
#include <boost/python/def_visitor.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/call_method.hpp>

namespace peer { namespace util {

template<typename Tcontainer_, typename TderivedPolicies_>
struct set_indexing_suite;

namespace detail
{
    template <class Tcontainer_>
    struct final_set_derived_policies 
        : public set_indexing_suite<Tcontainer_,
            final_set_derived_policies<Tcontainer_> > {};
}


template<typename Tcontainer_,
        typename TderivedPolicies_ = detail::final_set_derived_policies<Tcontainer_> >
struct set_indexing_suite: public boost::python::def_visitor<set_indexing_suite<Tcontainer_, TderivedPolicies_> >
{
    typedef typename Tcontainer_::value_type value_type;

    static std::size_t __len__(Tcontainer_& cntnr)
    {
        return TderivedPolicies_::size(cntnr);
    }

    static std::size_t size(Tcontainer_& cntnr)
    {
        return cntnr.size(); 
    }

    static bool add(Tcontainer_& cntnr, value_type const& item)
    {
        return TderivedPolicies_::insert(cntnr, item);
    }

    static bool insert(Tcontainer_& cntnr, value_type const& item)
    {
        return cntnr.insert(item).second;
    }

    static void remove(Tcontainer_& cntnr, value_type const& val)
    {
        if (!TderivedPolicies_::erase(cntnr, val)) {
            PyErr_SetObject(PyExc_KeyError,
                boost::python::incref((boost::python::object(val)).ptr()));
            boost::python::throw_error_already_set();
        }
    }

    static bool erase(Tcontainer_& cntnr, value_type const& val)
    {
        return cntnr.erase(val);
    }

    static bool __contains__(Tcontainer_& cntnr, value_type const& val)
    {
        return TderivedPolicies_::contains(cntnr, val);
    }

    static bool contains(Tcontainer_& cntnr, value_type const& val)
    {
        return cntnr.find(val) != cntnr.end();
    }

    template<typename Tclass_>
    void visit(Tclass_& klass) const
    {
        klass
            .def("__len__", &__len__)
            .def("add", &add)
            .def("remove", &remove)
            .def("__contains__", &__contains__)
            .def("__iter__", boost::python::iterator<Tcontainer_>())
            ;
    }
};

}} // namespace peer

#endif // OBJECTMATRIX_PEER_SET_INDEXING_SUITE_HPP
