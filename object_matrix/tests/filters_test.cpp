#define BOOST_TEST_MODULE "filters_test"

#include <boost/test/included/unit_test.hpp>
#include "object_container.hpp"
#include "sphere.hpp"
#include "filters.hpp"

template<typename Toc_>
struct collector
{
    void operator()(typename Toc_::iterator i,
            typename Toc_::position_type::value_type dist)
    {
        std::cout << (*i).second << ", " << dist << std::endl;
    }
};


BOOST_AUTO_TEST_CASE(basic)
{
    typedef object_container<sphere<double>, int> oc_type;
    typedef oc_type::position_type pos;
    oc_type oc(1.0, 10);

    oc.update(std::make_pair(0, oc_type::mapped_type(pos(0.2, 0.6, 0.4), 0.15)));
    oc.update(std::make_pair(1, oc_type::mapped_type(pos(0.2, 0.7, 0.5), 0.05)));
    oc.update(std::make_pair(2, oc_type::mapped_type(pos(0.9, 0.1, 0.4), 0.07)));
    oc.update(std::make_pair(3, oc_type::mapped_type(pos(0.9, 0.95, 0.4), 0.1)));

    collector<oc_type> col;
    oc_type::const_iterator f(oc.find(1));
    take_neighbor(oc, col, (*f).second);
}
