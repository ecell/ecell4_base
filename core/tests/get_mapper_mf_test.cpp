#define BOOST_TEST_MODULE "get_mapper_mf_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/type_traits.hpp>

#include "../get_mapper_mf.hpp"
#include <map>

using namespace ecell4;

BOOST_AUTO_TEST_CASE(get_mapper_mf_test_type)
{
  typedef utils::get_mapper_mf<std::string, std::string>::type string_map_type;
  BOOST_CHECK((boost::is_same<string_map_type, std::map<std::string, std::string> >::value));
}
