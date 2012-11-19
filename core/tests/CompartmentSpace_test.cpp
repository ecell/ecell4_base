#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

#include "../types.hpp"
#include "../CompartmentSpace.hpp"

using namespace ecell4;


template <typename Timpl_>
class CompartmentSpaceTest
    : public CppUnit::TestFixture
{
public:

    typedef Timpl_ implementation_type;

    CPPUNIT_TEST_SUITE(CompartmentSpaceTest);
    CPPUNIT_TEST(test_volume);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp()
    {
        Real const volume(1e-18);
        target = new implementation_type(volume);
    }

    void tearDown()
    {
        delete target;
    }

    void test_volume();

private:

    CompartmentSpace *target;
};

CPPUNIT_TEST_SUITE_REGISTRATION(
    CompartmentSpaceTest<CompartmentSpaceVectorImpl>);

template <typename Timpl_>
void CompartmentSpaceTest<Timpl_>::test_volume()
{
    target->set_volume(2 * target->volume());
}
