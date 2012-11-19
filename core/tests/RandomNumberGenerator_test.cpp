#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

#include "../RandomNumberGenerator.hpp"

using namespace ecell4;


template <typename Timpl_>
class RandomNumberGeneratorTest
    : public CppUnit::TestFixture
{
public:

    typedef Timpl_ implementation_type;

    CPPUNIT_TEST_SUITE(RandomNumberGeneratorTest);
    CPPUNIT_TEST(test_seed);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp()
    {
        target = new implementation_type();
    }

    void tearDown()
    {
        delete target;
    }

    void test_seed();

private:

    RandomNumberGenerator *target;
};

CPPUNIT_TEST_SUITE_REGISTRATION(
    RandomNumberGeneratorTest<GSLRandomNumberGenerator>);

template <typename Timpl_>
void RandomNumberGeneratorTest<Timpl_>::test_seed()
{
    target->seed(0);
    // CPPUNIT_ASSERT_EQUAL(1, 2);
}
