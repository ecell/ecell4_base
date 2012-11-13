#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

#include "../RandomNumberGenerator.hpp"


class GSLRandomNumberGeneratorTest
    : public CppUnit::TestFixture
{
public:

    CPPUNIT_TEST_SUITE(GSLRandomNumberGeneratorTest);
    CPPUNIT_TEST (test_seed);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp ()
    {
        target = new GSLRandomNumberGenerator();
    }

    void tearDown ()
    {
        delete target;
    }

    void test_seed();

private:

    GSLRandomNumberGenerator *target;
};

CPPUNIT_TEST_SUITE_REGISTRATION(GSLRandomNumberGeneratorTest);

void GSLRandomNumberGeneratorTest::test_seed()
{
    CPPUNIT_ASSERT_EQUAL(1, 2);
}
