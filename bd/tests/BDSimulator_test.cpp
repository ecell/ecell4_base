#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include "../BDSimulator.hpp"

using namespace ecell4;
using namespace ecell4::bd;


class BDSimulatorTest
    : public CppUnit::TestFixture
{
public:

    CPPUNIT_TEST_SUITE(BDSimulatorTest);
    CPPUNIT_TEST(test_step);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp()
    {
        Real const L(1e-6);
        Position3 const edge_lengths(L, L, L);

        boost::shared_ptr<Model> model(new NetworkModel());
        boost::shared_ptr<BDWorld> world(new BDWorld(edge_lengths));
        GSLRandomNumberGenerator rng;

        target = new BDSimulator(model, world, rng);
    }

    void tearDown()
    {
        delete target;
    }

    void test_step();

private:

    BDSimulator *target;
};

CPPUNIT_TEST_SUITE_REGISTRATION(BDSimulatorTest);

void BDSimulatorTest::test_step()
{
    target->step();
}
