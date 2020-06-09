#define BOOST_TEST_MODULE "ObjectIDContainer_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/tools/floating_point_comparison.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/ObjectIDContainer.hpp>
#include <ecell4/core/Identifier.hpp>
#include <ecell4/core/Particle.hpp>

constexpr std::size_t N = 1000;

BOOST_AUTO_TEST_CASE(ObjectIDContainer_add_update_remove)
{
    using namespace ecell4;
    ObjectIDContainer<ParticleID, Particle> container;
    SerialIDGenerator<ParticleID> idgen;

    std::vector<ParticleID> ids;
    for(std::size_t i=0; i<N; ++i)
    {
        ids.push_back(idgen());
        const bool newly_added = container.update(ids.back(), Particle());
        BOOST_CHECK(newly_added);
        BOOST_CHECK(container.diagnosis());
    }
    BOOST_CHECK_EQUAL(container.size(), N);

    for(const auto& id : ids)
    {
        BOOST_CHECK(container.has(id));
    }

    for(const auto& id : ids)
    {
        const bool newly_added = container.update(id, Particle());
        BOOST_CHECK(!newly_added);
        BOOST_CHECK(container.diagnosis());
    }

    for(const auto& id : ids)
    {
        container.remove(id);
        BOOST_CHECK(container.diagnosis());
    }
    BOOST_CHECK(container.diagnosis());
    BOOST_CHECK_EQUAL(container.size(), 0);
    BOOST_CHECK(container.empty());
}

