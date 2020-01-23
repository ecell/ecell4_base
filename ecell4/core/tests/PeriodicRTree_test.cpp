#define BOOST_TEST_MODULE "PeriodicRTree_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif


#include <boost/test/floating_point_comparison.hpp>

#include <ecell4/core/Particle.hpp>
#include <ecell4/core/PeriodicRTree.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/Identifier.hpp>
#include <tuple>
#include <random>

using namespace ecell4;

struct AlwaysTightAABBGetter
{
    AABB operator()(const Particle& p, const Real) const noexcept
    {
        const Real3 radius(p.radius(), p.radius(), p.radius());
        return AABB(p.position() - radius, p.position() + radius);
    }
};

struct FixedMarginAABBGetter
{
    AABB operator()(const Particle& p, const Real margin) const noexcept
    {
        const Real3 radius(p.radius() + margin,
                           p.radius() + margin,
                           p.radius() + margin);
        return AABB(p.position() - radius, p.position() + radius);
    }
};

struct ScaledMarginAABBGetter
{
    AABB operator()(const Particle& p, const Real margin) const noexcept
    {
        const Real3 radius(p.radius() + p.D() * margin,
                           p.radius() + p.D() * margin,
                           p.radius() + p.D() * margin);
        return AABB(p.position() - radius, p.position() + radius);
    }
};

using aabb_getters = std::tuple<
    AlwaysTightAABBGetter, FixedMarginAABBGetter, ScaledMarginAABBGetter>;

struct Query
{
    ParticleID ignore;
    Real3      center;
    Real       radius;

    bool operator()(const std::pair<ParticleID, Particle>& pidp,
                    const Real3& edges) const noexcept
    {
        if(pidp.first == ignore) {return false;}
        const auto rr = radius + pidp.second.radius();
        const auto rhs = this->periodic_transpose(pidp.second.position(),
                                                  center, edges);
        return length_sq(rhs - center) < rr * rr;
    }

    bool operator()(const AABB& box, const Real3& edge_lengths) const noexcept
    {
        return this->distance_sq(box, center, edge_lengths) < radius * radius;
    }

  private:

    // AABB-sphere intersection query under the PBC
    Real distance_sq(const AABB& box, Real3 pos, const Real3& edge_lengths) const noexcept
    {
        pos = periodic_transpose(pos, (box.upper() + box.lower()) * 0.5, edge_lengths);

        Real dist_sq = 0;
        for(std::size_t i=0; i<3; ++i)
        {
            const auto v = pos[i];
            if(v < box.lower()[i])
            {
                dist_sq += (v - box.lower()[i]) * (v - box.lower()[i]);
            }
            else if(box.upper()[i] < v)
            {
                dist_sq += (v - box.upper()[i]) * (v - box.upper()[i]);
            }
        }
        return dist_sq;
    }

    // transpose a position based on the periodic boundary condition.
    Real3 periodic_transpose(
        const Real3& pos1, const Real3& pos2, const Real3& edges) const
    {
        Real3 retval(pos1);
        for(Real3::size_type dim(0); dim < 3; ++dim)
        {
            const Real edge_length(edges[dim]);
            const Real diff(pos2[dim] - pos1[dim]), half(edge_length * 0.5);

            if (half < diff)
            {
                retval[dim] += edge_length;
            }
            else if (diff < -half)
            {
                retval[dim] -= edge_length;
            }
        }
        return retval;
    }
};

BOOST_AUTO_TEST_CASE_TEMPLATE(PeriodicRTree_query, AABBGetter, aabb_getters)
{
    constexpr std::size_t N = 500;
    constexpr Real        L = 1.0;
    const Real3 edge_lengths(L, 2*L, 3*L);
    std::mt19937 mt(123456789);
    std::uniform_real_distribution<Real> uni(0.0, L);

    const Species sp("A");
    const Real radius = 0.005;
    const Real D      = 1.0;

    PeriodicRTree<ParticleID, Particle, AABBGetter> tree(edge_lengths, 0.01);
    BOOST_TEST_MESSAGE("tree constructed");

    std::vector<std::pair<ParticleID, Particle>> full_list;
    SerialIDGenerator<ParticleID> pidgen;
    for(std::size_t i=0; i<N; ++i)
    {
        const Real3 pos(uni(mt), 2 * uni(mt), 3 * uni(mt));
        const auto  pid = pidgen();
        const auto  p   = Particle(sp, pos, radius, D);
        full_list.emplace_back(pid, p);

        tree.insert(pid, p);

        BOOST_REQUIRE(tree.diagnosis());
    }
    BOOST_TEST_MESSAGE("objects are inserted");

    // ----------------------------------------------------------------------
    // send a query and check particle can be found

    std::vector<std::pair<ParticleID, Particle>> query_results;

    const ParticleID nil = pidgen();
    for(const auto& pidp : full_list)
    {
        BOOST_REQUIRE(tree.has(pidp.first));
        BOOST_REQUIRE(tree.get(pidp.first) == pidp);

        const Query q{nil, pidp.second.position(), pidp.second.radius()};

        tree.query(q, std::back_inserter(query_results));
        BOOST_TEST((std::find_if(query_results.begin(), query_results.end(),
                    [&pidp](const std::pair<ParticleID, Particle>& lhs) -> bool {
                        return lhs.first == pidp.first;
                    }) != query_results.end()));
        query_results.clear();
    }

    // ----------------------------------------------------------------------
    // send a query and check all the possible collisions are detected

    for(std::size_t i=0; i<N; ++i)
    {
        const Real3 query_center(uni(mt), 2 * uni(mt), 3 * uni(mt));
        const Real  query_range = uni(mt) * L * 0.1;
        const Query query{full_list.front().first, query_center, query_range};

        tree.query(query, std::back_inserter(query_results));

        bool all_found = true;
        for(const auto& pidp : full_list)
        {
            if(query(pidp, edge_lengths))
            {
                const auto found = std::find_if(
                    query_results.begin(), query_results.end(),
                    [&pidp](const std::pair<ParticleID, Particle>& lhs) -> bool {
                        return lhs.first == pidp.first;
                    });
                if(found == query_results.end())
                {
                    all_found = false;
                }
            }
        }
        BOOST_TEST(all_found);
        query_results.clear();
    }
    BOOST_TEST_MESSAGE("query is tested");

    // ----------------------------------------------------------------------
    // check query results after erase/insert

    for(std::size_t i=0; i<N; ++i)
    {
        BOOST_TEST_MESSAGE("erase/insert particle " << i << ", "
                << full_list.at(i).first << ":" << full_list.at(i).second);

        const Real3 pos(uni(mt), 2 * uni(mt), 3 * uni(mt));
        const Particle p(sp, pos, radius, D);

        const auto old = full_list.at(i);
        tree.erase(old);

        BOOST_REQUIRE(tree.diagnosis());
        BOOST_TEST_MESSAGE("tree is still okay");

        {
            // make sure that the tree no longer contains full_list.at(i).
            const Query q{nil, old.second.position(), old.second.radius()};
            tree.query(q, std::back_inserter(query_results));
            BOOST_TEST((std::find_if(query_results.begin(), query_results.end(),
                        [&old](const std::pair<ParticleID, Particle>& lhs) -> bool {
                            return lhs.first == old.first;
                        }) == query_results.end()));
            query_results.clear();
            BOOST_TEST_MESSAGE("old value successfully removed");
        }

        BOOST_REQUIRE(tree.diagnosis());
        BOOST_TEST_MESSAGE("tree is still okay");

        full_list.at(i).second = p;
        tree.insert(full_list.at(i));

        const auto novel = full_list.at(i);

        {
            const Query q{nil, novel.second.position(), novel.second.radius()};
            tree.query(q, std::back_inserter(query_results));
            BOOST_TEST((std::find_if(query_results.begin(), query_results.end(),
                        [&novel](const std::pair<ParticleID, Particle>& lhs) -> bool {
                            return lhs.first == novel.first;
                        }) != query_results.end()));
            query_results.clear();
            BOOST_TEST_MESSAGE("new value successfully inserted");
        }

        BOOST_REQUIRE(tree.diagnosis());
        BOOST_TEST_MESSAGE("tree is still okay");
    }
    BOOST_TEST_MESSAGE("objects are updated");

    for(std::size_t i=0; i<N; ++i)
    {
        const Real3 query_center(uni(mt), 2 * uni(mt), 3 * uni(mt));
        const Real  query_range = uni(mt) * L * 0.1;
        const Query query{full_list.front().first, query_center, query_range};

        tree.query(query, std::back_inserter(query_results));

        for(const auto& pidp : full_list)
        {
            if(query(pidp, edge_lengths))
            {
                const auto found = std::find_if(
                    query_results.begin(), query_results.end(),
                    [&pidp](const std::pair<ParticleID, Particle>& lhs) -> bool {
                        return lhs.first == pidp.first;
                    });

                BOOST_CHECK(found != query_results.end());
            }
        }
        query_results.clear();
    }

    // ----------------------------------------------------------------------
    // update()

    for(std::size_t i=0; i<N; ++i)
    {
        BOOST_TEST_MESSAGE("updating particle " << i << ", "
                << full_list.at(i).first << ":" << full_list.at(i).second);

        const auto id = full_list.at(i).first;
        const Real3 pos(uni(mt), 2 * uni(mt), 3 * uni(mt));
        const Particle p(sp, pos, radius, D);

        full_list.at(i).second = p;
        tree.update(id, p);

        BOOST_REQUIRE(tree.diagnosis());
        BOOST_TEST_MESSAGE("tree is still okay");

        const auto novel = full_list.at(i);
        {
            const Query q{nil, novel.second.position(), novel.second.radius()};
            tree.query(q, std::back_inserter(query_results));
            BOOST_TEST((std::find_if(query_results.begin(), query_results.end(),
                        [&novel](const std::pair<ParticleID, Particle>& lhs) -> bool {
                            return lhs.first == novel.first;
                        }) != query_results.end()));
            query_results.clear();
            BOOST_TEST_MESSAGE("new value successfully inserted");
        }
    }
    BOOST_TEST_MESSAGE("objects are updated");

    for(std::size_t i=0; i<N; ++i)
    {
        const Real3 query_center(uni(mt), 2 * uni(mt), 3 * uni(mt));
        const Real  query_range = uni(mt) * L * 0.1;
        const Query query{full_list.front().first, query_center, query_range};

        tree.query(query, std::back_inserter(query_results));

        for(const auto& pidp : full_list)
        {
            if(query(pidp, edge_lengths))
            {
                const auto found = std::find_if(
                    query_results.begin(), query_results.end(),
                    [&pidp](const std::pair<ParticleID, Particle>& lhs) -> bool {
                        return lhs.first == pidp.first;
                    });

                BOOST_CHECK(found != query_results.end());
            }
        }
        query_results.clear();
    }
}
