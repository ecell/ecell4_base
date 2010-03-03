#define BOOST_TEST_MODULE "World_test"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/scoped_ptr.hpp>
#include "Defs.hpp"
#include "World.hpp"
#include "SerialIDGenerator.hpp"
#include "abstract_set.hpp"

BOOST_AUTO_TEST_CASE(add_species)
{
    typedef World<CyclicWorldTraits<Real, Real> > world_type;
    typedef world_type::species_id_type species_id;
    typedef world_type::species_type species;
    typedef SerialIDGenerator<species_id> id_generator;
    typedef world_type::particle_id_pair_and_distance_list particle_id_pair_and_distance_list;

    world_type i;
    id_generator gen;
    species s1(species(gen(), .3, .05));
    species s2(species(gen(), .2, .03));
    species s3(species(gen(), .1, .02));

    BOOST_CHECK_EQUAL(boost::size(i.get_species()), 0);
    i.add_species(s1);
    BOOST_CHECK_EQUAL(boost::size(i.get_species()), 1);
    BOOST_CHECK(contains(i.get_species(), s1));
    BOOST_CHECK(!contains(i.get_species(), s2));
    BOOST_CHECK(!contains(i.get_species(), s3));
    i.add_species(s2);
    BOOST_CHECK_EQUAL(boost::size(i.get_species()), 2);
    BOOST_CHECK(contains(i.get_species(), s1));
    BOOST_CHECK(contains(i.get_species(), s2));
    BOOST_CHECK(!contains(i.get_species(), s3));
}

BOOST_AUTO_TEST_CASE(new_particles)
{
    typedef World<CyclicWorldTraits<Real, Real> > world_type;
    typedef world_type::species_id_type species_id;
    typedef world_type::species_id_type species_id_type;
    typedef world_type::species_type species;
    typedef world_type::position_type position_type;
    typedef world_type::particle_id_pair particle_id_pair;
    typedef SerialIDGenerator<species_id> id_generator;
    typedef world_type::particle_id_pair_and_distance_list particle_id_pair_and_distance_list;

    world_type i;
    id_generator gen;
    species s1(species(gen(), .3, .05));
    species s2(species(gen(), .3, .08));
    i.add_species(s1);
    i.add_species(s2);

    BOOST_CHECK(contains(i.get_species(), s1));
    BOOST_CHECK(contains(i.get_species(), s2));

    particle_id_pair p1(i.new_particle(s1.id(), position_type(.2, .2, .2)));
    particle_id_pair p2(i.new_particle(s2.id(), position_type(.29, .27, .28)));
    BOOST_CHECK(p2.first != p1.first);

    BOOST_CHECK(!boost::scoped_ptr<particle_id_pair_and_distance_list>(i.check_overlap(p1.second.shape(), array_gen(p1.first))));
    BOOST_CHECK(!boost::scoped_ptr<particle_id_pair_and_distance_list>(i.check_overlap(p2.second.shape(), array_gen(p2.first))));

    BOOST_CHECK(!i.check_overlap(p1));
    BOOST_CHECK(!i.check_overlap(p2));

    particle_id_pair p3(i.new_particle(s1.id(), position_type(.35, .32, .34)));
    BOOST_CHECK(p3.first != p1.first);
    BOOST_CHECK(p3.first != p2.first);

    BOOST_CHECK(!boost::scoped_ptr<particle_id_pair_and_distance_list>(i.check_overlap(p1.second.shape(), array_gen(p1.first))));
    BOOST_CHECK(boost::scoped_ptr<particle_id_pair_and_distance_list>(i.check_overlap(p2.second.shape(), array_gen(p2.first))));
    BOOST_CHECK(boost::scoped_ptr<particle_id_pair_and_distance_list>(i.check_overlap(p2.second.shape(), array_gen(p1.first, p2.first))));
    BOOST_CHECK(!boost::scoped_ptr<particle_id_pair_and_distance_list>(i.check_overlap(p2.second.shape(), array_gen(p2.first, p3.first))));
    BOOST_CHECK(boost::scoped_ptr<particle_id_pair_and_distance_list>(i.check_overlap(p3.second.shape(), array_gen(p3.first))));
    BOOST_CHECK(!boost::scoped_ptr<particle_id_pair_and_distance_list>(i.check_overlap(p3.second.shape(), array_gen(p2.first, p3.first))));
    BOOST_CHECK(boost::scoped_ptr<particle_id_pair_and_distance_list>(i.check_overlap(p3.second.shape(), array_gen(p1.first, p3.first))));

    BOOST_CHECK(!i.check_overlap(p1));
    BOOST_CHECK(i.check_overlap(p2));
    BOOST_CHECK(boost::scoped_ptr<particle_id_pair_and_distance_list>(i.check_overlap(p2, array_gen(p1.first))));
    BOOST_CHECK(!boost::scoped_ptr<particle_id_pair_and_distance_list>(i.check_overlap(p2, array_gen(p3.first))));
    BOOST_CHECK(i.check_overlap(p3));
    BOOST_CHECK(!boost::scoped_ptr<particle_id_pair_and_distance_list>(i.check_overlap(p3, array_gen(p2.first))));
    BOOST_CHECK(boost::scoped_ptr<particle_id_pair_and_distance_list>(i.check_overlap(p3, array_gen(p1.first))));
}

BOOST_AUTO_TEST_CASE(get_particle)
{
    typedef World<CyclicWorldTraits<Real, Real> > world_type;
    typedef world_type::species_id_type species_id;
    typedef world_type::species_id_type species_id_type;
    typedef world_type::species_type species;
    typedef world_type::position_type position_type;
    typedef world_type::particle_id_pair particle_id_pair;
    typedef SerialIDGenerator<species_id> id_generator;
    typedef world_type::particle_id_pair_and_distance_list particle_id_pair_and_distance_list;

    world_type i;
    id_generator gen;
    species s1(species(gen(), .3, .05));
    species s2(species(gen(), .3, .08));
    i.add_species(s1);
    i.add_species(s2);

    particle_id_pair p1(i.new_particle(s1.id(), position_type(.2, .2, .2)));
    particle_id_pair p2(i.new_particle(s2.id(), position_type(.29, .27, .28)));

    BOOST_CHECK(i.get_particle(p1.first).second == p1.second);
    BOOST_CHECK(i.get_particle(p2.first).second == p2.second);
}

BOOST_AUTO_TEST_CASE(transaction_1)
{
    typedef World<CyclicWorldTraits<Real, Real> > world_type;
    typedef world_type::species_id_type species_id;
    typedef world_type::species_id_type species_id_type;
    typedef world_type::species_type species;
    typedef world_type::position_type position_type;
    typedef world_type::particle_id_pair particle_id_pair;
    typedef SerialIDGenerator<species_id> id_generator;
    typedef world_type::particle_id_pair_and_distance_list particle_id_pair_and_distance_list;

    world_type i;
    id_generator gen;
    species s1(species(gen(), .3, .05));
    species s2(species(gen(), .3, .08));
    i.add_species(s1);
    i.add_species(s2);

    particle_id_pair p1, p2;
    {
        boost::scoped_ptr<world_type::transaction_type> tx(i.create_transaction());

        new (&p1) particle_id_pair(tx->new_particle(s1.id(), position_type(.2, .2, .2)));
        new (&p2) particle_id_pair(tx->new_particle(s2.id(), position_type(.29, .27, .28)));

        BOOST_CHECK(tx->get_particle(p1.first).second == p1.second);
        BOOST_CHECK(tx->get_particle(p2.first).second == p2.second);
        BOOST_CHECK(i.get_particle(p1.first).second == p1.second);
        BOOST_CHECK(i.get_particle(p2.first).second == p2.second);

        BOOST_CHECK(cue(*boost::scoped_ptr<world_type::particle_id_pair_generator>(tx->get_added_particles()), p1));
        BOOST_CHECK(cue(*boost::scoped_ptr<world_type::particle_id_pair_generator>(tx->get_added_particles()), p2));
        tx->rollback();
    }

    BOOST_CHECK_THROW(i.get_particle(p1.first), not_found);
    BOOST_CHECK_THROW(i.get_particle(p2.first), not_found);
}

BOOST_AUTO_TEST_CASE(transaction_2)
{
    typedef World<CyclicWorldTraits<Real, Real> > world_type;
    typedef world_type::species_id_type species_id;
    typedef world_type::species_id_type species_id_type;
    typedef world_type::species_type species;
    typedef world_type::position_type position_type;
    typedef world_type::particle_id_pair particle_id_pair;
    typedef SerialIDGenerator<species_id> id_generator;
    typedef world_type::particle_id_pair_and_distance_list particle_id_pair_and_distance_list;

    world_type i;
    id_generator gen;
    species s1(species(gen(), .3, .05));
    species s2(species(gen(), .3, .08));
    i.add_species(s1);
    i.add_species(s2);

    particle_id_pair p1, p2;
    {
        boost::scoped_ptr<world_type::transaction_type> tx(i.create_transaction());

        new (&p1) particle_id_pair(tx->new_particle(s1.id(), position_type(.2, .2, .2)));
        new (&p2) particle_id_pair(tx->new_particle(s2.id(), position_type(.29, .27, .28)));

        BOOST_CHECK(tx->get_particle(p1.first).second == p1.second);
        BOOST_CHECK(tx->get_particle(p2.first).second == p2.second);
        BOOST_CHECK(i.get_particle(p1.first).second == p1.second);
        BOOST_CHECK(i.get_particle(p2.first).second == p2.second);
        BOOST_CHECK(cue(*boost::scoped_ptr<world_type::particle_id_pair_generator>(tx->get_added_particles()), p1));
        BOOST_CHECK(cue(*boost::scoped_ptr<world_type::particle_id_pair_generator>(tx->get_added_particles()), p2));
    }

    BOOST_CHECK(i.get_particle(p1.first).second == p1.second);
    BOOST_CHECK(i.get_particle(p2.first).second == p2.second);
}

BOOST_AUTO_TEST_CASE(transaction_3)
{
    typedef World<CyclicWorldTraits<Real, Real> > world_type;
    typedef world_type::species_id_type species_id;
    typedef world_type::species_id_type species_id_type;
    typedef world_type::species_type species;
    typedef world_type::position_type position_type;
    typedef world_type::particle_id_pair particle_id_pair;
    typedef SerialIDGenerator<species_id> id_generator;
    typedef world_type::particle_id_pair_and_distance_list particle_id_pair_and_distance_list;

    world_type i;
    id_generator gen;
    species s1(species(gen(), .3, .05));
    species s2(species(gen(), .3, .08));
    i.add_species(s1);
    i.add_species(s2);

    particle_id_pair p1(i.new_particle(s1.id(), position_type(.2, .2, .2)));
    particle_id_pair p2(i.new_particle(s2.id(), position_type(.29, .27, .28)));
    particle_id_pair p3;

    {
        boost::scoped_ptr<world_type::transaction_type> tx(i.create_transaction());

        new (&p3) particle_id_pair(tx->new_particle(s1.id(), position_type(.4, .2, .1)));

        BOOST_CHECK(tx->get_particle(p1.first).second == p1.second);
        BOOST_CHECK(tx->get_particle(p2.first).second == p2.second);
        BOOST_CHECK(tx->get_particle(p3.first).second == p3.second);
        BOOST_CHECK(i.get_particle(p1.first).second == p1.second);
        BOOST_CHECK(i.get_particle(p2.first).second == p2.second);
        BOOST_CHECK(i.get_particle(p3.first).second == p3.second);

        tx->remove_particle(p1.first);
        BOOST_CHECK_THROW(tx->get_particle(p1.first), not_found);
        BOOST_CHECK(tx->get_particle(p2.first).second == p2.second);
        BOOST_CHECK_THROW(i.get_particle(p1.first), not_found);

        BOOST_CHECK(!cue(*boost::scoped_ptr<world_type::particle_id_pair_generator>(tx->get_added_particles()), p1));
        BOOST_CHECK(!cue(*boost::scoped_ptr<world_type::particle_id_pair_generator>(tx->get_added_particles()), p2));
        BOOST_CHECK(cue(*boost::scoped_ptr<world_type::particle_id_pair_generator>(tx->get_added_particles()), p3));
        BOOST_CHECK(cue(*boost::scoped_ptr<world_type::particle_id_pair_generator>(tx->get_removed_particles()), p1));
        BOOST_CHECK(!cue(*boost::scoped_ptr<world_type::particle_id_pair_generator>(tx->get_removed_particles()), p2));
        BOOST_CHECK(!cue(*boost::scoped_ptr<world_type::particle_id_pair_generator>(tx->get_removed_particles()), p3));
        BOOST_CHECK(!cue(*boost::scoped_ptr<world_type::particle_id_pair_generator>(tx->get_modified_particles()), p1));
        BOOST_CHECK(!cue(*boost::scoped_ptr<world_type::particle_id_pair_generator>(tx->get_modified_particles()), p2));
        BOOST_CHECK(!cue(*boost::scoped_ptr<world_type::particle_id_pair_generator>(tx->get_modified_particles()), p3));

        tx->rollback();
    }

    BOOST_CHECK(i.get_particle(p1.first).second == p1.second);
    BOOST_CHECK(i.get_particle(p2.first).second == p2.second);
    BOOST_CHECK_THROW(i.get_particle(p3.first), not_found);
}


