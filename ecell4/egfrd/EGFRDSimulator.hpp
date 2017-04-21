#ifndef EGFRDSIMULATOR_HPP
#define EGFRDSIMULATOR_HPP

#include <boost/bind.hpp>
#include <boost/array.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/fusion/container/map.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/sequence/intrinsic/at_key.hpp>
#include <boost/fusion/sequence/intrinsic/value_at_key.hpp>
#include <boost/fusion/include/at_key.hpp>
#include <boost/none_t.hpp>
#include <boost/variant.hpp>

#include <gsl/gsl_sf_log.h>

#include <ecell4/core/get_mapper_mf.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/EventScheduler.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>

#include "utils/array_helper.hpp"
//#include "utils/get_mapper_mf.hpp"
#include "utils/fun_composition.hpp"
#include "utils/fun_wrappers.hpp"
#include "utils/pointer_as_ref.hpp"
#include "utils/pair.hpp"
#include "utils/math.hpp"
#include "utils/stringizer.hpp"
#include "ShellID.hpp"
#include "DomainID.hpp"
#include "Shell.hpp"
//#include "EventScheduler.hpp"
#include "ParticleSimulator.hpp"
#include "MatrixSpace.hpp"
#include "AnalyticalSingle.hpp"
#include "AnalyticalPair.hpp"
#include "Multi.hpp"

#include <greens_functions/PairGreensFunction.hpp>
#include <greens_functions/GreensFunction3DRadAbs.hpp>
#include <greens_functions/GreensFunction3DRadInf.hpp>
#include <greens_functions/GreensFunction3DAbsSym.hpp>
#include <greens_functions/GreensFunction3DAbs.hpp>
#include <greens_functions/GreensFunction3D.hpp>
// using namespace greens_functions;


template<typename Tworld_>
struct EGFRDSimulatorTraitsBase: public ParticleSimulatorTraitsBase<Tworld_>
{
    typedef ParticleSimulatorTraitsBase<Tworld_> base_type;
    typedef Tworld_ world_type;

    typedef ShellID shell_id_type;
    typedef DomainID domain_id_type;
    typedef ecell4::SerialIDGenerator<shell_id_type> shell_id_generator;
    typedef ecell4::SerialIDGenerator<domain_id_type> domain_id_generator;
    typedef Domain<EGFRDSimulatorTraitsBase> domain_type;
    typedef std::pair<const domain_id_type, boost::shared_ptr<domain_type> > domain_id_pair;
    typedef ecell4::EventScheduler event_scheduler_type; // base_type::time_type == ecell4::Real
    // typedef EventScheduler<typename base_type::time_type> event_scheduler_type;

    typedef typename event_scheduler_type::identifier_type event_id_type;
    typedef typename event_scheduler_type::value_type event_id_pair_type;
    typedef ecell4::Event event_type;

    template<typename Tshape_>
    struct shell_generator
    {
        typedef Shell<Tshape_, domain_id_type> type;
    };

    static const Real safety();
    static const Real single_shell_factor();
    static const Real default_dt_factor();
    static const Real cutoff_factor();

    static const Real SAFETY;
    static const Real SINGLE_SHELL_FACTOR;
    static const Real DEFAULT_DT_FACTOR;
    static const Real CUTOFF_FACTOR;
};

template<typename Tworld_>
const Real EGFRDSimulatorTraitsBase<Tworld_>::safety() { return 1. + 1e-5; }
template<typename Tworld_>
const Real EGFRDSimulatorTraitsBase<Tworld_>::single_shell_factor() { return 0.1; }
template<typename Tworld_>
const Real EGFRDSimulatorTraitsBase<Tworld_>::default_dt_factor() { return 1e-5; }
template<typename Tworld_>
const Real EGFRDSimulatorTraitsBase<Tworld_>::cutoff_factor() { return 5.6; }

template<typename Tworld_>
const Real EGFRDSimulatorTraitsBase<Tworld_>::SAFETY = EGFRDSimulatorTraitsBase<Tworld_>::safety();
template<typename Tworld_>
const Real EGFRDSimulatorTraitsBase<Tworld_>::SINGLE_SHELL_FACTOR = EGFRDSimulatorTraitsBase<Tworld_>::single_shell_factor();
template<typename Tworld_>
const Real EGFRDSimulatorTraitsBase<Tworld_>::DEFAULT_DT_FACTOR = EGFRDSimulatorTraitsBase<Tworld_>::default_dt_factor();
template<typename Tworld_>
const Real EGFRDSimulatorTraitsBase<Tworld_>::CUTOFF_FACTOR = EGFRDSimulatorTraitsBase<Tworld_>::cutoff_factor();

namespace detail {

template<typename T_>
struct get_greens_function {};

template<>
struct get_greens_function<ecell4::Sphere>
{
    typedef greens_functions::GreensFunction3DAbsSym type;
};

// template<>
// struct get_greens_function<Cylinder>
// {
//     typedef greens_functions::GreensFunction3DAbsSym type;
// };
// 
// template<>
// struct get_greens_function<Sphere>
// {
//     typedef greens_functions::GreensFunction3DAbsSym type;
// };

template<>
struct get_greens_function<ecell4::Cylinder>
{
    typedef greens_functions::GreensFunction3DAbsSym type;
};

template<typename T_>
struct get_pair_greens_function {};

template<>
struct get_pair_greens_function<ecell4::Sphere>
{
    typedef greens_functions::GreensFunction3DRadAbs iv_type;
    typedef greens_functions::GreensFunction3DAbsSym com_type;
};

template<>
struct get_pair_greens_function<ecell4::Cylinder>
{
    typedef greens_functions::GreensFunction3DRadAbs iv_type;
    typedef greens_functions::GreensFunction3DAbsSym com_type;
};

// template<>
// struct get_pair_greens_function<Sphere>
// {
//     typedef greens_functions::GreensFunction3DRadAbs iv_type;
//     typedef greens_functions::GreensFunction3DAbsSym com_type;
// };
// 
// template<>
// struct get_pair_greens_function<Cylinder>
// {
//     typedef greens_functions::GreensFunction3DRadAbs iv_type;
//     typedef greens_functions::GreensFunction3DAbsSym com_type;
// };

} // namespace detail

template<typename Ttraits_>
class EGFRDSimulator;

template<typename Ttraits_>
struct ImmutativeDomainVisitor
{
    typedef typename EGFRDSimulator<Ttraits_>::multi_type multi_type;
    typedef typename EGFRDSimulator<Ttraits_>::spherical_single_type spherical_single_type;
    typedef typename EGFRDSimulator<Ttraits_>::cylindrical_single_type cylindrical_single_type;
    typedef typename EGFRDSimulator<Ttraits_>::spherical_pair_type spherical_pair_type;
    typedef typename EGFRDSimulator<Ttraits_>::cylindrical_pair_type cylindrical_pair_type;

    virtual ~ImmutativeDomainVisitor() {}

    virtual void operator()(multi_type const&) const = 0;

    virtual void operator()(spherical_single_type const&) const = 0;

    virtual void operator()(cylindrical_single_type const&) const = 0;

    virtual void operator()(spherical_pair_type const&) const = 0;

    virtual void operator()(cylindrical_pair_type const&) const = 0;
};

template<typename Ttraits_>
struct MutativeDomainVisitor
{
    typedef typename EGFRDSimulator<Ttraits_>::multi_type multi_type;
    typedef typename EGFRDSimulator<Ttraits_>::spherical_single_type spherical_single_type;
    typedef typename EGFRDSimulator<Ttraits_>::cylindrical_single_type cylindrical_single_type;
    typedef typename EGFRDSimulator<Ttraits_>::spherical_pair_type spherical_pair_type;
    typedef typename EGFRDSimulator<Ttraits_>::cylindrical_pair_type cylindrical_pair_type;

    virtual ~MutativeDomainVisitor() {}

    virtual void operator()(multi_type&) const = 0;

    virtual void operator()(spherical_single_type&) const = 0;

    virtual void operator()(cylindrical_single_type&) const = 0;

    virtual void operator()(spherical_pair_type&) const = 0;

    virtual void operator()(cylindrical_pair_type&) const = 0;
};


#define CHECK(expr) \
    do \
    { \
        if (!(expr)) { retval = false; LOG_DEBUG(("checking [%s] failed", #expr)); } \
    } while (0)

template<typename Ttraits_>
class EGFRDSimulator: public ParticleSimulator<Ttraits_>
{
public:
    typedef Ttraits_ traits_type;
    typedef ParticleSimulator<Ttraits_> base_type;

    // typedef typename base_type::sphere_type sphere_type;
    // typedef typename base_type::cylinder_type cylinder_type;
    typedef typename base_type::model_type model_type;

    typedef typename traits_type::world_type world_type;
    typedef typename traits_type::domain_id_type domain_id_type;
    typedef typename traits_type::shell_id_type shell_id_type;
    typedef typename traits_type::template shell_generator<ecell4::Sphere>::type spherical_shell_type;
    typedef typename traits_type::template shell_generator<ecell4::Cylinder>::type cylindrical_shell_type;
    // typedef typename traits_type::template shell_generator<sphere_type>::type spherical_shell_type;
    // typedef typename traits_type::template shell_generator<cylinder_type>::type cylindrical_shell_type;
    typedef typename traits_type::domain_type domain_type;
    typedef typename traits_type::domain_id_pair domain_id_pair;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::shell_id_generator shell_id_generator;
    typedef typename traits_type::domain_id_generator domain_id_generator;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename traits_type::reaction_record_type reaction_record_type;
    typedef typename traits_type::reaction_recorder_type reaction_recorder_type;
    typedef typename traits_type::event_scheduler_type event_scheduler_type;
    typedef typename traits_type::event_type event_type;
    typedef typename traits_type::event_id_type event_id_type;
    typedef typename traits_type::event_id_pair_type event_id_pair_type;

    typedef typename world_type::traits_type::length_type length_type;
    typedef typename world_type::traits_type::position_type position_type;
    typedef typename world_type::traits_type::rng_type rng_type;
    typedef typename world_type::traits_type::particle_type particle_type;
    typedef typename world_type::traits_type::D_type D_type;
    typedef typename world_type::traits_type::molecule_info_type molecule_info_type;
    typedef typename world_type::traits_type::species_id_type species_id_type;
    typedef typename world_type::traits_type::structure_type structure_type;
    typedef typename world_type::particle_shape_type particle_shape_type;
    typedef typename world_type::traits_type::particle_id_type particle_id_type;
    typedef typename world_type::particle_id_pair particle_id_pair;
    typedef typename world_type::particle_id_pair_and_distance particle_id_pair_and_distance;
    typedef typename world_type::particle_id_pair_and_distance_list particle_id_pair_and_distance_list;

    typedef std::pair<const shell_id_type, spherical_shell_type> spherical_shell_id_pair;
    typedef std::pair<const shell_id_type, cylindrical_shell_type> cylindrical_shell_id_pair;


    typedef typename world_type::traits_type::particle_simulation_structure_type
        particle_simulation_structure_type;
    // typedef typename world_type::traits_type::spherical_surface_type
    //     spherical_surface_type;
    // typedef typename world_type::traits_type::cylindrical_surface_type
    //     cylindrical_surface_type;
    // typedef typename world_type::traits_type::planar_surface_type planar_surface_type;
    typedef typename world_type::traits_type::cuboidal_region_type cuboidal_region_type;

    typedef typename ReactionRecorderWrapper<reaction_record_type>::reaction_info_type reaction_info_type;

    typedef Single<traits_type> single_type;
    typedef Pair<traits_type> pair_type;
    typedef Multi<EGFRDSimulator> multi_type;
    typedef ShapedDomain<traits_type> shaped_domain_type;
    typedef AnalyticalSingle<traits_type, spherical_shell_type> spherical_single_type;
    typedef AnalyticalSingle<traits_type, cylindrical_shell_type> cylindrical_single_type;
    typedef AnalyticalPair<traits_type, spherical_shell_type> spherical_pair_type;
    typedef AnalyticalPair<traits_type, cylindrical_shell_type> cylindrical_pair_type;

    typedef boost::variant<boost::none_t, spherical_shell_type, cylindrical_shell_type> shell_variant_type;

    enum domain_kind
    {
        NONE = 0,
        SPHERICAL_SINGLE,
        CYLINDRICAL_SINGLE,
        SPHERICAL_PAIR,
        CYLINDRICAL_PAIR,
        MULTI,
        NUM_DOMAIN_KINDS
    };

    enum single_event_kind
    {
        SINGLE_EVENT_REACTION,
        SINGLE_EVENT_ESCAPE,
        NUM_SINGLE_EVENT_KINDS
    };

    enum pair_event_kind
    {
        PAIR_EVENT_SINGLE_REACTION_0,
        PAIR_EVENT_SINGLE_REACTION_1,
        PAIR_EVENT_COM_ESCAPE,
        PAIR_EVENT_IV_UNDETERMINED,
        PAIR_EVENT_IV_ESCAPE,
        PAIR_EVENT_IV_REACTION,
        NUM_PAIR_EVENT_KINDS
    };

protected:
    typedef boost::fusion::map<
        boost::fusion::pair<spherical_shell_type, 
                            MatrixSpace<spherical_shell_type,
                                        shell_id_type, ecell4::utils::get_mapper_mf>*>,
        boost::fusion::pair<cylindrical_shell_type, MatrixSpace<cylindrical_shell_type,
                                        shell_id_type, ecell4::utils::get_mapper_mf>*> >
            shell_matrix_map_type;
    typedef typename boost::remove_pointer<
        typename boost::fusion::result_of::value_at_key<
            shell_matrix_map_type,
            spherical_shell_type>::type>::type
                spherical_shell_matrix_type;
    typedef typename boost::remove_pointer<
        typename boost::fusion::result_of::value_at_key<
            shell_matrix_map_type,
            cylindrical_shell_type>::type>::type
                cylindrical_shell_matrix_type;
    typedef typename ecell4::utils::get_mapper_mf<domain_id_type, boost::shared_ptr<domain_type> >::type domain_map;
    typedef typename network_rules_type::reaction_rules reaction_rules;
    typedef typename network_rules_type::reaction_rule_type reaction_rule_type;
    typedef typename traits_type::rate_type rate_type;

    class birth_event: public event_type
    {
    public:
        birth_event(time_type time, const reaction_rule_type& rr)
            : event_type(time), rr_(rr)
        {
            ;
        }

        virtual ~birth_event() {}

        const reaction_rule_type& reaction_rule() const
        {
            return rr_;
        }

    protected:
        reaction_rule_type rr_;
    };

    struct domain_event_base: public event_type
    {
        domain_event_base(time_type time): event_type(time) {}

        virtual domain_type& domain() const = 0;
    };

    template<typename Tdomain_, typename TeventKind_>
    class domain_event: public domain_event_base
    {
    public:
        typedef domain_event_base base_type;
        typedef Tdomain_ domain_type;
        typedef TeventKind_ event_kind_type;

    public:
        virtual domain_type& domain() const { return domain_; }

        virtual ~domain_event() {}

        event_kind_type kind() const { return kind_; }

        domain_event(time_type time,
                     domain_type& domain,
                     event_kind_type kind)
            : base_type(time), domain_(domain), kind_(kind) {}

    private:
        domain_type& domain_;
        event_kind_type kind_;
    };

    typedef domain_event<single_type, single_event_kind> single_event;
    typedef domain_event<pair_type, pair_event_kind> pair_event;

    class multi_event: public domain_event_base
    {
    public:
        typedef domain_event_base base_type;
        typedef multi_type domain_type;

    public:
        virtual domain_type& domain() const { return domain_; }

        virtual ~multi_event() {}

        multi_event(time_type time,
                     domain_type& domain)
            : base_type(time), domain_(domain) {}

    private:
        domain_type& domain_;
    };

    struct intruder_collector
    {
        intruder_collector(world_type const& world,
                           particle_shape_type const& cmp,
                           domain_id_type const& ignore)
            : world(world), cmp(cmp), ignore(ignore),
              closest(domain_id_type(),
                      std::numeric_limits<length_type>::infinity()) {}

        template<typename Titer>
        void operator()(Titer const& i, position_type const& off)
        {
            domain_id_type const& did((*i).second.did());
            if (did == ignore)
                return;

            length_type const distance(
                    world.distance(
                        shape(offset((*i).second, off)), cmp.position()));
            if (distance > cmp.radius())
            {
                if (distance < closest.second)
                {
                    closest.first = did;
                    closest.second = distance;
                }
            }
            else
            {
                if (!intruders.container())
                {
                    intruders.container().set(new std::vector<domain_id_type>());
                }
                intruders.push_no_duplicate(did);
            }
        }

        world_type const& world;
        particle_shape_type cmp;
        domain_id_type ignore;
        std::pair<domain_id_type, length_type> closest;
        sorted_list<std::vector<domain_id_type>,
                    std::less<domain_id_type>,
                    pointer_as_ref<std::vector<domain_id_type> > > intruders;
    };

    struct no_filter
    {
        bool operator()(domain_id_type const&) const { return true; }
    };

    struct one_id_filter
    {
        bool operator()(domain_id_type const& did) const { return did != ignore; }
        one_id_filter(domain_id_type const& did): ignore(did) {}

        domain_id_type const ignore;
    };

    template<typename TfilterFn_>
    struct domain_collector 
    {
        typedef TfilterFn_ filter_function;

        domain_collector(world_type const& world,
                        particle_shape_type const& cmp,
                        filter_function const& filter)
            : world(world), cmp(cmp), filter(filter) {}

        template<typename Titer>
        void operator()(Titer const& i, position_type const& off)

        {
            domain_id_type const& did((*i).second.did());
            if (!filter(did))
                return;

            length_type const distance(world.distance(shape(offset((*i).second, off)), cmp.position()));
            if (distance < cmp.radius())
            {
                if (!neighbors.container())
                {
                    neighbors.container().set(new std::vector<domain_id_type>());
                }
                neighbors.push_no_duplicate(did);
            }
        }

        world_type const& world;
        particle_shape_type cmp;
        filter_function const& filter;

        std::pair<domain_id_type, length_type> closest;
        sorted_list<std::vector<domain_id_type>,
                    std::less<domain_id_type>,
                    pointer_as_ref<std::vector<domain_id_type> > > neighbors;
    };

    template<typename TdidSet_>
    struct closest_object_finder
    {
        typedef TdidSet_ domain_id_set;

        closest_object_finder(world_type const& world,
                              position_type const& cmp,
                              domain_id_set const& ignore)
            : world(world), cmp(cmp), ignore(ignore),
              closest(domain_id_type(),
                      std::numeric_limits<length_type>::infinity()) {}

        template<typename Titer>
        void operator()(Titer const& i, position_type const& off)
        {
            domain_id_type const& did((*i).second.did());
            if (collection_contains(ignore, did))
                return;

            length_type const distance(world.distance(shape(offset((*i).second, off)), cmp));
            if (distance < closest.second)
            {
                closest.first = did;
                closest.second = distance;
            }
        }

        world_type const& world;
        position_type cmp;
        domain_id_set const& ignore;
        std::pair<domain_id_type, length_type> closest;
    };

    template<typename Tcol_>
    struct shell_collector_applier
    {
        typedef Tcol_ collector_type;

        shell_collector_applier(collector_type& col,
                                   position_type const& pos)
            : col_(col), pos_(pos) {}

        template<typename T>
        void operator()(T const& smat) const
        {
            world_type::traits_type::each_neighbor(*smat.second, col_, pos_);
        }

    private:
        collector_type& col_;
        position_type pos_;
    };

    template<typename Tmap_>
    struct domain_shell_map_builder
    {
        typedef Tmap_ domain_shell_association;
        domain_shell_map_builder(world_type const& world,
                                 domain_shell_association& did_map)
            : world_(world), did_map_(did_map) {}

        template<typename T>
        void operator()(T const& smat) const
        {
            BOOST_ASSERT(world_.edge_lengths() == (*smat.second).edge_lengths());
            BOOST_FOREACH (typename boost::remove_pointer<typename T::second_type>::type::value_type pair, *smat.second)
            {
                did_map_[pair.second.did()].insert(pair.first);
            }
        }

    private:
        world_type const& world_;
        domain_shell_association& did_map_;
    };

    template<typename Tset_>
    struct shell_id_collector
    {
        shell_id_collector(Tset_& shell_ids)
            : shell_ids_(shell_ids) {}

        template<typename T>
        void operator()(T const& smat) const
        {
            std::for_each(
                boost::begin(*smat.second),
                boost::end(*smat.second),
                compose_unary(
                    boost::bind(&insert<Tset_>,
                                boost::reference_wrapper<Tset_>(shell_ids_),
                                _1),
                    select_first<typename boost::remove_pointer<
                        typename T::second_type>::type::value_type>()));
        }

    private:
        Tset_& shell_ids_;
    };

    struct shell_finder
    {
        shell_finder(shell_id_type const& id, shell_variant_type& result)
            : id(id), result(result) {}

        template<typename T>
        void operator()(T const& smat) const
        {
            typename boost::remove_pointer<typename T::second_type>::type::const_iterator i((*smat.second).find(id));
            if (i != (*smat.second).end())
            {
                result = (*i).second;
            }
        }

        shell_id_type id;
        shell_variant_type& result;
    };

    struct draw_on_com_escape
    {
        position_type draw_com(spherical_pair_type const& domain,
                               time_type dt) const
        {
            return add(
                domain.shell().second.position(),
                normalize(
                    create_vector<position_type>(
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.)),
                    domain.a_R()));
            //XXX: ERROR!!! Use rng_.direction3d
        }

        position_type draw_iv(spherical_pair_type const& domain,
                              time_type dt, position_type const& old_iv) const
        {
            boost::scoped_ptr<greens_functions::PairGreensFunction> const gf(
                choose_pair_greens_function(domain, dt));
            length_type const r(draw_r(
                rng_, *gf, dt, domain.a_r(), domain.sigma()));
            length_type const theta(draw_theta(rng_, *gf, dt, r));
            return adjust_iv_with_old_iv(
                spherical_to_cartesian(
                    position_type(r, theta, rng_.uniform(0., 1.) * 2 * M_PI)),
                old_iv);
        }

        position_type draw_com(cylindrical_pair_type const& domain,
                               time_type dt) const
        {
            throw not_implemented("unsupported pair type.");
            // boost::shared_ptr<structure_type> const _structure(
            //     world_.get_structure(
            //         world_.find_molecule_info(
            //             domain.particles()[0].second.sid())
            //         .structure_id));

            // cylindrical_surface_type const* const structure(
            //     dynamic_cast<cylindrical_surface_type*>(_structure.get()));

            // return add(
            //     domain.shell().second.position(),
            //     multiply(structure->shape().unit_z(), domain.a_R()));
        }

        position_type draw_iv(cylindrical_pair_type const& domain,
                              time_type dt, position_type const& old_iv) const
        {
            BOOST_ASSERT(::size(domain.reactions()) == 1);
            throw not_implemented("unsupported pair type.");
            // length_type const r(
            //     draw_r(rng_, greens_functions::GreensFunction3DRadAbs(domain.D_tot(),
            //         domain.reactions()[0].k(), domain.r0(),
            //         domain.sigma(), domain.a_r()),
            //        dt, domain.a_r(), domain.sigma()));
            // return multiply(normalize(old_iv), r);
        }

        draw_on_com_escape(rng_type& rng, world_type const& world)
            : rng_(rng), world_(world) {}

        rng_type& rng_;
        world_type const& world_;
    };

    // struct draw_on_single_reaction
    // {
    //     position_type draw_com(spherical_pair_type const& domain,
    //                            time_type dt) const
    //     {
    //         return add(
    //             domain.shell().second.position(),
    //             draw_r(rng_,
    //                     greens_functions::GreensFunction3DAbsSym(domain.D_R(), domain.a_R()),
    //                     dt, domain.a_R()));

    //     }

    //     position_type draw_iv(spherical_pair_type const& domain,
    //                           time_type dt, position_type const& old_iv) const
    //     {
    //         boost::scoped_ptr<greens_functions::PairGreensFunction> const gf(
    //             choose_pair_greens_function(domain, dt));
    //         length_type const r(draw_r(
    //             rng_, *gf, dt, domain.a_r(), domain.sigma()));
    //         length_type const theta(draw_theta(rng_, *gf, dt, r));
    //         return adjust_iv_with_old_iv(
    //             spherical_to_cartesian(
    //                 position_type(r, theta, rng_.uniform(0., 1.) * 2 * M_PI)),
    //             old_iv);
    //     }

    //     position_type draw_com(cylindrical_pair_type const& domain,
    //                            time_type dt) const
    //     {
    //         boost::shared_ptr<structure_type> const _structure(
    //             world_.find_molecule_info(
    //                 domain.particles()[0].second.sid())
    //             .structure_id);
    //
    //         cylindrical_surface_type const* const structure(
    //             dynamic_cast<cylindrical_surface_type*>(_structure.get()));

    //         BOOST_ASSERT(structure);

    //         return add(
    //             domain.shell().second.position(),
    //             multiply(structure->shape().unit_z(), domain.a_R()));
    //     }

    //     position_type draw_iv(cylindrical_pair_type const& domain,
    //                           time_type dt, position_type const& old_iv) const
    //     {
    //         BOOST_ASSERT(::size(domain.reactions()) == 1);
    //         length_type const r(
    //             draw_r(rng_, greens_functions::GreensFunction3DRadAbs(domain.D_tot(),
    //                 domain.reactions()[0].k(), domain.r0(),
    //                 domain.sigma(), domain().a_r()),
    //                dt, domain.a_r(), domain.sigma()));
    //         BOOST_ASSERT(r > domain.sigma() && r <= domain.a_r());
    //         return multiply(normalize(old_iv), r);
    //     }

    //     draw_on_single_reaction(rng_type& rng, world_type const& world)
    //         : rng_(rng), world_(world) {}

    //     rng_type& rng_;
    //     world_type const& world_;
    // };

    struct draw_on_iv_escape
    {
        position_type draw_com(spherical_pair_type const& domain,
                               time_type dt)
        {
            return add(
                domain.shell().second.position(),
                normalize(
                    create_vector<position_type>(
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.)),
                    draw_r(
                        rng_,
                        greens_functions::GreensFunction3DAbsSym(domain.D_R(), domain.a_R()),
                        dt, domain.a_R())));
        }

        position_type draw_iv(spherical_pair_type const& domain,
                              time_type dt, position_type const& old_iv)
        {
            boost::scoped_ptr<greens_functions::PairGreensFunction> const gf(
                choose_pair_greens_function(domain, dt));
            length_type const r(domain.a_r());
            length_type const theta(draw_theta(rng_, *gf, dt, r));
            return adjust_iv_with_old_iv(
                spherical_to_cartesian(
                    position_type(r, theta, rng_.uniform(0., 1.) * 2 * M_PI)),
                old_iv);
        }

        position_type draw_com(cylindrical_pair_type const& domain,
                               time_type dt)
        {
            throw not_implemented("unsupported pair type.");
            // boost::shared_ptr<structure_type> const _structure(
            //     world_.get_structure(
            //         world_.find_molecule_info(
            //             domain.particles()[0].second.sid())
            //         .structure_id));

            // cylindrical_surface_type const* const structure(
            //     dynamic_cast<cylindrical_surface_type*>(_structure.get()));

            // BOOST_ASSERT(structure);

            // length_type const r_R(draw_r(
            //     rng_,
            //     greens_functions::GreensFunction3DAbsSym(domain.D_R(), domain.a_R()),
            //     dt, domain.a_R()));
            // return add(
            //     domain.shell().second.position(),
            //     multiply(structure->shape().unit_z(), r_R));
        }

        position_type draw_iv(cylindrical_pair_type const& domain,
                              time_type dt, position_type const& old_iv)
        {
            throw not_implemented("unsupported pair type.");
            // return multiply(normalize(old_iv), domain.a_r());
        }

        draw_on_iv_escape(rng_type& rng, world_type const& world)
            : rng_(rng), world_(world) {}

        rng_type& rng_;
        world_type const& world_;
    };

    struct draw_on_iv_reaction
    {
        position_type draw_com(spherical_pair_type const& domain,
                               time_type dt)
        {
            return add(
                domain.shell().second.position(),
                normalize(
                    create_vector<position_type>(
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.)),
                    draw_r(
                        rng_,
                        greens_functions::GreensFunction3DAbsSym(domain.D_R(), domain.a_R()),
                        dt, domain.a_R())));
        }

        position_type draw_iv(spherical_pair_type const& domain,
                              time_type dt, position_type const& old_iv)
        {
            boost::scoped_ptr<greens_functions::PairGreensFunction> const gf(
                choose_pair_greens_function(domain, dt));
            length_type const r(domain.sigma());
            length_type const theta(draw_theta(rng_, *gf, dt, r));
            return adjust_iv_with_old_iv(
                spherical_to_cartesian(
                    position_type(r, theta, rng_.uniform(0., 1.) * 2 * M_PI)),
                old_iv);
        }

        position_type draw_com(cylindrical_pair_type const& domain,
                               time_type dt)
        {
            throw not_implemented("unsupported pair type.");
            // boost::shared_ptr<structure_type> const _structure(
            //     world_.get_structure(
            //         world_.find_molecule_info(
            //             domain.particles()[0].second.sid()).structure_id));
            //
            // cylindrical_surface_type const* const structure(
            //     dynamic_cast<cylindrical_surface_type*>(_structure.get()));

            // BOOST_ASSERT(structure);

            // length_type const r_R(draw_r(
            //     rng_,
            //     greens_functions::GreensFunction3DAbsSym(domain.D_R(), domain.a_R()),
            //     dt, domain.a_R()));
            // return add(
            //     domain.shell().second.position(),
            //     multiply(structure->shape().unit_z(), r_R));
        }

        position_type draw_iv(cylindrical_pair_type const& domain,
                              time_type dt, position_type const& old_iv)
        {
            throw not_implemented("unsupported pair type.");
            // return multiply(domain.sigma(), normalize(old_iv));
        }

        draw_on_iv_reaction(rng_type& rng, world_type const& world)
            : rng_(rng), world_(world) {}

        rng_type& rng_;
        world_type const& world_;
    };

    struct draw_on_burst
    {
        position_type draw_com(spherical_pair_type const& domain,
                               time_type dt)
        {
            return add(
                domain.shell().second.position(),
                normalize(
                    create_vector<position_type>(
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.)),
                    draw_r(
                        rng_,
                        greens_functions::GreensFunction3DAbsSym(domain.D_R(), domain.a_R()),
                        dt, domain.a_R())));
        }

        position_type draw_iv(spherical_pair_type const& domain,
                              time_type dt, position_type const& old_iv)
        {
            boost::scoped_ptr<greens_functions::PairGreensFunction> const gf(
                choose_pair_greens_function(domain, dt));
            length_type const r(draw_r(
                rng_, *gf, dt, domain.a_r(), domain.sigma()));
            length_type const theta(draw_theta(rng_, *gf, dt, r));
            return adjust_iv_with_old_iv(
                spherical_to_cartesian(
                    position_type(r, theta, rng_.uniform(0., 1.) * 2 * M_PI)),
                old_iv);
        }

        position_type draw_com(cylindrical_pair_type const& domain,
                               time_type dt)
        {
            throw not_implemented("unsupported pair type.");
            // boost::shared_ptr<structure_type> const _structure(
            //     world_.get_structure(
            //         world_.find_molecule_info(
            //             domain.particles()[0].second.sid())
            //         .structure_id));

            // cylindrical_surface_type const* const structure(
            //     dynamic_cast<cylindrical_surface_type*>(_structure.get()));

            // BOOST_ASSERT(structure);

            // length_type const r_R(draw_r(
            //     rng_,
            //     greens_functions::GreensFunction3DAbsSym(domain.D_R(), domain.a_R()),
            //     dt, domain.a_R()));
            // return add(
            //     domain.shell().second.position(),
            //     multiply(structure->shape().unit_z(), r_R));
        }

        position_type draw_iv(cylindrical_pair_type const& domain,
                              time_type dt, position_type const& old_iv)
        {
            throw not_implemented("unsupported pair type.");
            // BOOST_ASSERT(::size(domain.reactions()) == 1);
            // length_type const r(
            //     draw_r(rng_,
            //         greens_functions::GreensFunction3DRadAbs(
            //             domain.D_tot(),
            //             domain.reactions()[0].k(), domain.r0(),
            //             domain.sigma(), domain.a_r()),
            //         dt, domain.a_r(), domain.sigma()));
            // return multiply(normalize(old_iv), r);
        }

        draw_on_burst(rng_type& rng, world_type const& world)
            : rng_(rng), world_(world) {}

        rng_type& rng_;
        world_type const& world_;
    };
public:
    typedef abstract_limited_generator<domain_id_pair> domain_id_pair_generator;

public:
    virtual ~EGFRDSimulator()
    {
        //std::for_each(domains_.begin(), domains_.end(),
        //    compose_unary(delete_ptr<domain_type>(),
        //                  select_second<typename domain_map::value_type>()));
    }

    EGFRDSimulator(
        const boost::shared_ptr<world_type>& world,
        const boost::shared_ptr<model_type>& ecell4_model,
        Real bd_dt_factor = 1e-5, int dissociation_retry_moves = 1,
        length_type user_max_shell_size = std::numeric_limits<length_type>::infinity())
        : base_type(world, ecell4_model),
          bd_dt_factor_(bd_dt_factor),
          num_retries_(dissociation_retry_moves),
          user_max_shell_size_(user_max_shell_size),
          ssmat_(new spherical_shell_matrix_type((*world).edge_lengths(), (*world).matrix_sizes())),
          csmat_(new cylindrical_shell_matrix_type((*world).edge_lengths(), (*world).matrix_sizes())),
          smatm_(boost::fusion::pair<spherical_shell_type,
                                     spherical_shell_matrix_type*>(ssmat_.get()),
                 boost::fusion::pair<cylindrical_shell_type,
                                     cylindrical_shell_matrix_type*>(csmat_.get())),
          single_shell_factor_(.1),
          multi_shell_factor_(.05),
          rejected_moves_(0), zero_step_count_(0), dirty_(true)
    {
        std::fill(domain_count_per_type_.begin(), domain_count_per_type_.end(), 0);
        std::fill(single_step_count_.begin(), single_step_count_.end(), 0);
        std::fill(pair_step_count_.begin(), pair_step_count_.end(), 0);
        std::fill(multi_step_count_.begin(), multi_step_count_.end(), 0);
    }

    EGFRDSimulator(
        const boost::shared_ptr<world_type>& world,
        Real bd_dt_factor = 1e-5, int dissociation_retry_moves = 1,
        length_type user_max_shell_size = std::numeric_limits<length_type>::infinity())
        : base_type(world),
          bd_dt_factor_(bd_dt_factor),
          num_retries_(dissociation_retry_moves),
          user_max_shell_size_(user_max_shell_size),
          ssmat_(new spherical_shell_matrix_type((*world).edge_lengths(), (*world).matrix_sizes())),
          csmat_(new cylindrical_shell_matrix_type((*world).edge_lengths(), (*world).matrix_sizes())),
          smatm_(boost::fusion::pair<spherical_shell_type,
                                     spherical_shell_matrix_type*>(ssmat_.get()),
                 boost::fusion::pair<cylindrical_shell_type,
                                     cylindrical_shell_matrix_type*>(csmat_.get())),
          single_shell_factor_(.1),
          multi_shell_factor_(.05),
          rejected_moves_(0), zero_step_count_(0), dirty_(true)
    {
        std::fill(domain_count_per_type_.begin(), domain_count_per_type_.end(), 0);
        std::fill(single_step_count_.begin(), single_step_count_.end(), 0);
        std::fill(pair_step_count_.begin(), pair_step_count_.end(), 0);
        std::fill(multi_step_count_.begin(), multi_step_count_.end(), 0);
    }

    length_type user_max_shell_size() const
    {
        return user_max_shell_size_;
    }

    length_type max_shell_size() const
    {
        const position_type& cell_sizes((*base_type::world_).cell_sizes());
        const length_type min_cell_size(
            std::min(cell_sizes[0], std::min(cell_sizes[1], cell_sizes[2])));
        return std::min(min_cell_size / 2 /
                        traits_type::SAFETY,
                   user_max_shell_size_);
    }

    template<typename Tshape>
    std::pair<const shell_id_type,
              typename traits_type::template shell_generator<Tshape>::type>
    new_shell(domain_id_type const& did, Tshape const& shape)
    {
        typedef typename traits_type::template shell_generator<Tshape>::type shell_type;
        std::pair<const shell_id_type, shell_type> const retval(shidgen_(), shell_type(did, shape));
        boost::fusion::at_key<shell_type>(smatm_)->update(retval);
        return retval;
    }

    template<typename T>
    std::pair<const shell_id_type, T> const& get_shell(shell_id_type const& id) const
    {
        typedef typename boost::remove_pointer<
            typename boost::fusion::result_of::value_at_key<
                shell_matrix_map_type, T>::type>::type shell_matrix_type;

        shell_matrix_type const& smat(*boost::fusion::at_key<T>(smatm_));
        
        typename shell_matrix_type::const_iterator i(smat.find(id));
        if (i == smat.end())
        {
            throw not_found(
                (boost::format("shell id #%s not found") % boost::lexical_cast<std::string>(id)).str());
        }

        return *i;
    }

    std::pair<shell_id_type, shell_variant_type> get_shell(shell_id_type const& id)
    {
        shell_variant_type result;
        boost::fusion::for_each(smatm_, shell_finder(id, result));
        return std::make_pair(id, result);
    }

    boost::shared_ptr<domain_type> get_domain(domain_id_type const& id) const
    {
        typename domain_map::const_iterator i(domains_.find(id));

        if (i == domains_.end())
        {
            throw not_found(
                (boost::format("domain id #%s not found") % boost::lexical_cast<std::string>(id)).str());
        }

        return (*i).second;
    }

    domain_id_pair_generator* get_domains() const
    {
        return make_range_generator<domain_id_pair>(domains_);
    }

    typename domain_map::size_type num_domains() const
    {
        return domains_.size();
    }

    int num_domains_per_type(domain_kind kind) const
    {
        return domain_count_per_type_[kind];
    }

    int num_single_steps_per_type(single_event_kind kind) const
    {
        return single_step_count_[kind];
    }

    int num_pair_steps_per_type(pair_event_kind kind) const
    {
        return pair_step_count_[kind];
    }

    int num_multi_steps_per_type(typename multi_type::event_kind kind) const
    {
        return multi_step_count_[kind];
    }

    std::vector<domain_id_type>*
    get_neighbor_domains(particle_shape_type const& p)
    {
        typedef domain_collector<no_filter> collector_type;
        no_filter f;
        collector_type col((*base_type::world_), p, f);
        boost::fusion::for_each(smatm_, shell_collector_applier<collector_type>(col, p.position()));
        return col.neighbors.container().get();
    }

    std::vector<domain_id_type>*
    get_neighbor_domains(particle_shape_type const& p, domain_id_type const& ignore)
    {
        typedef domain_collector<one_id_filter> collector_type;
        one_id_filter f(ignore);
        collector_type col((*base_type::world_), p, f);
        boost::fusion::for_each(smatm_, shell_collector_applier<collector_type>(col, p.position()));
        return col.neighbors.container().get();
    }

    virtual void initialize()
    {
        const position_type& edge_lengths((*base_type::world_).edge_lengths());
        const typename world_type::matrix_sizes_type&
            matrix_sizes((*base_type::world_).matrix_sizes());

        domains_.clear();
        (*ssmat_).clear();
        (*csmat_).clear();
        scheduler_.clear();

        if (edge_lengths != (*ssmat_).edge_lengths()
            || matrix_sizes != (*ssmat_).matrix_sizes())
        {
            boost::scoped_ptr<spherical_shell_matrix_type>
                newssmat(new spherical_shell_matrix_type(edge_lengths, matrix_sizes));
            boost::scoped_ptr<cylindrical_shell_matrix_type>
                newcsmat(new cylindrical_shell_matrix_type(edge_lengths, matrix_sizes));
            ssmat_.swap(newssmat);
            csmat_.swap(newcsmat);
            boost::fusion::at_key<spherical_shell_type>(smatm_) = ssmat_.get();
            boost::fusion::at_key<cylindrical_shell_type>(smatm_) = csmat_.get();
        }

        BOOST_FOREACH (particle_id_pair const& pp,
                       (*base_type::world_).get_particles_range())
        {
            boost::shared_ptr<single_type> single(create_single(pp));
            add_event(*single, SINGLE_EVENT_ESCAPE);
        }

        BOOST_FOREACH (reaction_rule_type const& rr,
                       (*base_type::network_rules_).zeroth_order_reaction_rules())
        {
            add_event(rr);
        }
        dirty_ = false;
    }

    /**
     * override
     * HERE
     */

    virtual Real next_time() const
    {
        return scheduler_.next_time();
    }

    virtual Real dt() const
    {
        return scheduler_.next_time() - base_type::t();
    }

    /**
     * override
     * THERE
     */

    virtual void step()
    {
        if (dirty_)
        {
            initialize();
        }

        _step();
    }

    void finalize()
    {
        std::vector<domain_id_type> non_singles;

        // first burst all Singles.
        BOOST_FOREACH (event_id_pair_type const& event, scheduler_.events())
        {
            {
                single_event const* single_ev(
                        dynamic_cast<single_event const*>(event.second.get()));
                if (single_ev)
                {
                    burst(single_ev->domain());
                    continue;
                }
            }
            {
                birth_event const* birth_ev(
                        dynamic_cast<birth_event const*>(event.second.get()));
                if (birth_ev)
                {
                    continue;
                }
            }
            {
                domain_event_base const* domain_ev(
                    dynamic_cast<domain_event_base const*>(event.second.get()));
                BOOST_ASSERT(domain_ev);
                non_singles.push_back(domain_ev->domain().id());
                // continue;
            }
        }

        // then burst all Pairs and Multis.
        burst_domains(non_singles);

        base_type::dt_ = 0.;
    }

    // virtual bool step(time_type upto)
    virtual bool step(const time_type& upto)
    {
        if (dirty_)
        {
            initialize();
        }

        if (upto <= this->t())
        {
            return false;
        }

        // if (upto >= scheduler_.top().second->time())
        if (upto >= scheduler_.next_time())
        {
            _step();
            return true;
        }

        LOG_INFO(("stop at %.16g", upto));

        this->set_t(upto);
        this->finalize();
        return false;
    }

    // {{{ clear_volume
    // called by Multi
    void clear_volume(particle_shape_type const& p)
    {
        boost::scoped_ptr<std::vector<domain_id_type> > domains(
            get_neighbor_domains(p));
        if (domains)
        {
            burst_domains(*domains);
        }
    }

    void clear_volume(particle_shape_type const& p,
                      domain_id_type const& ignore)
    {
        boost::scoped_ptr<std::vector<domain_id_type> > domains(
            get_neighbor_domains(p, ignore));

        if (domains)
        {
            burst_domains(*domains);
        }
    }
    // }}}

    bool check() const
    {
        LOG_INFO(("checking overall consistency"));

        bool retval(true);

        CHECK(this->t() >= 0.0);
        CHECK(base_type::dt_ >= 0.0);

        typedef std::map<domain_id_type, std::set<shell_id_type> >
            domain_shell_association;

        domain_shell_association did_map;

        boost::fusion::for_each(smatm_,
                domain_shell_map_builder<domain_shell_association>(
                    (*base_type::world_), did_map));

        std::set<domain_id_type> scheduled_domains;
        typename domain_type::size_type shells_correspond_to_domains(0);
        std::size_t particles_correspond_to_domains(0);

        BOOST_FOREACH (typename event_scheduler_type::value_type const& value,
                       scheduler_.events())
        {
            domain_type const& domain(dynamic_cast<domain_event_base&>(*value.second).domain());
            CHECK(check_domain(domain));

            if (!scheduled_domains.insert(domain.id()).second)
            {
                LOG_WARNING(("domain id %s is doubly scheduled!", boost::lexical_cast<std::string>(domain.id()).c_str()));
            }

            CHECK(domain.event() == value);

            typename domain_type::size_type const num_shells(domain.num_shells());

            shells_correspond_to_domains += num_shells;
            particles_correspond_to_domains += domain.multiplicity();

            typename std::set<shell_id_type> const& shell_ids(
                did_map[domain.id()]);
            CHECK(static_cast<typename domain_type::size_type>(
                    ::size(shell_ids)) == num_shells);
        }

        CHECK((*base_type::world_).num_particles() == particles_correspond_to_domains);

        {
            std::vector<domain_id_type> diff;
            ::difference(make_select_first_range(did_map), scheduled_domains,
                    std::back_inserter(diff));

            if (diff.size() != 0)
            {
                LOG_WARNING(("domains not scheduled: %s",
                    stringize_and_join(diff, ", ").c_str()));
                BOOST_FOREACH (domain_id_type const& domain_id, diff)
                {
                    LOG_WARNING(("  shells that belong to unscheduled domain %s: %s",
                        boost::lexical_cast<std::string>(domain_id).c_str(),
                        stringize_and_join(did_map[domain_id], ", ").c_str()));
                }
                retval = false;
            }
        }

        std::set<shell_id_type> all_shell_ids;
        boost::fusion::for_each(smatm_,
                shell_id_collector<std::set<shell_id_type> >(all_shell_ids));

        if (shells_correspond_to_domains != static_cast<std::size_t>(::size(all_shell_ids)))
        {
            LOG_WARNING(("shells_correspond_to_domains=%zu, shell_population=%zu", shells_correspond_to_domains, static_cast<std::size_t>(::size(all_shell_ids))));
            dump_events();
            retval = false;
        }
        return retval;
    }

    virtual bool check_reaction() const
    {
        return last_reactions().size() > 0;
    }

    std::vector<std::pair<ecell4::ReactionRule, reaction_info_type> > last_reactions() const
    {
        return (*dynamic_cast<ReactionRecorderWrapper<reaction_record_type>*>(
            base_type::rrec_.get())).last_reactions();
    }

protected:
    template<typename Tshell>
    void move_shell(std::pair<const shell_id_type, Tshell> const& shell)
    {
        typedef Tshell shell_type;
        boost::fusion::at_key<shell_type>(smatm_)->update(shell);
    }

    template<typename T>
    void update_shell_matrix(AnalyticalSingle<traits_type, T> const& domain)
    {
        move_shell(domain.shell());
    }

    template<typename T>
    void update_shell_matrix(AnalyticalPair<traits_type, T> const& domain)
    {
        move_shell(domain.shell());
    }

    void update_shell_matrix(shaped_domain_type const& domain)
    {
        {
            spherical_single_type const* _domain(dynamic_cast<spherical_single_type const*>(&domain));
            if (_domain) {
                update_shell_matrix(*_domain);
                return;
            }
        }
        {
            cylindrical_single_type const* _domain(dynamic_cast<cylindrical_single_type const*>(&domain));
            if (_domain) {
                update_shell_matrix(*_domain);
                return;
            }
        }
        {
            spherical_pair_type const* _domain(dynamic_cast<spherical_pair_type const*>(&domain));
            if (_domain) {
                update_shell_matrix(*_domain);
                return;
            }
        }
        {
            cylindrical_pair_type const* _domain(dynamic_cast<cylindrical_pair_type const*>(&domain));
            if (_domain) {
                update_shell_matrix(*_domain);
                return;
            }
        }
        throw not_implemented("unsupported domain type");
    }

    // remove_domain_but_shell {{{
    template<typename T>
    void remove_domain_but_shell(AnalyticalSingle<traits_type, T>& domain)
    {
        --domain_count_per_type_[get_domain_kind(domain)];
        _remove_domain_but_shell(domain);
    }

    template<typename T>
    void remove_domain_but_shell(AnalyticalPair<traits_type, T>& domain)
    {
        --domain_count_per_type_[get_domain_kind(domain)];
        _remove_domain_but_shell(domain);
    }

    void remove_domain_but_shell(multi_type& domain)
    {
        --domain_count_per_type_[get_domain_kind(domain)];
        _remove_domain_but_shell(domain);
    }

    void remove_domain_but_shell(domain_type& domain)
    {
        --domain_count_per_type_[get_domain_kind(domain)];
        _remove_domain_but_shell(domain);
    }

    void _remove_domain_but_shell(domain_type& domain)
    {
        LOG_DEBUG(("remove_domain_but_shell: %s", boost::lexical_cast<std::string>(domain.id()).c_str()));
        event_id_type const event_id(domain.event().first);

        // domains_.erase(domain.id()); // this hits a bug in gcc 4.4 (at least)'s unordered_map.
        typename domain_map::iterator domain_to_be_removed(domains_.find(domain.id()));
        if (base_type::paranoiac_)
        {
            BOOST_ASSERT(domain_to_be_removed != domains_.end());
        }
        domains_.erase(domain_to_be_removed);

        try
        {
            remove_event(event_id);
        }
        catch (std::out_of_range const&)
        {
            LOG_DEBUG(("event %s already removed; ignoring.", boost::lexical_cast<std::string>(event_id).c_str()));
        }
    }

    // }}}

    // remove_domain {{{
    template<typename T>
    void remove_domain(AnalyticalSingle<traits_type, T>& domain)
    {
        typedef T shell_type;
        boost::fusion::at_key<shell_type>(smatm_)->erase(domain.shell().first);
        remove_domain_but_shell(domain);
    }

    template<typename T>
    void remove_domain(AnalyticalPair<traits_type, T>& domain)
    {
        typedef T shell_type;
        boost::fusion::at_key<shell_type>(smatm_)->erase(domain.shell().first);
        remove_domain_but_shell(domain);
    }

    void remove_domain(multi_type& domain)
    {
        BOOST_FOREACH (spherical_shell_id_pair const& shell, domain.get_shells())
        {
            boost::fusion::at_key<spherical_shell_type>(smatm_)->erase(shell.first);
        }
        remove_domain_but_shell(domain);
    }

    void remove_domain(domain_type& domain)
    {
        {
            spherical_single_type* _domain(dynamic_cast<spherical_single_type*>(&domain));
            if (_domain)
            {
                remove_domain(*_domain);
                return;
            }
        }
        {
            cylindrical_single_type* _domain(dynamic_cast<cylindrical_single_type*>(&domain));
            if (_domain)
            {
                remove_domain(*_domain);
                return;
            }
        }
        {
            spherical_pair_type* _domain(dynamic_cast<spherical_pair_type*>(&domain));
            if (_domain)
            {
                remove_domain(*_domain);
                return;
            }
        }
        {
            cylindrical_pair_type* _domain(dynamic_cast<cylindrical_pair_type*>(&domain));
            if (_domain)
            {
                remove_domain(*_domain);
                return;
            }
        }
        {
            multi_type* _domain(dynamic_cast<multi_type*>(&domain));
            if (_domain)
            {
                remove_domain(*_domain);
                return;
            }
        }
        throw not_implemented(std::string("unsupported domain type"));
    }
    // }}}

    void add_event(single_type& domain, single_event_kind const& kind)
    {
        if (base_type::paranoiac_)
            BOOST_ASSERT(domains_.find(domain.id()) != domains_.end());

        boost::shared_ptr<event_type> new_event(
            new single_event(this->t() + domain.dt(), domain, kind));
        domain.event() = std::make_pair(scheduler_.add(new_event), new_event);
        LOG_DEBUG(("add_event: #%d - %s", domain.event().first, boost::lexical_cast<std::string>(domain).c_str()));
    }

    void add_event(pair_type& domain, pair_event_kind const& kind)
    {
        if (base_type::paranoiac_)
            BOOST_ASSERT(domains_.find(domain.id()) != domains_.end());

        boost::shared_ptr<event_type> new_event(
            new pair_event(this->t() + domain.dt(), domain, kind));
        domain.event() = std::make_pair(scheduler_.add(new_event), new_event);
        LOG_DEBUG(("add_event: #%d - %s", domain.event().first, boost::lexical_cast<std::string>(domain).c_str()));
    }

    void add_event(multi_type& domain)
    {
        if (base_type::paranoiac_)
            BOOST_ASSERT(domains_.find(domain.id()) != domains_.end());

        boost::shared_ptr<event_type> new_event(
            new multi_event(this->t() + domain.dt(), domain));
        domain.event() = std::make_pair(scheduler_.add(new_event), new_event);
        LOG_DEBUG(("add_event: #%d - %s", domain.event().first, boost::lexical_cast<std::string>(domain).c_str()));
    }

    /**
     * The following add_event function is for birth_event.
     */
    void add_event(reaction_rule_type const& rr)
    {
        const double rnd(this->rng().uniform(0, 1));
        const double dt(gsl_sf_log(1.0 / rnd) / double(rr.k() * (*base_type::world_).volume()));
        boost::shared_ptr<event_type> new_event(new birth_event(this->t() + dt, rr));
        scheduler_.add(new_event);
    }

    void remove_event(event_id_type const& id)
    {
        LOG_DEBUG(("remove_event: #%d", id));
        scheduler_.remove(id);
    }

    void remove_event(domain_type const& domain)
    {
        remove_event(domain.event().first);
    }

    // create_single {{{
    boost::shared_ptr<single_type> create_single(particle_id_pair const& p)
    {
        domain_kind kind(NONE);
        single_type* new_single(0);
        domain_id_type did(didgen_());

        struct factory: ImmutativeStructureVisitor<typename world_type::traits_type>
        {
            virtual ~factory() {}

            // virtual void operator()(spherical_surface_type const& structure) const
            // {
            //     throw not_implemented(
            //         (boost::format("unsupported structure type: %s") %
            //             boost::lexical_cast<std::string>(structure)).str());
            // }

            // virtual void operator()(cylindrical_surface_type const& structure) const
            // {
            //     // Heads up. The cylinder's *size*, not radius, is changed when
            //     // you make the cylinder bigger, because of the redefinition of
            //     // set_radius.
            //     // The radius of a rod is not more than it has to be (namely
            //     // the radius of the particle), so if the particle undergoes an
            //     // unbinding reaction, we still have to clear the target volume
            //     // and the move may be rejected (NoSpace error).
            //     const cylindrical_shell_id_pair new_shell(
            //         _this->new_shell(
            //             did, typename cylindrical_shell_type::shape_type(
            //                 p.second.position(), p.second.radius(),
            //                 structure.shape().unit_z(),
            //                 p.second.radius())));
            //     new_single = new cylindrical_single_type(did, p, new_shell);
            //     kind = CYLINDRICAL_SINGLE;
            // }

            // virtual void operator()(planar_surface_type const& structure) const
            // {
            //     cylindrical_shell_id_pair const new_shell(
            //         _this->new_shell(did, typename cylindrical_shell_type::shape_type(
            //             p.second.position(), p.second.radius(),
            //             normalize(cross_product(
            //                 structure.shape().unit_x(),
            //                 structure.shape().unit_y())),
            //             p.second.radius())));
            //     new_single = new cylindrical_single_type(did, p, new_shell);
            //     kind = CYLINDRICAL_SINGLE;
            // }

            virtual void operator()(cuboidal_region_type const& structure) const
            {
                spherical_shell_id_pair new_shell(
                    _this->new_shell(
                        did, typename spherical_shell_type::shape_type(
                            p.second.position(), p.second.radius())));
                    // _this->new_shell(did, ::shape(p.second)));
                new_single = new spherical_single_type(did, p, new_shell);
                kind = SPHERICAL_SINGLE;
            }

            factory(EGFRDSimulator* _this, particle_id_pair const& p,
                    domain_id_type const& did, single_type*& new_single,
                    domain_kind& kind)
                : _this(_this), p(p), did(did), new_single(new_single),
                  kind(kind) {}

            EGFRDSimulator* _this;
            particle_id_pair const& p;
            domain_id_type const& did;
            single_type*& new_single;
            domain_kind& kind;
        };

        molecule_info_type const species((*base_type::world_).get_molecule_info(p.second.species()));
        // molecule_info_type const& species((*base_type::world_).find_molecule_info(p.second.species()));
        dynamic_cast<particle_simulation_structure_type const&>(*(*base_type::world_).get_structure(species.structure_id)).accept(factory(this, p, did, new_single, kind));
        boost::shared_ptr<domain_type> const retval(new_single);
        domains_.insert(std::make_pair(did, retval));
        BOOST_ASSERT(kind != NONE);
        ++domain_count_per_type_[kind];
        return boost::dynamic_pointer_cast<single_type>(retval);
    }
    // }}}

    // create_pair {{{
    boost::shared_ptr<pair_type> create_pair(particle_id_pair const& p0,
                                             particle_id_pair const& p1,
                                             position_type const& com,
                                             position_type const& iv,
                                             length_type shell_size)
    {
        domain_kind kind(NONE);
        pair_type* new_pair(0);
        domain_id_type did(didgen_());

        struct factory: ImmutativeStructureVisitor<typename world_type::traits_type>
        {
            // virtual void operator()(spherical_surface_type const& structure) const
            // {
            //     throw not_implemented(
            //         (boost::format("unsupported structure type: %s") %
            //             boost::lexical_cast<std::string>(structure)).str());
            // }

            // virtual void operator()(cylindrical_surface_type const& structure) const
            // {
            //     // The radius of a rod is not more than it has to be (namely
            //     // the radius of the biggest particle), so if the particle
            //     // undergoes an unbinding reaction we still have to clear the
            //     // target volume and the move may be rejected (NoSpace error).
            //     cylindrical_shell_id_pair const new_shell(
            //         _this->new_shell(did, typename cylindrical_shell_type::shape_type(
            //             com,
            //             shell_size,
            //             shape(structure).unit_z(),
            //             std::max(p0.second.radius(), p1.second.radius()))));
            //     new_pair = new cylindrical_pair_type(did, p0, p1, new_shell,
            //                                          iv, rules);
            //     kind = CYLINDRICAL_PAIR;
            // }

        
            // virtual void operator()(planar_surface_type const& structure) const
            // {
            //     cylindrical_shell_id_pair const new_shell(
            //         _this->new_shell(did, typename cylindrical_shell_type::shape_type(
            //             com,
            //             shell_size,
            //             normalize(cross_product(
            //                 shape(structure).unit_x(),
            //                 shape(structure).unit_y())),
            //             std::max(p0.second.radius(), p1.second.radius()))));
            //     new_pair = new cylindrical_pair_type(did, p0, p1, new_shell,
            //                                            iv, rules);
            //     kind = CYLINDRICAL_PAIR;
            // }

            virtual void operator()(cuboidal_region_type const& structure) const
            {
                spherical_shell_id_pair new_shell(
                    _this->new_shell(did,
                        typename spherical_shell_type::shape_type(com, shell_size)));
                new_pair = new spherical_pair_type(did, p0, p1, new_shell,
                                                   iv, rules);
                kind = SPHERICAL_PAIR;
            }

            factory(EGFRDSimulator* _this, particle_id_pair const& p0,
                    particle_id_pair const& p1, position_type const& com,
                    position_type const& iv, length_type shell_size,
                    domain_id_type const& did, pair_type*& new_pair,
                    domain_kind& kind)
                : _this(_this), p0(p0), p1(p1), com(com), iv(iv),
                  shell_size(shell_size), did(did),
                  rules((*_this->network_rules_).query_reaction_rule(
                        p0.second.species(), p1.second.species())),
                  new_pair(new_pair), kind(kind) {}

            EGFRDSimulator* _this;
            particle_id_pair const& p0;
            particle_id_pair const& p1;
            position_type const& com;
            position_type const& iv;
            const length_type shell_size;
            domain_id_type const& did;
            typename network_rules_type::reaction_rule_vector const& rules;
            pair_type*& new_pair;
            domain_kind& kind;
        };

        molecule_info_type const species((*base_type::world_).get_molecule_info(p0.second.species()));
        // molecule_info_type const& species((*base_type::world_).find_molecule_info(p0.second.species()));
        dynamic_cast<particle_simulation_structure_type&>(*(*base_type::world_).get_structure(species.structure_id)).accept(factory(this, p0, p1, com, iv, shell_size, did, new_pair, kind));

        boost::shared_ptr<domain_type> const retval(new_pair);
        domains_.insert(std::make_pair(did, retval));
        BOOST_ASSERT(kind != NONE);
        ++domain_count_per_type_[kind];
        return boost::dynamic_pointer_cast<pair_type>(retval);
    }
    // }}}

    // create_multi {{{
    boost::shared_ptr<multi_type> create_multi()
    {
        domain_id_type did(didgen_());
        multi_type* new_multi(new multi_type(did, *this, bd_dt_factor_));
        boost::shared_ptr<domain_type> const retval(new_multi);
        domains_.insert(std::make_pair(did, retval));
        ++domain_count_per_type_[MULTI];
        return boost::dynamic_pointer_cast<multi_type>(retval);
    }
    // }}}

    // draw_r {{{
    template<typename Tgf>
    static length_type draw_r(rng_type& rng,
                              Tgf const& gf,
                              time_type dt,
                              length_type a,
                              length_type sigma = -1.)
    {
        LOG_DEBUG(("draw_r: dt=%.16g, a=%.16g, sigma=%.16g", dt, a, sigma));
        BOOST_ASSERT(a > sigma && a >= 0.);
        length_type r(0.);
        double rnd(0.);
        try
        {
            do
            {
                rnd = rng.uniform(0., 1.);
                r = gf.drawR(rnd, dt);
            } while (r > a || r <= sigma);
        }
        catch (std::exception const& e)
        {
            throw propagation_error(
                (boost::format(
                    "gf.drawR() failed: %s, gf=%s, rnd=%.16g, dt=%.16g, a=%.16g, sigma=%.16g: %s") % e.what() % gf.getName() % rnd % dt % a % sigma % gf.dump()).str());
        }

        return r;
    }
    // }}}

    // draw_theta {{{
    template<typename Tgf>
    static length_type draw_theta(rng_type& rng,
                              Tgf const& gf,
                              time_type dt,
                              length_type r)
    {
        length_type theta(0.);
        double rnd(0.);
        try
        {
            rnd = rng.uniform(0., 1.);
            theta = gf.drawTheta(rnd, r, dt);
        }
        catch (std::exception const& e)
        {
            throw propagation_error(
                (boost::format(
                    "gf.drawTheta() failed: %s, gf=%s, rnd=%.16g, dt=%.16g, r=%.16g: %s") % e.what() % gf.getName() % rnd % dt % r % gf.dump()).str());
        }

        return theta;
    }
    // }}}

    // draw_displacement {{{
    position_type draw_displacement(
        AnalyticalSingle<traits_type, spherical_shell_type> const& domain,
        length_type r)
    {
        // const double cos_theta(this->rng().uniform(-1., 1.));
        // const double sin_theta(sqrt(1 - cos_theta * cos_theta));
        // double sin_phi, cos_phi;
        // sincos(this->rng().uniform(0., 2 * M_PI), &sin_phi, &cos_phi);
        // return normalize(
        //     create_vector<position_type>(
        //         sin_theta * cos_phi, sin_theta * sin_phi, cos_theta), r);

        // double x, y, z;
        // this->rng().dir_3d(&x, &y, &z);
        // return normalize(
        //     create_vector<position_type>(x, y, z), r);

        // return this->rng().direction3d(r);
        return normalize(this->rng().direction3d(1), r);
    }

    position_type draw_displacement(
        AnalyticalSingle<traits_type, cylindrical_shell_type> const& domain,
        length_type r)
    {
        // return multiply(shape(domain.shell().second).unit_z(), r);
        return multiply(shape(domain.shell().second).axis(), r);
    }
    // }}}

    // draw_new_position {{{
    template<typename Tshell>
    position_type draw_new_position(
            AnalyticalSingle<traits_type, Tshell> const& domain,
            time_type dt)
    {
        typedef Tshell shell_type;
        typedef typename shell_type::shape_type shape_type;
        typedef typename detail::get_greens_function<shape_type>::type greens_function;
        length_type const r(
            draw_r(
                this->rng(),
                greens_function(
                    domain.particle().second.D(),
                    domain.mobility_radius()),
                dt,
                domain.mobility_radius()));
        position_type const displacement(draw_displacement(domain, r));
        LOG_DEBUG(("draw_new_position(domain=%s, dt=%.16g): mobility_radius=%.16g, r=%.16g, displacement=%s (%.16g)",
                boost::lexical_cast<std::string>(domain).c_str(), dt,
                domain.mobility_radius(),
                r, boost::lexical_cast<std::string>(displacement).c_str(),
                length(displacement)));
        if (base_type::paranoiac_)
        {
            BOOST_ASSERT(r <= domain.mobility_radius());
            // length_type const scale(domain.particle().second.radius());
            // BOOST_ASSERT(feq(length(displacement), std::abs(r), scale));
        }
        return (*base_type::world_).apply_boundary(add(domain.particle().second.position(), displacement));
    }

    position_type draw_new_position(single_type& domain, time_type dt)
    {
        {
            spherical_single_type* _domain(dynamic_cast<spherical_single_type*>(&domain));
            if (_domain)
            {
                return draw_new_position(*_domain, dt);
            }
        }
        {
            cylindrical_single_type* _domain(dynamic_cast<cylindrical_single_type*>(&domain));
            if (_domain)
            {
                return draw_new_position(*_domain, dt);
            }
        }
        throw not_implemented(std::string("unsupported domain type"));
    }
    // }}}

    template<typename Tshell>
    position_type draw_escape_position(
            AnalyticalSingle<traits_type, Tshell> const& domain)
    {
        position_type const displacement(draw_displacement(domain, domain.mobility_radius()));
        LOG_DEBUG(("draw_escape_position(domain=%s): mobility_radius=%.16g, displacement=%s (%.16g)",
                boost::lexical_cast<std::string>(domain).c_str(),
                domain.mobility_radius(),
                boost::lexical_cast<std::string>(displacement).c_str(),
                length(displacement)));
        if (base_type::paranoiac_)
        {
            ; // do nothing
            // BOOST_ASSERT(feq(length(displacement)) <= domain.mobility_radius());
            // length_type const scale(domain.particle().second.radius());
            // BOOST_ASSERT(feq(length(displacement), std::abs(domain.mobility_radius()), scale));
        }
        return (*base_type::world_).apply_boundary(add(domain.particle().second.position(), displacement));
    }

    position_type draw_escape_position(single_type& domain)
    {
        {
            spherical_single_type* _domain(dynamic_cast<spherical_single_type*>(&domain));
            if (_domain)
            {
                return draw_escape_position(*_domain);
            }
        }
        {
            cylindrical_single_type* _domain(dynamic_cast<cylindrical_single_type*>(&domain));
            if (_domain)
            {
                return draw_escape_position(*_domain);
            }
        }
        throw not_implemented(std::string("unsupported domain type"));
    }

    // draw_new_positions {{{
    template<typename Tdraw, typename T>
    boost::array<position_type, 2> draw_new_positions(
        AnalyticalPair<traits_type, T> const& domain, time_type dt)
    {
        Tdraw d(this->rng(), *base_type::world_);
        position_type const new_com(d.draw_com(domain, dt));
        position_type const new_iv(d.draw_iv(domain, dt, domain.iv()));
        D_type const D0(domain.particles()[0].second.D());
        D_type const D1(domain.particles()[1].second.D());
        return array_gen(
            (*base_type::world_).apply_boundary(subtract(new_com, multiply(new_iv, D0 / (D0 + D1)))),
            (*base_type::world_).apply_boundary(add(new_com, multiply(new_iv, D1 / (D0 + D1)))));
    }
    // }}}

    // propagate {{{
    /**
     * The difference between a burst and a propagate is that a burst 
     * always takes place before the actual scheduled event for the single, 
     * while propagate_single can be called for an escape event.
     *
     * Another subtle difference is that burst_single always reschedules 
     * (update_event) the single, while just calling propagate does not. 
     * So whoever calls propagate_single directly should reschedule the single 
     * afterwards.
     */
    //template<typename T>
    //void propagate(AnalyticalSingle<traits_type, T>& domain, position_type const& new_pos)
    void propagate(single_type& domain, position_type const& new_pos,
                   bool do_update_shell_matrix)
    {
        LOG_DEBUG(("propagate: domain=%s, new_pos=%s, do_update_shell_matrix=%d",
                boost::lexical_cast<std::string>(domain).c_str(),
                boost::lexical_cast<std::string>(new_pos).c_str(),
                do_update_shell_matrix));

        if (base_type::paranoiac_)
        {
            particle_shape_type const new_particle(new_pos, domain.particle().second.radius());
            BOOST_ASSERT(check_overlap(new_particle, domain.particle().first));
        }

        particle_type const& old(domain.particle().second);
        domain.particle().second = particle_type(old.sid(),
                new_pos, old.radius(), old.D());
        (*base_type::world_).update_particle(
            domain.particle().first, domain.particle().second);

        domain.position() = new_pos;
        domain.size() = domain.particle().second.radius();
        if (do_update_shell_matrix)
            update_shell_matrix(domain);
    }

    template<typename T>
    boost::array<boost::shared_ptr<single_type>, 2>
    propagate(AnalyticalPair<traits_type, T>& domain,
              boost::array<position_type, 2> const& new_pos)
    {
        boost::array<particle_id_pair, 2> const& particles(domain.particles());
        boost::array<particle_id_pair, 2> new_particles(particles);
        new_particles[0].second.position() = new_pos[0];
        new_particles[1].second.position() = new_pos[1];

        if (base_type::paranoiac_)
        {
            BOOST_ASSERT(distance(domain, new_particles[0].second.position()) <= -new_particles[0].second.radius());
            BOOST_ASSERT(distance(domain, new_particles[1].second.position()) <= -new_particles[1].second.radius());
            BOOST_ASSERT(check_overlap(
                shape(new_particles[0].second),
                new_particles[0].first, new_particles[1].first));
            BOOST_ASSERT(check_overlap(
                shape(new_particles[1].second),
                new_particles[0].first, new_particles[1].first));
            BOOST_ASSERT(check_pair_pos(domain, new_particles));
        }

        (*base_type::world_).update_particle(
            new_particles[0].first, new_particles[0].second);
        (*base_type::world_).update_particle(
            new_particles[1].first, new_particles[1].second);

        remove_domain(domain);

        boost::array<boost::shared_ptr<single_type>, 2> const singles = { {
            create_single(new_particles[0]),
            create_single(new_particles[1]) 
        } };

        if (log_.level() == Logger::L_DEBUG)
        {
            for (int i = 0; i < 2; i++)
            {
                LOG_DEBUG(("propagate: #%d: %s => %s", i,
                    boost::lexical_cast<std::string>(particles[i].second.position()).c_str(),
                    boost::lexical_cast<std::string>(*singles[i]).c_str()));
            }
        }

        return singles;
    }
    // }}}


    template<typename Trange>
    void burst_domains(Trange const& domain_ids, boost::optional<std::vector<boost::shared_ptr<domain_type> >&> const& result = boost::optional<std::vector<boost::shared_ptr<domain_type> >&>())
    {
        BOOST_FOREACH(domain_id_type id, domain_ids)
        {
            boost::shared_ptr<domain_type> domain(get_domain(id));
            burst(domain, result);
        }
    }

    // burst {{{
    template<typename T>
    void burst(AnalyticalSingle<traits_type, T>& domain)
    {
        position_type const old_pos(domain.position());
        //length_type const old_shell_size(domain.size()); 
        length_type const particle_radius(domain.particle().second.radius());

        // Override dt, burst happens before single's scheduled event.
        domain.dt() = this->t() - domain.last_time();
        LOG_DEBUG(("t=%.16g, domain.last_time=%.16g", this->t(), domain.last_time()));

        position_type const new_pos(draw_new_position(domain, domain.dt()));

        propagate(domain, new_pos, true);

        domain.last_time() = this->t();
        domain.dt() = 0.;
        try
        {
            remove_event(domain);
            add_event(domain, SINGLE_EVENT_ESCAPE);
        }
        catch (std::out_of_range const&)
        {
            // event may have been removed.
            LOG_DEBUG(("event %s already removed; ignoring.", boost::lexical_cast<std::string>(domain.event().first).c_str()));
        }

        // Displacement check is in draw_new_position.
        // BOOST_ASSERT(
        //     (*base_type::world_).distance(new_pos, old_pos)
        //         <= old_shell_size - particle_radius);

        BOOST_ASSERT(domain.size() == particle_radius);
    }

    template<typename T>
    boost::array<boost::shared_ptr<single_type>, 2> burst(AnalyticalPair<traits_type, T>& domain)
    {
        length_type const dt(this->t() - domain.last_time());

        boost::array<boost::shared_ptr<single_type>, 2> const singles(
            propagate(domain, draw_new_positions<draw_on_burst>(domain, dt)));

        add_event(*singles[0], SINGLE_EVENT_ESCAPE);
        add_event(*singles[1], SINGLE_EVENT_ESCAPE);

        return singles;
    }

    void burst(multi_type& domain, boost::optional<std::vector<boost::shared_ptr<domain_type> >&> const& result = boost::optional<std::vector<boost::shared_ptr<domain_type> >&>())
    {
        BOOST_FOREACH(particle_id_pair p, domain.get_particles_range())
        {
            boost::shared_ptr<single_type> s(create_single(p));
            add_event(*s, SINGLE_EVENT_ESCAPE);
            if (result)
            {
                result.get().push_back(boost::dynamic_pointer_cast<domain_type>(s));
            }
        }
        remove_domain(domain);
    }

    void burst(single_type& domain)
    {
        LOG_DEBUG(("burst: bursting %s", boost::lexical_cast<std::string>(domain).c_str()));
        BOOST_ASSERT(this->t() >= domain.last_time());
        BOOST_ASSERT(this->t() <= domain.last_time() + domain.dt());
        {
            spherical_single_type* _domain(dynamic_cast<spherical_single_type*>(&domain));
            if (_domain)
            {
                burst(*_domain);
                return;
            }
        }
        {
            cylindrical_single_type* _domain(dynamic_cast<cylindrical_single_type*>(&domain));
            if (_domain)
            {
                burst(*_domain);
                return;
            }
        }
        throw not_implemented("?");
    }

    void burst(boost::shared_ptr<domain_type> domain, boost::optional<std::vector<boost::shared_ptr<domain_type> >&> const& result = boost::optional<std::vector<boost::shared_ptr<domain_type> >&>())
    {
        LOG_DEBUG(("burst: bursting %s", boost::lexical_cast<std::string>(*domain).c_str()));
        {
            spherical_single_type* _domain(dynamic_cast<spherical_single_type*>(domain.get()));
            if (_domain)
            {
                burst(*_domain);
                if (result)
                    result.get().push_back(domain);
                return;
            }
        }
        {
            cylindrical_single_type* _domain(dynamic_cast<cylindrical_single_type*>(domain.get()));
            if (_domain)
            {
                burst(*_domain);
                if (result)
                    result.get().push_back(domain);
                return;
            }
        }
        {
            spherical_pair_type* _domain(dynamic_cast<spherical_pair_type*>(domain.get()));
            if (_domain)
            {
                boost::array<boost::shared_ptr<single_type>, 2> bursted(burst(*_domain));
                if (result)
                {
                    result.get().push_back(boost::dynamic_pointer_cast<domain_type>(bursted[0]));
                    result.get().push_back(boost::dynamic_pointer_cast<domain_type>(bursted[1]));
                }
                return;
            }
        }
        {
            cylindrical_pair_type* _domain(dynamic_cast<cylindrical_pair_type*>(domain.get()));
            if (_domain)
            {
                boost::array<boost::shared_ptr<single_type>, 2> bursted(burst(*_domain));
                if (result)
                {
                    result.get().push_back(boost::dynamic_pointer_cast<domain_type>(bursted[0]));
                    result.get().push_back(boost::dynamic_pointer_cast<domain_type>(bursted[1]));
                }
                return;
            }
        }
        {
            multi_type* _domain(dynamic_cast<multi_type*>(domain.get()));
            if (_domain)
            {
                burst(*_domain, result);
                return;
            }
        }
        throw not_implemented("?");
    }
    // }}}

    // attempt_single_reaction {{{
    bool attempt_single_reaction(single_type& domain)
    {
        const particle_id_pair reactant(domain.particle());
        const molecule_info_type reactant_species((*base_type::world_).get_molecule_info(reactant.second.species()));
        // const molecule_info_type reactant_species((*base_type::world_).find_molecule_info(reactant.second.species()));
        reaction_rules const& rules((*base_type::network_rules_).query_reaction_rule(reactant.second.species()));
        if (::size(rules) == 0)
        {
            return false;
        }

        reaction_rule_type const& r(draw_reaction_rule(rules));
        LOG_DEBUG(("attempt_single_reaction: reactant=%s, products=[%s]",
                boost::lexical_cast<std::string>(reactant.second.sid()).c_str(),
                stringize_and_join(r.get_products(), ", ").c_str()));

        switch (::size(r.get_products()))
        {
        case 0:
            remove_domain(domain);
            (*base_type::world_).remove_particle(reactant.first);
            if (base_type::rrec_)
            {
                // (*base_type::rrec_)(reaction_record_type(
                //     r.id(), array_gen<particle_id_type>(), reactant.first));
                (*base_type::rrec_)(reaction_record_type(
                    r.id(), array_gen<particle_id_pair>(), reactant));
            }
            break;
        case 1: 
            {
                species_id_type const& product_id0(r.get_products()[0]);
                molecule_info_type const product_species(
                    (*base_type::world_).get_molecule_info(product_id0));

                if (reactant_species.radius < product_species.radius)
                    clear_volume(::shape(reactant.second), domain.id());

                if (!(*base_type::world_).no_overlap(
                    ::shape(reactant.second), reactant.first))
                {
                    LOG_INFO(("no space for product particle."));
                    throw no_space();
                }

                remove_domain(domain);
                (*base_type::world_).remove_particle(reactant.first);
                // particle_id_pair product(
                //     (*base_type::world_).new_particle(
                //         product_species.id(), reactant.second.position()).first);
                particle_id_pair product(
                    (*base_type::world_).new_particle(
                        product_id0, reactant.second.position()).first);
                boost::shared_ptr<single_type> new_domain(create_single(product));
                add_event(*new_domain, SINGLE_EVENT_ESCAPE);
                if (base_type::rrec_)
                {
                    // (*base_type::rrec_)(reaction_record_type(
                    //     r.id(), array_gen(product.first), reactant.first));
                    (*base_type::rrec_)(reaction_record_type(
                        r.id(), array_gen(product), reactant));
                }
            }
            break;
        case 2:
            {
                species_id_type const& product_id0(r.get_products()[0]);
                species_id_type const& product_id1(r.get_products()[1]);
                // molecule_info_type const* const product_species[] = {
                //     &(*base_type::world_).get_molecule_info(product_id0),
                //     &(*base_type::world_).get_molecule_info(product_id1)
                // };

                // D_type const D0(product_species[0]->D), D1(product_species[1]->D);
                // length_type const radius0(product_species[0]->radius),
                //     radius1(product_species[1]->radius);
                molecule_info_type const product_species[] = {
                    (*base_type::world_).get_molecule_info(product_id0),
                    (*base_type::world_).get_molecule_info(product_id1)
                };

                D_type const D0(product_species[0].D), D1(product_species[1].D);
                length_type const radius0(product_species[0].radius),
                    radius1(product_species[1].radius);
                D_type const D01(D0 + D1);
                length_type r01(radius0 + radius1);
                Real const rad(std::max(
                        r01 * (D0 / D01) + radius0, r01 * (D1 / D01) + radius1));
                clear_volume(particle_shape_type(reactant.second.position(), rad), domain.id());

                particle_shape_type new_particles[2];

                int i = num_retries_;
                while (--i >= 0)
                {
                    boost::shared_ptr<structure_type> structure(
                        (*base_type::world_).get_structure(
                            reactant_species.structure_id));
                    position_type vector(
                        structure->random_vector(
                            r01 * traits_type::MINIMAL_SEPARATION_FACTOR,
                            this->rng()));
                    // place particles according to the ratio D1:D2
                    // this way, species with D=0 doesn't move.
                    // FIXME: what if D1 == D2 == 0?
                    for (;;) {
                        new_particles[0] = particle_shape_type(
                            (*base_type::world_).apply_boundary(
                                add(reactant.second.position(),
                                    multiply(vector, D0 / D01))),
                            radius0);
                        new_particles[1] = particle_shape_type(
                            (*base_type::world_).apply_boundary(
                                add(reactant.second.position(),
                                    multiply(vector, -D1 / D01))),
                            radius1);

                        length_type const distance_between_new_particles(
                            (*base_type::world_).distance(
                                new_particles[0].position(),
                                new_particles[1].position()));
                        if (distance_between_new_particles >= r01)
                            break;

                        vector = multiply(vector, 1.0 + 1e-7);
                    }

                    // accept the new positions if there is enough space.
                    if (((*base_type::world_).no_overlap(
                            new_particles[0], reactant.first)) &&
                        ((*base_type::world_).no_overlap(
                            new_particles[1], reactant.first)))
                        break;
                }
                if (i < 0)
                {
                    LOG_INFO(("no space for product particles."));
                    throw no_space();
                }

                remove_domain(domain);
                (*base_type::world_).remove_particle(reactant.first);

                particle_id_pair const pp[] = {
                    (*base_type::world_).new_particle(
                        product_id0, new_particles[0].position()).first,
                    (*base_type::world_).new_particle(
                        product_id1, new_particles[1].position()).first
                };
                // create domains for two particles and add them to
                // the event queue
                add_event(*create_single(pp[0]), SINGLE_EVENT_ESCAPE);
                add_event(*create_single(pp[1]), SINGLE_EVENT_ESCAPE);

                if (base_type::rrec_)
                {
                    // (*base_type::rrec_)(reaction_record_type(
                    //     r.id(), array_gen(pp[0].first, pp[1].first),
                    //     reactant.first));
                    (*base_type::rrec_)(reaction_record_type(
                        r.id(), array_gen(pp[0], pp[1]), reactant));
                }
            }
            break;
        default:
            throw not_implemented("reactions that produces more than two products are not supported.");
        }
        return true;
    }
    // }}}

    time_type draw_single_reaction_time(species_id_type const& sid)
    {
        reaction_rules const& rules(
            (*base_type::network_rules_).query_reaction_rule(sid));
        rate_type const k_tot(calculate_k_tot(rules));
        if (k_tot <= 0.)
        {
            return std::numeric_limits<time_type>::infinity();
        }
        else if (k_tot == std::numeric_limits<rate_type>::infinity())
        {
            return 0.;
        }
        else
        {
            const double rnd(this->rng().uniform(0., 1.));
            if(rnd <= 0.)
            {
                return std::numeric_limits<time_type>::infinity();
            }
            else
            {
                return (1. / k_tot) * (- std::log(rnd)); // log(1/x) == - log(x)
            }
        }
    }

    template<typename Tshell>
    time_type draw_escape_or_interaction_time(AnalyticalSingle<traits_type, Tshell> const& domain)
    {
        if (domain.particle().second.D() == 0.)
        {
            return std::numeric_limits<time_type>::infinity();
        }
        else
        {
            typedef Tshell shell_type;
            typedef typename shell_type::shape_type shape_type;
            typedef typename detail::get_greens_function<shape_type>::type greens_function;
            return greens_function(domain.particle().second.D(),
                            domain.mobility_radius())
                .drawTime(this->rng().uniform(0., 1.));
        }
    }

    template<typename Tshell>
    std::pair<time_type, pair_event_kind>
    draw_com_escape_or_iv_event_time(AnalyticalPair<traits_type, Tshell> const& domain)
    {
        typedef Tshell shell_type;
        typedef typename shell_type::shape_type shape_type;
        typedef typename detail::get_pair_greens_function<shape_type> pair_greens_functions;
        typedef typename pair_greens_functions::iv_type iv_greens_function;
        typedef typename pair_greens_functions::com_type com_greens_function;
        BOOST_ASSERT(::size(domain.reactions()) == 1);
        time_type const dt_com(
            com_greens_function(domain.D_R(), domain.a_R()).drawTime(this->rng().uniform(0., 1.)));
        time_type const dt_iv(
            iv_greens_function(domain.D_tot(), domain.reactions()[0].k(),
                           domain.r0(), domain.sigma(), domain.a_r()).drawTime(this->rng().uniform(0., 1.)));
        if (dt_com < dt_iv)
        {
            return std::make_pair(dt_com, PAIR_EVENT_COM_ESCAPE);
        }
        else
        {
            return std::make_pair(dt_iv, PAIR_EVENT_IV_UNDETERMINED);
        }
    }

    template<typename Tshell>
    std::pair<time_type, pair_event_kind>
    draw_single_reaction_time(AnalyticalPair<traits_type, Tshell> const& domain)
    {
        time_type const dt[2] = {
            draw_single_reaction_time(domain.particles()[0].second.species()),
            draw_single_reaction_time(domain.particles()[1].second.species())
        };
        if (dt[0] < dt[1])
        {
            return std::make_pair(dt[0], PAIR_EVENT_SINGLE_REACTION_0);
        }
        else
        {
            return std::make_pair(dt[1], PAIR_EVENT_SINGLE_REACTION_1);
        }
    }

    // {{{ determine_next_event
    template<typename Tshell>
    void determine_next_event(AnalyticalSingle<traits_type, Tshell>& domain)
    {
        typedef Tshell shell_type;
        typedef typename shell_type::shape_type shape_type;
        typedef typename detail::get_greens_function<shape_type>::type greens_function;
        time_type const dt_reaction(draw_single_reaction_time(domain.particle().second.species()));
        time_type const dt_escape_or_interaction(draw_escape_or_interaction_time(domain));
        LOG_DEBUG(("determine_next_event: %s => dt_reaction=%.16g, "
                   "dt_escape_or_interaction=%.16g",
                   boost::lexical_cast<std::string>(domain).c_str(),
                   dt_reaction, dt_escape_or_interaction));
        single_event_kind event_kind;
        if (dt_reaction < dt_escape_or_interaction)
        {
            domain.dt() = dt_reaction;
            event_kind = SINGLE_EVENT_REACTION;
        }
        else
        {
            domain.dt() = dt_escape_or_interaction;
            event_kind = SINGLE_EVENT_ESCAPE;
        }

        domain.last_time() = this->t();
        add_event(domain, event_kind);
    }

    void determine_next_event(single_type& domain)
    {
        {
            spherical_single_type* _domain(
                dynamic_cast<spherical_single_type*>(&domain));
            if (_domain)
            {
                determine_next_event(*_domain);
                return;
            }
        }
        {
            cylindrical_single_type* _domain(
                dynamic_cast<cylindrical_single_type*>(&domain));
            if (_domain)
            {
                determine_next_event(*_domain);
                return;
            }
        }
        throw not_implemented("unsupported domain type");
    }

    template<typename Tshell>
    void determine_next_event(AnalyticalPair<traits_type, Tshell>& domain)
    {
        std::pair<time_type, pair_event_kind> const dt_reaction(draw_single_reaction_time(domain));
        std::pair<time_type, pair_event_kind> const dt_com_escape_or_iv_event(draw_com_escape_or_iv_event_time(domain));
        std::pair<time_type, pair_event_kind> const dt_and_event_pair(
            dt_reaction.first < dt_com_escape_or_iv_event.first ?
                dt_reaction: dt_com_escape_or_iv_event);
        domain.dt() = dt_and_event_pair.first;
        domain.last_time() = this->t();
        add_event(domain, dt_and_event_pair.second);
    }

    void determine_next_event(pair_type& domain)
    {
        {
            spherical_pair_type* _domain(
                dynamic_cast<spherical_pair_type*>(&domain));
            if (_domain)
            {
                determine_next_event(*_domain);
                return;
            }
        }
        {
            cylindrical_pair_type* _domain(
                dynamic_cast<cylindrical_pair_type*>(&domain));
            if (_domain)
            {
                determine_next_event(*_domain);
                return;
            }
        }
        throw not_implemented("unsupported domain type");
    }
    // }}}

    // get_intruders {{{ 
    std::pair<std::vector<domain_id_type>*,
              std::pair<domain_id_type, length_type> >
    get_intruders(particle_shape_type const& p,
                  domain_id_type const& ignore) const
    {
        typedef intruder_collector collector_type;

        collector_type col((*base_type::world_), p, ignore);
        boost::fusion::for_each(smatm_, shell_collector_applier<collector_type>(col, p.position()));
        return std::make_pair(col.intruders.container().get(), col.closest);
    }
    // }}}

    template<typename TdidSet>
    std::pair<domain_id_type, length_type>
    get_closest_domain(position_type const& p, TdidSet const& ignore) const
    {
        typedef closest_object_finder<TdidSet> collector_type;

        collector_type col((*base_type::world_), p, ignore);
        boost::fusion::for_each(smatm_, shell_collector_applier<collector_type>(col, p));
        return col.closest;
    }

    void restore_domain(single_type& domain)
    {
        std::pair<domain_id_type, length_type> const closest(
            get_closest_domain(
                domain.position(), 
                array_gen(domain.id())));
        restore_domain(domain, closest);
    }

    template<typename T>
    void restore_domain(AnalyticalSingle<traits_type, T>& domain,
                        std::pair<domain_id_type, length_type> const& closest)
    {
        typedef typename AnalyticalSingle<traits_type, T>::shell_type shell_type;
        domain_type const* closest_domain(
            closest.second == std::numeric_limits<length_type>::infinity() ?
                (domain_type const*)0: get_domain(closest.first).get());
        length_type new_shell_size(0.);

        if (closest_domain)
        {
            single_type const* const _closest_domain(
                dynamic_cast<single_type const*>(closest_domain));
            if (_closest_domain)
            {
                length_type const distance_to_closest(
                    (*base_type::world_).distance(
                        domain.position(), _closest_domain->position()));
                new_shell_size = calculate_single_shell_size(
                        domain, *_closest_domain,
                        distance_to_closest,
                        closest.second);
            } else {
                new_shell_size = closest.second / traits_type::SAFETY;
            }
            new_shell_size = std::min(max_shell_size(), 
                std::max(domain.particle().second.radius(), new_shell_size));
        }
        else
        {
            new_shell_size = max_shell_size();
        }
        LOG_DEBUG(("restore domain: %s (shell_size=%.16g, dt=%.16g) closest=%s (distance=%.16g)",
            boost::lexical_cast<std::string>(domain).c_str(),
            new_shell_size,
            domain.dt(),
            closest_domain ?
                boost::lexical_cast<std::string>(*closest_domain).c_str():
                "(none)",
            closest.second));
        if (base_type::paranoiac_)
        {
            BOOST_ASSERT(check_overlap(
                particle_shape_type(domain.position(), new_shell_size),
                domain.particle().first));
        }
        domain.size() = new_shell_size;
        update_shell_matrix(domain);
        BOOST_ASSERT(domain.size() == new_shell_size);
    }

    void restore_domain(single_type& domain,
                        std::pair<domain_id_type, length_type> const& closest)
    {
        {
            spherical_single_type *_domain(
                dynamic_cast<spherical_single_type*>(&domain));
            if (_domain)
                return restore_domain(*_domain, closest);
        }
        {
            cylindrical_single_type *_domain(
                dynamic_cast<cylindrical_single_type*>(&domain));
            if (_domain)
                return restore_domain(*_domain, closest);
        }
        throw not_implemented(std::string("unsupported domain type"));
    }

    template<typename Trange>
    void burst_non_multis(Trange const& domain_ids,
                          std::vector<boost::shared_ptr<domain_type> >& bursted)
    {
        BOOST_FOREACH (domain_id_type id, domain_ids)
        {
            boost::shared_ptr<domain_type> domain(get_domain(id));
            if (dynamic_cast<multi_type*>(domain.get()))
            {
                bursted.push_back(domain);
            }
            else
            {
                burst(domain, bursted);
            }
        }
    }

    template<typename T>
    length_type distance(AnalyticalSingle<traits_type, T> const& domain,
                         position_type const& pos) const
    {
        return (*base_type::world_).distance(shape(domain.shell().second), pos);
    }

    template<typename T>
    length_type distance(AnalyticalPair<traits_type, T> const& domain,
                         position_type const& pos) const
    {
        return (*base_type::world_).distance(shape(domain.shell().second), pos);
    }

    length_type distance(multi_type const& domain, position_type const& pos) const
    {
        length_type retval(std::numeric_limits<length_type>::infinity());
        BOOST_FOREACH (spherical_shell_id_pair const& shell,
                       domain.get_shells())
        {
            length_type const dist((*base_type::world_).distance(
                    shape(shell.second), pos));
            if (retval > dist)
            {
                retval = dist;
            }
        }
        return retval;
    }

    length_type distance(domain_type const& domain, position_type const& pos) const
    {
        length_type retval;

        struct distance_visitor: ImmutativeDomainVisitor<traits_type>
        {
            virtual ~distance_visitor() {}

            virtual void operator()(multi_type const& domain) const
            {
                retval = outer.distance(domain, pos);
            }

            virtual void operator()(spherical_single_type const& domain) const
            {
                retval = outer.distance(domain, pos);
            }

            virtual void operator()(cylindrical_single_type const& domain) const
            {
                retval = outer.distance(domain, pos);
            }

            virtual void operator()(spherical_pair_type const& domain) const
            {
                retval = outer.distance(domain, pos);
            }

            virtual void operator()(cylindrical_pair_type const& domain) const
            {
                retval = outer.distance(domain, pos);
            }

            distance_visitor(EGFRDSimulator const& outer, position_type const& pos,
                             length_type& retval)
                : outer(outer), pos(pos), retval(retval) {}

            EGFRDSimulator const& outer;
            position_type const& pos;
            length_type& retval;
        };

        domain.accept(distance_visitor(*this, pos, retval));
        return retval;
    }

    boost::optional<pair_type&>
    form_pair(single_type& domain, single_type& possible_partner,
              std::vector<boost::shared_ptr<domain_type> > const& neighbors)
    {
        LOG_DEBUG(("trying to form Pair(%s, %s)",
                    boost::lexical_cast<std::string>(domain).c_str(),
                    boost::lexical_cast<std::string>(possible_partner).c_str()));
        // 1. Determine min shell size.
        length_type const r[] = {
            domain.particle().second.radius(),
            possible_partner.particle().second.radius()
        };
        length_type const sigma(r[0] + r[1]);

        D_type const D[] = {
            domain.particle().second.D(),
            possible_partner.particle().second.D()
        };
        D_type const D01(D[0] + D[1]);

        BOOST_ASSERT(domain.particle().second.position() ==
                     domain.position());

        BOOST_ASSERT(possible_partner.particle().second.position() ==
                     possible_partner.position());

        position_type iv(
                subtract(domain.position(),
                    (*base_type::world_).periodic_transpose(
                        possible_partner.position(),
                        domain.position())));
        length_type const r0(length(iv));
        length_type const distance_from_sigma(r0 - sigma);
        BOOST_ASSERT(distance_from_sigma >= 0);

        length_type const shell_size[] = {
           r0 * D[0] / D01 + r[0], r0 * D[1] / D01 + r[1]
        };
        length_type const shell_size_margin[] = {
            r[0] * 2,
            r[1] * 2
        };
        size_t const larger_shell_index(
            shell_size[0] + shell_size_margin[0] >=
                shell_size[1] + shell_size_margin[1] ? 0: 1);
        length_type const min_shell_size(shell_size[larger_shell_index]);
        length_type const min_shell_size_margin(shell_size_margin[larger_shell_index]);

        // 2. Check if min shell size not larger than max shell size or 
        // sim cell size.
        position_type com((*base_type::world_).apply_boundary(
            (*base_type::world_).calculate_pair_CoM(
                domain.position(), possible_partner.position(),
                D[0], D[1])));
        length_type const min_shell_size_with_margin(
            min_shell_size + min_shell_size_margin);
        length_type const max_shell_size(
            std::min(this->max_shell_size(),
                     distance_from_sigma * 100
                        + sigma + min_shell_size_margin));

        if (min_shell_size_with_margin >= max_shell_size)
        {
            LOG_DEBUG(("Pair(%s, %s) not formed: min_shell_size %.16g >="
                       "max_shell_size %.16g",
                       boost::lexical_cast<std::string>(domain).c_str(),
                       boost::lexical_cast<std::string>(possible_partner).c_str(),
                       min_shell_size_with_margin, max_shell_size));
            return boost::optional<pair_type&>();
        }

        // 3. Check if bursted Singles not too close.
        // The simple check for closest below could miss
        // some of them, because sizes of these Singles for this
        // distance check has to include SINGLE_SHELL_FACTOR, while
        // these burst objects have zero mobility radii.  This is not
        // beautiful, a cleaner framework may be possible.

        domain_type* closest_domain (0);
        length_type closest_shell_distance(std::numeric_limits<length_type>::infinity());
        BOOST_FOREACH (boost::shared_ptr<domain_type> _neighbor, neighbors)
        {
            single_type* const neighbor(
                dynamic_cast<single_type*>(_neighbor.get()));
            if (neighbor && neighbor->id() != possible_partner.id())
            {
                length_type const shell_distance(
                    (*base_type::world_).distance(com, neighbor->position()) -
                        neighbor->particle().second.radius() *
                            (1.0 + traits_type::SINGLE_SHELL_FACTOR));
                if (shell_distance < closest_shell_distance)
                {
                    closest_domain = neighbor;
                    closest_shell_distance = shell_distance;
                }
            }
        }

        if (closest_domain)
        {
            BOOST_ASSERT(closest_shell_distance > 0);

            if (closest_shell_distance <= min_shell_size_with_margin)
            {
                LOG_DEBUG(("Pair(%s, %s) not formed: squeezed by burst neighbor %s",
                           boost::lexical_cast<std::string>(domain).c_str(),
                           boost::lexical_cast<std::string>(possible_partner).c_str(),
                           boost::lexical_cast<std::string>(*closest_domain).c_str()));
                return boost::optional<pair_type&>();
            }
        }

        // 4. Determine shell size and check if closest object not too 
        // close (squeezing).
        {
            std::pair<domain_id_type, length_type> possible_closest(
                get_closest_domain(com, array_gen(domain.id(),
                                                  possible_partner.id())));
            if (possible_closest.second < closest_shell_distance)
            {
                domain_type* const _closest_domain(
                        get_domain(possible_closest.first).get());
                closest_domain = _closest_domain;
                closest_shell_distance = possible_closest.second;
            }
        }

        if (closest_domain)
        {
            LOG_DEBUG(("Pair closest neighbor: %s %.16g, "
                       "min_shell_with_margin=%.16g",
                       boost::lexical_cast<std::string>(*closest_domain).c_str(),
                       closest_shell_distance,
                       min_shell_size_with_margin));
            BOOST_ASSERT(closest_shell_distance > 0);
        }

        length_type new_shell_size(0.);

        {
            single_type* const _closest_domain(
                    dynamic_cast<single_type*>(closest_domain));
            if (_closest_domain)
            {
                particle_type const& closest_domain_particle(
                        _closest_domain->particle().second);
                D_type const D_tot(closest_domain_particle.D() + D01);
                length_type const closest_particle_distance(
                    (*base_type::world_).distance(
                        com, closest_domain_particle.position()));
                length_type const closest_min_shell(
                        closest_domain_particle.radius() *
                            (traits_type::SINGLE_SHELL_FACTOR + 1.0));

                // options for shell size:
                // a. ideal shell size
                // b. closest shell is from a bursted single
                // c. closest shell is closer than ideal shell size 
                new_shell_size = std::min(
                    std::min(
                        (D01 / D_tot) * (
                            closest_particle_distance - min_shell_size 
                            - closest_domain_particle.radius())
                        + min_shell_size,
                        closest_particle_distance - closest_min_shell),
                    closest_shell_distance);
            }
            else
            {
                new_shell_size = closest_shell_distance;
            }
            new_shell_size /= traits_type::SAFETY;

            if (new_shell_size <= min_shell_size_with_margin)
            {
                LOG_DEBUG(("Pair(%s, %s) not formed%s%s",
                    boost::lexical_cast<std::string>(domain).c_str(),
                    boost::lexical_cast<std::string>(possible_partner).c_str(),
                    closest_domain ? ": squeezed by ": "",
                    closest_domain ? boost::lexical_cast<std::string>(*closest_domain).c_str(): ""));
                return boost::optional<pair_type&>();
            }
        }

        // 5. Check if singles would not be better.
        {
            length_type const dist[] = {
                (*base_type::world_).distance(com, domain.position()),
                (*base_type::world_).distance(com, domain.position())
            };

            if (new_shell_size < std::max(
                    dist[0] + r[0] *
                        (1.0 + traits_type::SINGLE_SHELL_FACTOR),
                    dist[1] + r[1] *
                        (1.0 + traits_type::SINGLE_SHELL_FACTOR)) * 1.3)
            {
                LOG_DEBUG(("Pair(%s, %s) not formed: leaving singles are better",
                            boost::lexical_cast<std::string>(domain).c_str(),
                            boost::lexical_cast<std::string>(possible_partner).c_str()));
                return boost::optional<pair_type&>();
            }
        }

        // 6. Ok, Pair makes sense. Create one.
        new_shell_size = std::min(new_shell_size, max_shell_size);

        boost::shared_ptr<pair_type> new_pair(
            create_pair(
                domain.particle(),
                possible_partner.particle(),
                com, iv, new_shell_size));

        determine_next_event(*new_pair);
        BOOST_ASSERT(new_pair->dt() >= 0);

        new_pair->last_time() = this->t();

        remove_domain(domain);
        remove_domain(possible_partner);

        BOOST_ASSERT(
                (closest_domain && closest_shell_distance ==
                    std::numeric_limits<length_type>::infinity())
                || new_shell_size < closest_shell_distance);
        BOOST_ASSERT(new_shell_size >= min_shell_size_with_margin);
        BOOST_ASSERT(new_shell_size <= max_shell_size);

        LOG_INFO(("new_pair=%s, closest_shell_distance=%.16g, closest=%s",
                  boost::lexical_cast<std::string>(*new_pair).c_str(),
                  closest_shell_distance,
                  closest_domain ? boost::lexical_cast<std::string>(closest_domain).c_str(): "(none)"));

        return *new_pair;
    }

    boost::optional<multi_type&>
    form_multi(single_type& domain,
               std::vector<boost::shared_ptr<domain_type> > const& neighbors,
               std::pair<domain_type*, length_type> closest)
    {
        LOG_DEBUG(("form multi: neighbors=[%s], closest=%s",
                stringize_and_join(
                    make_transform_iterator_range(neighbors,
                        dereference<boost::shared_ptr<domain_type> >()),
                    ", ").c_str(),
                boost::lexical_cast<std::string>(*closest.first).c_str()));
        length_type const min_shell_size(
                domain.particle().second.radius() *
                    (1.0 + multi_shell_factor_));

        // Multis shells need to be contiguous.
        if (closest.second > min_shell_size)
        {
            LOG_DEBUG(("multi shells aren't close enough to each other (closest distance=%.16g, min_shell_size=%.16g)", closest.second, min_shell_size));
            return boost::optional<multi_type&>();
        }

        // If there's a multi neighbor, merge others into it.
        // Otherwise, create a new multi and let it hold them all.
        multi_type* retval(0);
        retval = dynamic_cast<multi_type*>(closest.first);
        if (!retval)
        {
            retval = create_multi().get();
            add_event(*retval);
            LOG_DEBUG(("form multi: created a new multi %s",
                    boost::lexical_cast<std::string>(*retval).c_str()));
        }

        position_type const single_pos(domain.position());
        add_to_multi(*retval, domain);

        BOOST_FOREACH (boost::shared_ptr<domain_type> neighbor, neighbors)
        {
            length_type const dist(distance(*neighbor, single_pos));
            if (dist < min_shell_size)
                add_to_multi_recursive(*retval, *neighbor); 
        }

        return *retval;
    }

    bool add_to_multi(multi_type& multi, single_type& single)
    {
        LOG_DEBUG(("adding single to multi: %s => %s",
                boost::lexical_cast<std::string>(single).c_str(),
                boost::lexical_cast<std::string>(multi).c_str()));

        if (!multi.add_particle(single.particle()))
        {
            LOG_DEBUG(("particle %s is already added to %s",
                boost::lexical_cast<std::string>(single.particle().first).c_str(),
                boost::lexical_cast<std::string>(multi).c_str()));
            return false;
        }

        spherical_shell_id_pair sid_shell_pair(
            new_shell(
                multi.id(),
                typename spherical_shell_type::shape_type(
                    single.particle().second.position(),
                    single.particle().second.radius() *
                        (1. + multi_shell_factor_))));
        multi.add_shell(sid_shell_pair);
        remove_domain(single);
        return true;
    }

    void add_to_multi(multi_type& multi, multi_type& other_multi)
    {
        if (multi.id() == other_multi.id())
        {
            LOG_DEBUG(("add multi to multi: given two multis are the same; do nothing"));
            return;
        }
        if (multi.has_particle(other_multi.get_particles_range().front().first))
        {
            LOG_DEBUG(("add multi to multi: given multi already added."));
            return;
        }

        LOG_DEBUG(("adding multi to multi: %s => %s",
                boost::lexical_cast<std::string>(other_multi).c_str(),
                boost::lexical_cast<std::string>(multi).c_str()));

        // merge other_multi into multi. other_multi will be removed.
        spherical_shell_matrix_type& mat(
            *boost::fusion::at_key<spherical_shell_type>(smatm_));
        BOOST_FOREACH (spherical_shell_id_pair const& _shell,
                       other_multi.get_shells())
        {
            typename spherical_shell_matrix_type::iterator const i(
                mat.find(_shell.first));
            BOOST_ASSERT(i != mat.end());
            spherical_shell_type& shell((*i).second);
            shell.did() = multi.id();
            multi.add_shell(spherical_shell_id_pair(_shell.first, shell));
        }

        BOOST_FOREACH (particle_id_pair const& particle,
                       other_multi.get_particles_range())
        {
            multi.add_particle(particle);
        }

        remove_domain_but_shell(other_multi);
    }

    void add_to_multi_recursive(multi_type& multi, domain_type& domain)
    {
        LOG_DEBUG(("add_to_multi_recursive: multi=%s, domain=%s",
                boost::lexical_cast<std::string>(multi).c_str(),
                boost::lexical_cast<std::string>(domain).c_str()));
        {
            single_type* single(dynamic_cast<single_type*>(&domain));
            if (single)
            {
                particle_shape_type const new_shell(
                    single->particle().second.position(),
                    single->particle().second.radius() *
                        (1.0 + multi_shell_factor_));

                if (!add_to_multi(multi, *single))
                {
                    return;
                }

                boost::scoped_ptr<std::vector<domain_id_type> > neighbors(
                    get_neighbor_domains(new_shell, single->id()));

                std::vector<boost::shared_ptr<domain_type> > bursted;
                burst_non_multis(*neighbors, bursted);

                LOG_DEBUG(("add_to_multi_recursive: bursted=[%s]",
                        stringize_and_join(
                            make_transform_iterator_range(
                                bursted,
                                dereference<boost::shared_ptr<domain_type> >()),
                            ", ").c_str()));

                BOOST_FOREACH (boost::shared_ptr<domain_type> neighbor, bursted)
                {
                    length_type const dist(distance(*neighbor, single->position()));
                    if (dist < new_shell.radius())
                        add_to_multi_recursive(multi, *neighbor); 
                }
                return;
            }
        }
        {
            multi_type* other_multi(dynamic_cast<multi_type*>(&domain));
            if (other_multi)
            {
                add_to_multi(multi, *other_multi);
            }
        }
    }

    boost::optional<domain_type&> form_pair_or_multi(
        single_type& domain,
        std::vector<boost::shared_ptr<domain_type> > const& neighbors)
    {
        BOOST_ASSERT(!neighbors.empty());

        domain_type* possible_partner(0);
        length_type length_to_possible_partner(
                std::numeric_limits<length_type>::infinity());
        BOOST_FOREACH (boost::shared_ptr<domain_type> neighbor, neighbors)
        {
            length_type const dist(distance(*neighbor, domain.position()));
            if (dist < length_to_possible_partner)
            {
                possible_partner = neighbor.get();
                length_to_possible_partner = dist;
            }
        }

        // First, try forming a Pair.
        {
            single_type* const _possible_partner(
                    dynamic_cast<single_type*>(possible_partner));
            if (_possible_partner)
            {
                boost::optional<pair_type&> new_pair(
                    form_pair(domain, *_possible_partner, neighbors));
                if (new_pair)
                {
                    return new_pair.get();
                }
            }
        }

        // If a Pair is not formed, then try forming a Multi.
        {
            boost::optional<multi_type&> new_multi(
                    form_multi(domain, neighbors,
                               std::pair<domain_type*, length_type>(
                                    possible_partner,
                                    length_to_possible_partner)));
            if (new_multi)
            {
                return new_multi.get();
            }
        }
        return boost::optional<domain_type&>();
    }
  
    void fire_event(single_event const& event)
    {
        single_type& domain(event.domain());
#if 0
        BOOST_ASSERT(
            std::abs(domain.dt() + domain.last_time() - this->t())
                <= 1e-18 * this->t());
#endif
        ++single_step_count_[event.kind()];
        switch (event.kind())
        {
        default: /* never get here */ BOOST_ASSERT(0); break;
        case SINGLE_EVENT_REACTION:
            LOG_DEBUG(("fire_single: single reaction (%s)", boost::lexical_cast<std::string>(domain).c_str()));
            propagate(domain, draw_new_position(domain, domain.dt()), false);
            try
            {
                attempt_single_reaction(domain);
            }
            catch (no_space const&)
            {
                LOG_DEBUG(("single reaction rejected"));
                ++rejected_moves_;
                domain.dt() = 0.;
                domain.last_time() = this->t();
                add_event(domain, SINGLE_EVENT_ESCAPE);
            }
            break;

        case SINGLE_EVENT_ESCAPE:
            LOG_DEBUG(("fire_single: single escape (%s)", boost::lexical_cast<std::string>(domain).c_str()));

            // handle immobile case
            if (domain.D() == 0.)
            {
                determine_next_event(domain);
                domain.last_time() = this->t();
                return;
            }

            if (domain.dt() != 0.)
                // Heads up: shell matrix will be updated later in restore_domain().
                // propagate(domain, draw_new_position(domain, domain.dt()), false);
                propagate(domain, draw_escape_position(domain), false);
            length_type const min_shell_radius(domain.particle().second.radius() * (1. + single_shell_factor_));
            {
                std::vector<domain_id_type>* intruders;
                std::pair<domain_id_type, length_type> closest;

                // boost::tie(intruders, closest) = get_intruders(
                //     particle_shape_type(
                //         domain.position(), min_shell_radius), domain.id());
                {
                    std::pair<std::vector<domain_id_type>*, 
                        std::pair<domain_id_type, length_type> > 
                        res(get_intruders(particle_shape_type(
                                              domain.position(), 
                                              min_shell_radius), 
                                          domain.id()));
                    intruders = res.first;
                    closest = res.second;
                }

                boost::scoped_ptr<std::vector<domain_id_type> > _(intruders);

                LOG_DEBUG(("intruders: %s, closest: %s (dist=%.16g)",
                    intruders ?
                        stringize_and_join(*intruders, ", ").c_str():
                        "(none)",
                    boost::lexical_cast<std::string>(closest.first).c_str(),
                    closest.second));
                if (intruders)
                {
                    std::vector<boost::shared_ptr<domain_type> > bursted;
                    burst_non_multis(*intruders, bursted);
                    if (form_pair_or_multi(domain, bursted))
                        return;
                    // if nothing was formed, recheck closest and restore shells.
                    restore_domain(domain);
                    BOOST_FOREACH (boost::shared_ptr<domain_type> _single, bursted)
                    {
                        boost::shared_ptr<single_type> single(
                            boost::dynamic_pointer_cast<single_type>(_single));
                        if (!single)
                            continue;
                        restore_domain(*single);
                        // reschedule events for the restored domains
                        remove_event(*single);
                        determine_next_event(*single);
                    }
                } else {
                    restore_domain(domain, closest);
                }
                determine_next_event(domain);
                LOG_DEBUG(("%s (dt=%.16g)",
                    boost::lexical_cast<std::string>(domain).c_str(),
                    domain.dt()));
            }
        }
    }

    template<typename Tshell>
    greens_functions::GreensFunction3DRadAbs::EventKind
    draw_iv_event_type(AnalyticalPair<traits_type, Tshell> const& domain)
    {
        typedef Tshell shell_type;
        typedef typename shell_type::shape_type shape_type;
        typedef typename detail::get_pair_greens_function<shape_type>::iv_type iv_greens_function;
        // Draw actual pair event for iv at very last minute.
        BOOST_ASSERT(::size(domain.reactions()) == 1);
        reaction_rule_type const& r(domain.reactions()[0]);
        iv_greens_function const gf(domain.D_tot(), r.k(), domain.r0(), domain.sigma(), domain.a_r());

        double const rnd(this->rng().uniform(0, 1.));
        return gf.drawEventType(rnd, domain.dt());
    }

    void fire_event(pair_event const& event)
    {
        {
            spherical_pair_type* _domain(dynamic_cast<spherical_pair_type*>(&event.domain()));
            if (_domain)
            {
                fire_event(*_domain, event.kind());
                return;
            }
        }
        {
            cylindrical_pair_type* _domain(dynamic_cast<cylindrical_pair_type*>(&event.domain()));
            if (_domain)
            {
                fire_event(*_domain, event.kind());
                return;
            }
        }
    }

    template<typename T>
    void fire_event(AnalyticalPair<traits_type, T>& domain, pair_event_kind kind)
    {
        typedef AnalyticalSingle<traits_type, T> corresponding_single_type;

        if (kind == PAIR_EVENT_IV_UNDETERMINED)
        {
            // Draw actual pair event for iv at very last minute.
            switch (draw_iv_event_type(domain))
            {
            case greens_functions::GreensFunction3DRadAbs::IV_ESCAPE:
                kind = PAIR_EVENT_IV_ESCAPE;
                break;
            case greens_functions::GreensFunction3DRadAbs::IV_REACTION:
                kind = PAIR_EVENT_IV_REACTION;
                break;
            }
        }

        ++pair_step_count_[kind];
        LOG_DEBUG(("fire_pair: %s", stringize_event_kind(kind).c_str()));

        //  1. Single reaction
        //  2. Pair reaction
        //  3a. IV escape
        //  3b. com escape

        switch (kind)
        {
        default: /* never get here */ BOOST_ASSERT(0); break;
        case PAIR_EVENT_SINGLE_REACTION_0: 
        case PAIR_EVENT_SINGLE_REACTION_1:
            {
                int const index(kind == PAIR_EVENT_SINGLE_REACTION_0 ? 0 : 1);
                // TODO.
                //int const theother_index(1 - index);
                position_type const old_CoM(domain.position());
                LOG_DEBUG(("pair: single reaction %s", boost::lexical_cast<std::string>(domain.particles()[index].first).c_str()));

                boost::array<boost::shared_ptr<single_type>, 2> const new_single(burst(domain));

                try
                {
                    attempt_single_reaction(*new_single[index]);
                }
                catch (no_space const&)
                {
                    LOG_DEBUG(("pair event single reaction rejected"));
                    ++rejected_moves_;
                }
            }
            break;

        case PAIR_EVENT_COM_ESCAPE:
            {
                LOG_DEBUG(("=> com_escape"));
                time_type const dt(domain.dt());
                boost::array<position_type, 2> const new_pos(
                    draw_new_positions<draw_on_com_escape>(
                        domain, dt));
                boost::array<boost::shared_ptr<single_type>, 2> const new_single(
                    propagate(domain, new_pos));

                add_event(*new_single[0], SINGLE_EVENT_ESCAPE);
                add_event(*new_single[1], SINGLE_EVENT_ESCAPE);
            }
            break;
        
        case PAIR_EVENT_IV_REACTION:
            {
                LOG_DEBUG(("=> iv_reaction"));
                BOOST_ASSERT(::size(domain.reactions()) == 1);
                reaction_rule_type const& r(domain.reactions()[0]);

                switch (::size(r.get_products()))
                {
                case 0:
                    {
                        (*base_type::world_).remove_particle(domain.particles()[0].first);
                        (*base_type::world_).remove_particle(domain.particles()[1].first);
                        if (base_type::rrec_)
                        {
                            (*base_type::rrec_)(reaction_record_type(
                                r.id(),
                                array_gen<particle_id_pair>(),
                                domain.particles()[0],
                                domain.particles()[1]));
                        }
                    }
                    break;
                case 1:
                    {
                        species_id_type const& new_species_id(r.get_products()[0]);
                        molecule_info_type const new_species(
                            (*base_type::world_).get_molecule_info(new_species_id));

                        // calculate new R
                        position_type const new_com(
                            (*base_type::world_).apply_boundary(
                                draw_on_iv_reaction(
                                    this->rng(),
                                    *base_type::world_).draw_com(
                                        domain, domain.dt())));
                   
                        BOOST_ASSERT(
                            (*base_type::world_).distance(
                                domain.shell().second.position(),
                                new_com) + new_species.radius
                            < shape(domain.shell().second).radius());

                        (*base_type::world_).remove_particle(domain.particles()[0].first);
                        (*base_type::world_).remove_particle(domain.particles()[1].first);

                        particle_id_pair const new_particle(
                            (*base_type::world_).new_particle(
                                new_species_id, new_com).first);
                        boost::shared_ptr<single_type> new_single(
                            create_single(new_particle));
                        add_event(*new_single, SINGLE_EVENT_ESCAPE);

                        if (base_type::rrec_)
                        {
                            // (*base_type::rrec_)(reaction_record_type(
                            //     r.id(),
                            //     array_gen(new_particle.first),
                            //     domain.particles()[0].first,
                            //     domain.particles()[1].first));
                            (*base_type::rrec_)(reaction_record_type(
                                r.id(),
                                array_gen(new_particle),
                                domain.particles()[0],
                                domain.particles()[1]));
                        }
                    }
                    break;
                default:
                    throw not_implemented("num products >= 2 not supported.");
                }
                remove_domain(domain);
            }
            break;
        case PAIR_EVENT_IV_ESCAPE:
            {
                LOG_DEBUG(("=> iv_escape"));
                time_type const dt(domain.dt());
                boost::array<position_type, 2> const new_pos(
                    draw_new_positions<draw_on_iv_escape>(
                        domain, dt));
                boost::array<boost::shared_ptr<single_type>, 2> const new_single(
                    propagate(domain, new_pos));

                add_event(*new_single[0], SINGLE_EVENT_ESCAPE);
                add_event(*new_single[1], SINGLE_EVENT_ESCAPE);
            }
            break;
        }
    }

    void fire_event(multi_event& event)
    {
        multi_type& domain(event.domain());
        domain.step();
        LOG_DEBUG(("fire_multi: last_event=%s", boost::lexical_cast<std::string>(domain.last_event()).c_str()));
        multi_step_count_[domain.last_event()]++; 
        switch (domain.last_event())
        {
        default: /* never get here */ BOOST_ASSERT(0); break;
        case multi_type::REACTION:
            if (base_type::rrec_)
                (*base_type::rrec_)(domain.last_reaction());
            burst(domain);
            break;
        case multi_type::ESCAPE:
            burst(domain);
            break;
        case multi_type::NONE:
            add_event(domain);
            break;
        }
    }

    void fire_event(birth_event& event)
    {
        const reaction_rule_type& rr(event.reaction_rule());
        BOOST_ASSERT(::size(rr.get_products()));
        species_id_type const& sp(rr.get_products()[0]);
        LOG_DEBUG(("fire_birth: product=%s", boost::lexical_cast<std::string>(sp).c_str()));

        try
        {
            molecule_info_type const minfo(
                (*base_type::world_).get_molecule_info(sp));

            //XXX: A cuboidal region is expected here.
            const position_type new_pos(
                this->rng().uniform(0, (*base_type::world_).edge_lengths()[0]),
                this->rng().uniform(0, (*base_type::world_).edge_lengths()[1]),
                this->rng().uniform(0, (*base_type::world_).edge_lengths()[2]));

            const particle_shape_type new_particle(new_pos, minfo.radius);

            clear_volume(new_particle);

            if (!(*base_type::world_).no_overlap(new_particle))
            {
                LOG_INFO(("no space for product particle."));
                throw no_space();
            }

            particle_id_pair pp(
                (*base_type::world_).new_particle(sp, new_pos).first);

            if (base_type::rrec_)
            {
                (*base_type::rrec_)(
                    reaction_record_type(rr.id(), array_gen(pp)));
            }

            boost::shared_ptr<single_type> single(create_single(pp));
            add_event(*single, SINGLE_EVENT_ESCAPE);
        }
        catch (no_space const&)
        {
            LOG_DEBUG(("birth reaction rejected."));
            ++rejected_moves_;
        }

        add_event(rr);
    }

    void fire_event(event_type& event)
    {
        {
            single_event* _event(dynamic_cast<single_event*>(&event));
            if (_event)
            {
                fire_event(*_event);
                return;
            }
        }
        {
            pair_event* _event(dynamic_cast<pair_event*>(&event));
            if (_event)
            {
                fire_event(*_event);
                return;
            }
        }
        {
            multi_event* _event(dynamic_cast<multi_event*>(&event));
            if (_event)
            {
                fire_event(*_event);
                return;
            }
        }
        {
            birth_event* _event(dynamic_cast<birth_event*>(&event));
            if (_event)
            {
                fire_event(*_event);
                return;
            }
        }
        throw not_implemented(std::string("unsupported domain type"));
    }

    void _step()
    {
        if (base_type::paranoiac_)
            BOOST_ASSERT(check());

        ++base_type::num_steps_;

        (*dynamic_cast<ReactionRecorderWrapper<reaction_record_type>*>(
            base_type::rrec_.get())).clear();

        if (scheduler_.size() == 0)
        {
            this->set_t(scheduler_.next_time());
            return;
        }

        event_id_pair_type ev(scheduler_.pop());
        this->set_t(ev.second->time());

        LOG_INFO(("%d: t=%.16g dt=%.16g domain=%s rejectedmoves=%d",
                  base_type::num_steps_, this->t(), base_type::dt_,
                  boost::lexical_cast<std::string>(dynamic_cast<domain_event_base const*>(ev.second.get())->domain()).c_str(),
                  rejected_moves_));

        fire_event(*ev.second);

        time_type const next_time(scheduler_.top().second->time());
        base_type::dt_ = next_time - this->t();

        if (base_type::dt_ == 0.)
        {
            ++zero_step_count_;
            if (zero_step_count_ >= std::max(scheduler_.size() * 3, static_cast<std::size_t>(10u)))
            {
                throw illegal_state("too many dt=zero steps. simulator halted?");
            }
        }
        else
        {
            zero_step_count_ = 0;
        }
    }

    static domain_kind get_domain_kind(domain_type const& domain)
    {
        struct domain_kind_visitor: ImmutativeDomainVisitor<traits_type>
        {
            virtual ~domain_kind_visitor() {}

            virtual void operator()(multi_type const&) const
            {
                retval = MULTI;
            }

            virtual void operator()(spherical_single_type const&) const
            {
                retval = SPHERICAL_SINGLE;
            }

            virtual void operator()(cylindrical_single_type const&) const
            {
                retval = CYLINDRICAL_SINGLE;
            }

            virtual void operator()(spherical_pair_type const&) const
            {
                retval = SPHERICAL_PAIR;
            }

            virtual void operator()(cylindrical_pair_type const&) const
            {
                retval = CYLINDRICAL_PAIR;
            }

            domain_kind_visitor(domain_kind& retval): retval(retval) {}

            domain_kind& retval;
        };

        domain_kind retval = NONE;
        domain.accept(domain_kind_visitor(retval));
        return retval;
    }

    static domain_kind get_domain_kind(spherical_single_type const&)
    {
        return SPHERICAL_SINGLE;
    }

    static domain_kind get_domain_kind(cylindrical_single_type const&)
    {
        return CYLINDRICAL_SINGLE;
    }

    static domain_kind get_domain_kind(spherical_pair_type const&)
    {
        return SPHERICAL_PAIR;
    }

    static domain_kind get_domain_kind(cylindrical_pair_type const&)
    {
        return CYLINDRICAL_PAIR;
    }

    static domain_kind get_domain_kind(multi_type const&)
    {
        return MULTI;
    }

    void dump_events() const
    {
        LOG_INFO(("QUEUED EVENTS:"));
        BOOST_FOREACH (event_id_pair_type const& ev, scheduler_.events())
        {
            LOG_INFO(("  #%d: %s", ev.first, stringize_event(*ev.second).c_str()));
        }
    }

    static std::string stringize_event(event_type const& ev)
    {
        {
            single_event const* _ev(dynamic_cast<single_event const*>(&ev));
            if (_ev)
            {
                return stringize_event(*_ev);
            }
        }
        {
            pair_event const* _ev(dynamic_cast<pair_event const*>(&ev));
            if (_ev)
            {
                return stringize_event(*_ev);
            }
        }
        {
            multi_event const* _ev(dynamic_cast<multi_event const*>(&ev));
            if (_ev)
            {
                return stringize_event(*_ev);
            }
        }
        return (boost::format("Event(t=%.16g)") % ev.time()).str();
    }

    static std::string stringize_event_kind(enum single_event_kind kind)
    {
        switch (kind)
        {
        default: /* never get here */ BOOST_ASSERT(0); break;
        case SINGLE_EVENT_ESCAPE:
            return "escape";

        case SINGLE_EVENT_REACTION:
            return "reaction";
        }
    }

    static std::string stringize_event_kind(enum pair_event_kind kind)
    {
        switch (kind)
        {
        default: /* never get here */ BOOST_ASSERT(0); break;
        case PAIR_EVENT_SINGLE_REACTION_0:
            return "reaction(0)";

        case PAIR_EVENT_SINGLE_REACTION_1:
            return "reaction(1)";

        case PAIR_EVENT_COM_ESCAPE:
            return "com_escape";

        case PAIR_EVENT_IV_UNDETERMINED:
            return "iv_undetermined";

        case PAIR_EVENT_IV_ESCAPE:
            return "iv_escape";

        case PAIR_EVENT_IV_REACTION:
            return "iv_reaction";
        }
    }

    static std::string stringize_event(single_event const& ev)
    {
        return (boost::format("SingleEvent(t=%.16g, kind=%s, domain=%s)") %
            ev.time() % stringize_event_kind(ev.kind()) %
            boost::lexical_cast<std::string>(ev.domain())).str();
    }

    static std::string stringize_event(pair_event const& ev)
    {
        return (boost::format("PairEvent(t=%.16g, kind=%s, domain=%s)") %
            ev.time() % stringize_event_kind(ev.kind()) %
            boost::lexical_cast<std::string>(ev.domain())).str();
    }

    static std::string stringize_event(multi_event const& ev)
    {
        return (boost::format("MultiEvent(t=%.16g, domain=%s)") %
            ev.time() % boost::lexical_cast<std::string>(ev.domain())).str();
    }

    template<typename T>
    bool check_domain(AnalyticalSingle<traits_type, T> const& domain) const
    {
        LOG_DEBUG(("checking domain %s", boost::lexical_cast<std::string>(domain).c_str()));
        bool retval(true);
        std::pair<domain_id_type, length_type> closest(
            get_closest_domain(domain.position(), array_gen(domain.id())));
        CHECK(shape_size(shape(domain.shell().second)) <= user_max_shell_size_);
        CHECK(shape_size(shape(domain.shell().second)) <= max_shell_size());
        CHECK(closest.second > shape_size(shape(domain.shell().second)));
        return retval;
    }

    template<typename T>
    bool check_domain(AnalyticalPair<traits_type, T> const& domain) const
    {
        LOG_DEBUG(("checking domain %s", boost::lexical_cast<std::string>(domain).c_str()));
        bool retval(true);
        std::pair<domain_id_type, length_type> closest(
            get_closest_domain(domain.position(), array_gen(domain.id())));
        CHECK(shape_size(shape(domain.shell().second)) <= user_max_shell_size_);
        CHECK(shape_size(shape(domain.shell().second)) <= max_shell_size());
        CHECK(closest.second > shape_size(shape(domain.shell().second)));
        return retval;
    }

    bool check_domain(multi_type const& domain) const
    {
        LOG_DEBUG(("checking domain %s", boost::lexical_cast<std::string>(domain).c_str()));
        bool retval(true);
        BOOST_FOREACH (typename multi_type::spherical_shell_id_pair const& shell,
                       domain.get_shells())
        {
            std::pair<domain_id_type, length_type> closest(
                get_closest_domain(shape_position(shape(shell.second)),
                                   array_gen(domain.id())));
            CHECK(shape_size(shape(shell.second)) <= user_max_shell_size_);
            CHECK(shape_size(shape(shell.second)) <= max_shell_size());
            CHECK(closest.second > shape_size(shape(shell.second)));
        }
        return retval;
    }

    bool check_domain(domain_type const& domain) const
    {
        struct visitor: public ImmutativeDomainVisitor<traits_type>
        {
            virtual ~visitor() {}

            virtual void operator()(multi_type const& domain) const
            {
                retval = self.check_domain(domain);
            }

            virtual void operator()(spherical_single_type const& domain) const
            {
                retval = self.check_domain(domain);
            }

            virtual void operator()(cylindrical_single_type const& domain) const
            {
                retval = self.check_domain(domain);
            }

            virtual void operator()(spherical_pair_type const& domain) const
            {
                retval = self.check_domain(domain);
            }

            virtual void operator()(cylindrical_pair_type const& domain) const
            {
                retval = self.check_domain(domain);
            }

            visitor(EGFRDSimulator const& self, bool& retval)
                : self(self), retval(retval) {}

            EGFRDSimulator const& self;
            bool& retval;
        };

        bool retval;
        domain.accept(visitor(*this, retval));
        return retval;
    }

    bool check_overlap(particle_shape_type const& s) const
    {
        const particle_id_pair_and_distance_list overlapped(
            (*base_type::world_).check_overlap(s));

        if (overlapped.size() > 0)
        {
            LOG_DEBUG(("check_overlap %s failed:",
                boost::lexical_cast<std::string>(s).c_str()));
            dump_overlapped(overlapped);
            return false;
        }
        return true;
    }

    bool check_overlap(particle_shape_type const& s, particle_id_type const& ignore) const
    {
        const particle_id_pair_and_distance_list overlapped(
            (*base_type::world_).check_overlap(s, ignore));

        if (overlapped.size() > 0)
        {
            LOG_DEBUG(("check_overlap %s failed:",
                boost::lexical_cast<std::string>(s).c_str());
            dump_overlapped(overlapped));
            return false;
        }
        return true;
    }

    bool check_overlap(particle_shape_type const& s, particle_id_type const& ignore1, particle_id_type const& ignore2) const
    {
        const particle_id_pair_and_distance_list overlapped(
            (*base_type::world_).check_overlap(s, ignore1, ignore2));

        if (overlapped.size() > 0)
        {
            LOG_DEBUG(("check_overlap %s failed:",
                boost::lexical_cast<std::string>(s).c_str()));
            dump_overlapped(overlapped);
            return false;
        }
        return true;
    }

    void dump_overlapped(particle_id_pair_and_distance_list const& list)const
    {
        if (log_.level() == Logger::L_DEBUG)
        {
            BOOST_FOREACH (particle_id_pair_and_distance const& i, list)
            {
                log_.debug("  (%s:%s) %.16g",
                    boost::lexical_cast<std::string>(i.first.first).c_str(),
                    boost::lexical_cast<std::string>(i.first.second).c_str(),
                    i.second);
            }
        }
    }

    static rate_type calculate_k_tot(reaction_rules const& rules)
    {
        rate_type k_tot(0.);
        BOOST_FOREACH (reaction_rule_type const& rule, rules)
        {
            k_tot += rule.k();
        }
        return k_tot;
    }

    reaction_rule_type const& draw_reaction_rule(reaction_rules const& rules)
    {
        const rate_type k_tot(calculate_k_tot(rules));
        if(k_tot == std::numeric_limits<rate_type>::infinity())
        {
            LOG_WARNING(("k_tot == infinite: first reaction type applied."));
            return rules[0];
        }

        const rate_type t(this->rng().uniform(0., 1.) * k_tot);
        rate_type a(0.);
        BOOST_FOREACH(reaction_rule_type const& r, rules)
        {
            a += r.k();
            if (a > t)
                return r;
        }

        BOOST_ASSERT(false); // should never happen
    }

    //template<typename T1, typename T2>
    static position_type
    //adjust_iv_with_old_iv(T1 const& new_iv, T2 const& old_iv)
    adjust_iv_with_old_iv(position_type const& new_iv, position_type const& old_iv)
    {
        length_type const angle(std::acos(old_iv[2] / length(old_iv)));
        if (std::fmod(angle, M_PI) != 0.0)
        {
            position_type const rotation_axis(
                normalize(position_type(-old_iv[1], old_iv[0], 0.)));
            return rotate_vector(new_iv, rotation_axis, angle);
        }
        else if (angle == 0.)
        {
            return new_iv;
        }
        else
        {
            return position_type(new_iv[0], new_iv[1], -new_iv[1]); 
        }
    }

    template<typename T>
    bool check_pair_pos(AnalyticalPair<traits_type, T> const& domain,
                        boost::array<particle_id_pair, 2> const& new_particles)
    {
        length_type const new_distance(
            (*base_type::world_).distance(new_particles[0].second.position(),
                                          new_particles[1].second.position()));
        length_type const r01(new_particles[0].second.radius() +
                              new_particles[1].second.radius());

        if (new_distance <= r01)
        {
            log_.warn(
                "rejected move: pair=%s, radii=%.16g, particle_distance=%.16g",
                boost::lexical_cast<std::string>(domain).c_str(),
                r01, new_distance);
            return false;
        }

        // particles within mobility radius.
        position_type const& com(shape(domain.shell().second).position());
        length_type const radius(shape(domain.shell().second).radius());
        length_type const d[2] = {
            (*base_type::world_).distance(com, new_particles[0].second.position()) + new_particles[0].second.radius(),
            (*base_type::world_).distance(com, new_particles[1].second.position()) + new_particles[1].second.radius()
        };
        if (d[0] > radius || d[1] > radius)
        {
            log_.warn(
                "rejected move: new particle(s) out of protective sphere: pair=%s, radii=%.16g, d0=%.16g, d1=%.16g",
                boost::lexical_cast<std::string>(domain).c_str(),
                d[0], d[1]);
            return false;
        }
        return true;
    }

    template<typename T>
    static greens_functions::PairGreensFunction* choose_pair_greens_function(
            AnalyticalPair<traits_type, T> const& domain, time_type t)
    {
        length_type const r0(domain.r0());
        length_type const distance_from_sigma(r0 - domain.sigma());
        length_type const distance_from_shell(domain.a_r() - r0);
        length_type const threshold_distance(
            traits_type::CUTOFF_FACTOR * std::sqrt(6. * domain.D_tot() * t));

        BOOST_ASSERT(::size(domain.reactions()) == 1);
        if (distance_from_sigma < threshold_distance)
        {
            if (distance_from_shell < threshold_distance)
            {
                // near both a and sigma;
                // use greens_functions::GreensFunction3DRadAbs
                LOG_DEBUG(("GF: normal"));
                return new greens_functions::GreensFunction3DRadAbs(
                    domain.D_tot(), domain.reactions()[0].k(),
                    r0, domain.sigma(), domain.a_r());
            }
            else
            {
                // near sigma; use greens_functions::GreensFunction3DRadInf
                LOG_DEBUG(("GF: only sigma"));
                return new greens_functions::GreensFunction3DRadInf(
                    domain.D_tot(), domain.reactions()[0].k(),
                    r0, domain.sigma());
            }
        }
        else
        {
            if (distance_from_shell < threshold_distance)
            {
                // near a;
                LOG_DEBUG(("GF: only a"));
                return new greens_functions::GreensFunction3DAbs(
                    domain.D_tot(), r0, domain.a_r());
            }
            else
            {
                // distant from both a and sigma; 
                LOG_DEBUG(("GF: free"));
                return new greens_functions::GreensFunction3D(domain.D_tot(), r0);
            }
        }
    }

    static length_type calculate_single_shell_size(
            single_type const& single,
            single_type const& closest,
            length_type distance,
            length_type shell_distance)
    {
        length_type const min_radius0(single.particle().second.radius());
        D_type const D0(single.particle().second.D());
        if (D0 == 0)
            return min_radius0;
        length_type const min_radius1(closest.particle().second.radius());
        D_type const D1(closest.particle().second.D());
        length_type const min_radius01(min_radius0 + min_radius1);
        length_type const sqrtD0(std::sqrt(D0));

        return std::max(std::min(sqrtD0 / (sqrtD0 + std::sqrt(D1))
                            * (distance - min_radius01) + min_radius0,
                            shell_distance / traits_type::SAFETY),
                        min_radius0);
    }

protected:
    double const bd_dt_factor_;
    int const num_retries_;
    length_type const user_max_shell_size_;

    domain_map domains_;
    boost::scoped_ptr<spherical_shell_matrix_type> ssmat_;
    boost::scoped_ptr<cylindrical_shell_matrix_type> csmat_;
    shell_matrix_map_type smatm_;
    shell_id_generator shidgen_;
    domain_id_generator didgen_;
    event_scheduler_type scheduler_;
    boost::array<int, NUM_SINGLE_EVENT_KINDS> single_step_count_;
    boost::array<int, NUM_PAIR_EVENT_KINDS> pair_step_count_;
    boost::array<int, multi_type::NUM_MULTI_EVENT_KINDS> multi_step_count_;
    boost::array<int, NUM_DOMAIN_KINDS> domain_count_per_type_;
    length_type single_shell_factor_;
    length_type multi_shell_factor_;
    unsigned int rejected_moves_;
    unsigned int zero_step_count_;
    bool dirty_;
    static Logger& log_;
};
#undef CHECK

template<typename Ttraits>
inline char const* retrieve_domain_type_name(
    typename EGFRDSimulator<Ttraits>::spherical_single_type const&)
{
    return "SphericalSingle";
}

template<typename Ttraits>
inline char const* retrieve_domain_type_name(
    typename EGFRDSimulator<Ttraits>::cylindrical_single_type const&)
{
    return "CylindricalSingle";
}

template<typename Ttraits>
inline char const* retrieve_domain_type_name(
    typename EGFRDSimulator<Ttraits>::spherical_pair_type const&)
{
    return "SphericalPair";
}

template<typename Ttraits>
inline char const* retrieve_domain_type_name(
    typename EGFRDSimulator<Ttraits>::cylindrical_pair_type const&)
{
    return "CylindricalPair";
}



template<typename Ttraits_>
Logger& EGFRDSimulator<Ttraits_>::log_(Logger::get_logger("ecell.EGFRDSimulator"));

#endif /* EGFRDSIMULATOR_HPP */
