#ifndef EGFRDSIMULATOR_HPP
#define EGFRDSIMULATOR_HPP

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include "utils/get_mapper_mf.hpp"
#include "utils/fun_composition.hpp"
#include "utils/fun_wrappers.hpp"
#include "utils/pair.hpp"
#include "ShellID.hpp"
#include "DomainID.hpp"
#include "Shell.hpp"
#include "PairGreensFunction.hpp"
#include "ParticleSimulator.hpp"
#include "MatrixSpace.hpp"
#include "AnalyticalSingle.hpp"
#include "AnalyticalPair.hpp"
#include "Multi.hpp"

template<typename Tworld_>
struct EGFRDSimulatorTraitsBase: public ParticleSimulatorTraitsBase<Tworld_>
{
    typedef Tworld_ world_type;
    typedef ShellID shell_id_type;
    typedef DomainID domain_id_type;
    typedef typename ParticleSimulatorTraitsBase<Tworld_>::sphere_type sphere_type;
    typedef typename ParticleSimulatorTraitsBase<Tworld_>::cylinder_type cylinder_type;
    typedef SerialIDGenerator<shell_id_type> shell_id_generator;
    typedef SerialIDGenerator<domain_id_type> domain_id_generator;
    typedef Shell<sphere_type, domain_id_type> spherical_shell_type;
    typedef Shell<cylinder_type, domain_id_type> cylindrical_shell_type;
    typedef Domain<EGFRDSimulatorTraitsBase> domain_type;
    typedef std::pair<const domain_id_type, domain_type*> domain_id_pair;
    typedef std::pair<const shell_id_type, spherical_shell_type> spherical_shell_id_pair;
    typedef std::pair<const shell_id_type, cylindrical_shell_type> cylindrical_shell_id_pair;
    typedef int event_id_type;
    typedef EventType event_kind_type;
};

template<typename Ttraits_>
class EGFRDSimulator: public ParticleSimulator<Ttraits_>
{
public:
    typedef Ttraits_ traits_type;
    typedef ParticleSimulator<Ttraits_> base_type;
    typedef typename traits_type::world_type world_type;
    typedef typename traits_type::domain_id_type domain_id_type;
    typedef typename traits_type::shell_id_type shell_id_type;
    typedef typename traits_type::spherical_shell_id_pair spherical_shell_id_pair;
    typedef typename traits_type::cylindrical_shell_id_pair cylindrical_shell_id_pair;
    typedef typename traits_type::shell_id_generator shell_id_generator;
    typedef typename traits_type::domain_id_generator domain_id_generator;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename traits_type::rng_type rng_type;

    typedef typename traits_type::spherical_shell_type spherical_shell_type;
    typedef typename traits_type::cylindrical_shell_type cylindrical_shell_type;
    typedef typename traits_type::domain_type domain_type;
    typedef typename traits_type::domain_id_pair domain_id_pair;

    typedef AnalyticalSingle<traits_type, spherical_shell_type> spherical_single_type;
    typedef AnalyticalSingle<traits_type, cylindrical_shell_type> cylindrical_single_type;
    typedef AnalyticalPair<traits_type, spherical_shell_type> spherical_pair_type;
    typedef AnalyticalPair<traits_type, cylindrical_shell_type> cylindrical_pair_type;
    typedef Multi<EGFRDSimulator> multi_type;

protected:
    typedef MatrixSpace<spherical_shell_type, shell_id_type, get_mapper_mf> spherical_shell_matrix_type;
    typedef MatrixSpace<cylindrical_shell_type, shell_id_type, get_mapper_mf> cylindrical_shell_matrix_type;
    typedef typename get_mapper_mf<domain_id_type, domain_type*>::type domain_map;

public:
    typedef abstract_limited_generator<domain_id_pair> domain_id_pair_generator;

public:
    virtual ~EGFRDSimulator()
    {
        //std::for_each(domains_.begin(), domains_.end(),
        //    compose_unary(delete_ptr<domain_type>(),
        //                  select_second<typename domain_map::value_type>()));
    }

    EGFRDSimulator(world_type& world, rng_type& rng,
                   network_rules_type const& network_rules)
        : base_type(world, rng, network_rules),
          ssmat_(world.world_size(), world.matrix_size()),
          csmat_(world.world_size(), world.matrix_size()) {}

    spherical_shell_id_pair
    new_spherical_shell(domain_id_type const& did,
                        typename spherical_shell_type::shape_type const& shape)
    {
        spherical_shell_id_pair retval(sidgen_(),
            spherical_shell_type(did, shape));
        ssmat_.update(retval);
        return retval;
    }

    cylindrical_shell_id_pair
    new_cylindrical_shell(domain_id_type const& did,
                        typename cylindrical_shell_type::shape_type const& shape)
    {
        cylindrical_shell_id_pair retval(sidgen_(),
            cylindrical_shell_type(did, shape));
        ssmat_.update(retval);
        return retval;
    }

    spherical_shell_id_pair const& get_spherical_shell(shell_id_type const& id) const
    {
        return ssmat_[id];
    }

    cylindrical_shell_id_pair const& get_cylindrical_shell(shell_id_type const& id) const
    {
        return csmat_[id];
    }

    domain_type* get_domain(domain_id_type const& id) const
    {
        typename domain_map::const_iterator i(domains_.find(id));

        if (i == domains_.end())
        {
            throw not_found(
                (boost::format("domain id #%s not found") % boost::lexical_cast<std::string>(id)).str());
        }

        return (*i).second;
    }

    bool update_domain(domain_id_pair const& dp)
    {
        typename domain_map::iterator i(domains_.find(dp.first));
        if (domains_.end() != i)
        {
            (*i).second = dp.second;
            return false;
        }
        else
        {
            domains_.insert(dp);
            return true;
        }
    }

    domain_id_pair_generator* get_domains() const
    {
        return make_range_generator<domain_id_pair>(domains_);
    }

    virtual void initialize()
    {
        domains_.clear();
        ssmat_.clear();
        csmat_.clear();
    }

protected:
    domain_map domains_;
    spherical_shell_matrix_type ssmat_;
    cylindrical_shell_matrix_type csmat_;
    shell_id_generator sidgen_;
    domain_id_generator didgen_;
};


#endif /* EGFRDSIMULATOR_HPP */
