#ifndef ECELL4_ODE_ODE_SIMULATOR_NEW_HPP
#define ECELL4_ODE_ODE_SIMULATOR_NEW_HPP

#include <cstring>
#include <vector>
#include <numeric>
#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/get_mapper_mf.hpp>

#include <ecell4/core/SimulatorBase.hpp>
#include <ecell4/core/NetworkModel.hpp>
//#include "ODEWorld.hpp"
#include <ecell4/core/Species.hpp>
#include <ecell4/core/Context.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Space.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/NetfreeModel.hpp>

#include "ODESimulator.hpp"

namespace ecell4
{
namespace ode
{

class ODEWorld_New
    : public Space
{
public:
    typedef std::vector<Real> num_molecules_container_type;
    typedef std::vector<Species> species_container_type;
    typedef utils::get_mapper_mf<
        Species, num_molecules_container_type::size_type>::type species_map_type;

public:
    ODEWorld_New(const Real3& edge_lengths = Real3(1, 1, 1))
        : t_(0.0)
    {
        reset(edge_lengths);
    }

    ODEWorld_New(const std::string& filename)
        : t_(0.0)
    {
        reset(Real3(1, 1, 1));
        //this->load(filename);
    }

    // Space Traints
    const Real t() const
    {
        return t_;
    }

    void set_t(const Real& t)
    {
        if (t < 0.0)
        {
            throw std::invalid_argument("the time must be positive.");
        }
        t_ = t;
    }

    const Real3& edge_lengths() const
    {
        return edge_lengths_;
    }

    void reset(const Real3& edge_lengths)
    {
        t_ = 0.0;
        index_map_.clear();
        num_molecules_.clear();
        species_.clear();

        for (Real3::size_type dim(0); dim < 3; ++dim)
        {
            if (edge_lengths[dim] <= 0)
            {
                throw std::invalid_argument("the edge length must be positive.");
            }
        }

        edge_lengths_ = edge_lengths;
        volume_ = edge_lengths[0] * edge_lengths[1] * edge_lengths[2];
    }

    const Real volume() const
    {
        return volume_;
    }

    void set_volume(const Real& volume)
    {
        if (volume <= 0.0)
        {
            throw std::invalid_argument("The volume must be positive.");
        }

        volume_ = volume;
        const Real L(cbrt(volume));
        edge_lengths_ = Real3(L, L, L);
    }

    Integer num_molecules(const Species& sp) const
    {
        return static_cast<Integer>(get_value(sp));
    }

    Integer num_molecules_exact(const Species& sp) const
    {
        return static_cast<Integer>(get_value_exact(sp));
    }

    std::vector<Species> list_species() const
    {
        return species_;
    }

    // CompartmentSpace member functions

    void add_molecules(const Species& sp, const Real& num)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            reserve_species(sp);
            i = index_map_.find(sp);
        }

        num_molecules_[(*i).second] += num;
    }

    void remove_molecules(const Species& sp, const Real& num)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            throw NotFound("Species not found");
        }

        num_molecules_[(*i).second] -= num;
    }

    // Optional members

    Real get_value(const Species& sp) const
    {
        SpeciesExpressionMatcher sexp(sp);
        Real retval(0);
        for (species_map_type::const_iterator i(index_map_.begin());
            i != index_map_.end(); ++i)
        {
            if (sexp.match((*i).first))
            {
                do
                {
                    retval += num_molecules_[(*i).second];
                } while (sexp.next());
            }
        }
        return retval;
    }

    Real get_value_exact(const Species& sp) const
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            // throw NotFound("Species not found");
            return 0.0;
        }

        return num_molecules_[(*i).second];
    }

    void set_value(const Species& sp, const Real& num)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            reserve_species(sp);
            i = index_map_.find(sp);
        }

        num_molecules_[(*i).second] = num;
    }

    //FIXME Following 2 member functions are not implemented for now!
    void save(const std::string& filename) const;
    void load(const std::string& filename);

    bool has_species(const Species& sp) const
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        return (i != index_map_.end());
    }

    void reserve_species(const Species& sp)
    {
        species_map_type::const_iterator i(index_map_.find(sp));
        if (i != index_map_.end())
        {
            throw AlreadyExists("Species already exists");
        }

        index_map_.insert(std::make_pair(sp, num_molecules_.size()));
        species_.push_back(sp);
        num_molecules_.push_back(0);
    }

    void release_species(const Species& sp)
    {
        species_map_type::iterator i(index_map_.find(sp));
        if (i == index_map_.end())
        {
            throw NotFound("Species not found");
        }

        species_map_type::mapped_type
            idx((*i).second), last_idx(num_molecules_.size() - 1);
        if (idx != last_idx)
        {
            const species_container_type::size_type
                idx_(static_cast<species_container_type::size_type>(idx)),
                last_idx_(
                    static_cast<species_container_type::size_type>(last_idx));
            const Species& last_sp(species_[last_idx_]);
            species_[idx_] = last_sp;
            num_molecules_[idx] = num_molecules_[last_idx];
            index_map_[last_sp] = idx;
        }

        species_.pop_back();
        num_molecules_.pop_back();
        index_map_.erase(sp);
    }

    void bind_to(boost::shared_ptr<Model> model);
    //void bind_to(boost::shared_ptr<NetworkModel> model);

    //boost::shared_ptr<NetworkModel> lock_model() const
    boost::shared_ptr<Model> lock_model() const
    {
        //if (generated_)
        //{
        //    return generated_;
        //}
        //else
        //{
        //    return model_.lock();
        //}
        return model_.lock();
    }
    //
    void add_molecules(const Species& sp, const Integer& num,
        const boost::shared_ptr<Shape> shape)
    {
        add_molecules(sp, num);
    }

    std::pair<std::pair<ParticleID, Particle>, bool> new_particle(const Particle& p)
    {
        add_molecules(p.species(), 1);
        return std::make_pair(std::make_pair(ParticleID(), p), true);
    }

    std::pair<std::pair<ParticleID, Particle>, bool> new_particle(
        const Species& sp, const Real3& pos)
    {
        add_molecules(sp, 1);
        return std::make_pair(
            std::make_pair(ParticleID(), Particle(sp, pos, 0.0, 0.0)), true);
    }

    //Real evaluate(const ReactionRule& rr) const
    //{
    //    ODEReactionRule oderr(rr);
    //    return evaluate(oderr);
    //}

    //Real evaluate(ODEReactionRule& rr) const;

protected:
    Real3 edge_lengths_;
    Real volume_;
    Real t_;

    num_molecules_container_type num_molecules_;
    species_container_type species_;
    species_map_type index_map_;

    boost::weak_ptr<Model> model_;
    //boost::shared_ptr<NetworkModel> generated_;
};

class ODESimulator_New
    : public SimulatorBase<Model, ODEWorld_New>
{
public:
    typedef SimulatorBase<Model, ODEWorld_New> base_type;
public:
    typedef boost::numeric::ublas::vector<Real> state_type;
    typedef boost::numeric::ublas::matrix<Real> matrix_type;
    typedef std::vector<state_type::size_type> index_container_type;
    typedef std::vector<Real> coefficient_container_type;
    typedef ReactionRule reaction_container_type;
    typedef Model::reaction_rule_container_type reaction_rule_container_type;

    struct mapped_reaction_type
    {
        index_container_type reactants;
        coefficient_container_type reactant_coefficients;
        index_container_type products;
        coefficient_container_type product_coefficients;
        Real k;
        //boost::weak_ptr<ODERatelaw> ratelaw;
        const ReactionRule *raw;
        boost::weak_ptr<ReactionRuleDescriptor> rr_descriptor;
    };
    typedef std::vector<mapped_reaction_type> mapped_reaction_container_type;

    class deriv_func 
    {
    public:
        deriv_func(const mapped_reaction_container_type &reactions, const Real &volume)
            : reactions_(reactions), volume_(volume), vinv_(1./volume)
        {;}
        void operator()(const state_type &x, state_type &dxdt, const double &t)
        {
            std::fill(dxdt.begin(), dxdt.end(), 0.);
            for(mapped_reaction_container_type::const_iterator mapped_rr_it(reactions_.begin());
                    mapped_rr_it != reactions_.end(); mapped_rr_it++) 
            {   // Calculate the flux of each ReactionRule
                //  FIXME.  Currently, this can handle only the mass action
                Real flux = (mapped_rr_it->k * volume_);
                //for (index_container_type::const_iterator r = mapped_rr_it->reactants.begin(); 
                //        r != mapped_rr_it->reactants.end(); r++) {
                //    //flux *= std::pow( x[*r] * vinv_, 1. );
                //    flux *=  x[*r] * vinv_;
                //}
                for( index_container_type::size_type i = 0; i < mapped_rr_it->reactants.size(); i++) {
                    index_container_type::value_type reactant_index = mapped_rr_it->reactants[i];
                    flux *= std::pow(x[reactant_index] * vinv_, mapped_rr_it->reactant_coefficients[i]);
                }

                // Merge each reactions's flux 
                //for(index_container_type::const_iterator r = mapped_rr_it->reactants.begin();
                //        r != mapped_rr_it->reactants.end(); r++) 
                //{
                //    dxdt[*r] -= (flux); // 
                //    //dxdt[*j] -= (flux * (double)i->reactant_coefficients[nth]);
                //}
                //for(index_container_type::const_iterator p = mapped_rr_it->products.begin();
                //        p != mapped_rr_it->products.end(); p++)
                //{
                //    dxdt[*p] += (flux);
                //}
                for( index_container_type::size_type i = 0; i < mapped_rr_it->reactants.size(); i++) {
                    index_container_type::value_type reactant_index = mapped_rr_it->reactants[i];
                    dxdt[reactant_index] -= (flux * (double)mapped_rr_it->reactant_coefficients[i]);
                }
                for( index_container_type::size_type i = 0; i < mapped_rr_it->products.size(); i++) {
                    index_container_type::value_type product_index = mapped_rr_it->products[i];
                    dxdt[product_index] += (flux * (double)mapped_rr_it->product_coefficients[i]);
                }
            }
            //std::cout << "========================================" << std::endl;
            //std::cout << "dxdt dumping at t = " << t << std::endl;
            //for (state_type::const_iterator it = dxdt.begin();
            //        it != dxdt.end(); it++) 
            //{
            //    std::cout << *it << std::endl;
            //}
            return;
        }
    protected:
        const mapped_reaction_container_type reactions_;
        const Real volume_;
        const Real vinv_;
    };

    class jacobi_func
    {
    public:
        jacobi_func(
            //const reaction_container_type &reactions, const Real& volume,
            const mapped_reaction_container_type &reactions, const Real& volume,
            const Real& abs_tol, const Real& rel_tol)
            : reactions_(reactions), volume_(volume), vinv_(1.0 / volume), abs_tol_(abs_tol), rel_tol_(rel_tol)
        {
            ;
        }
        void operator()(
                const state_type& x, matrix_type& jacobi, const double &t, state_type &dfdt) const
        {
            //fill 0 into jacobi and dfdt
            std::fill(dfdt.begin(), dfdt.end(), 0.0);
            std::fill(jacobi.data().begin(), jacobi.data().end(), 0.0);

            // const Real ETA(2.2204460492503131e-16);
            const Real SQRTETA(1.4901161193847656e-08);
            const Real r0(1.0);
            // XXX Following comment-outed code is derived from previous ODESimulator.hpp
            //          Should be erased??????
            // Real fac(0.0);
            // for (std::size_t k(0); k < dfdt.size(); ++k)
            // {
            //     const Real ewtk(atol + rtol * x[k]);
            //     fac = std::max(fac, dfdt[k] * ewtk);
            // }
            // const Real r0(1000.0 * h * ETA * dfdt.size() * fac);  //XXX: h means the step interval
            // {
            //     const Real ewtj(atol + rtol * x[j]);
            //     const Real dyj(std::max(SQRTETA * abs(x[j]), r0 * ewtj));
            // }

            // const Real h(1.0e-8);
            // const Real ht(1.0e-10);
            const Real ht(1.0e-10);

            //for(mapped_reaction_container_type::const_iterator mapped_rr_it = (reactions_.begin());
            //        mapped_rr_it != reactions_.end(); mapped_rr_it++) 
            //{
            //    // 1. Calculate the current flux. This is the same as the flux calculated in deriv_func.
            //    Real flux0 = (mapped_rr_it->k * volume_);
            //    for (index_container_type::const_iterator r = mapped_rr_it->reactants.begin(); 
            //            r != mapped_rr_it->reactants.end(); r++)
            //    {
            //        //flux *= std::pow( x[*r] * vinv_, 1. );
            //        flux0 *=  x[*r] * vinv_;
            //    }
            //    // 2. Differentiate by time
            //    {
            //        // In the mass reactio, this become ZERO.
            //        //  Currently, not implemented.
            //    }
            //    // Differentiate by each reactants
            //    {
            //        for(std::size_t i(0); i < mapped_rr_it->reactants.size(); i++) {

            //        }
            //            const Real ewt = abs_tol_ + rel_tol_ * abs(reactants_states[j]);
            //            const Real h = std::max(SQRTETA * abs(reactants_states[j]), r0 * ewt);
            //    }
            //    // Differentiate by each products
            //}
        }
    protected:
        const mapped_reaction_container_type reactions_;
        const Real volume_;
        const Real vinv_;
        const Real abs_tol_, rel_tol_;
    };

    struct StateAndTimeBackInserter
    {
        typedef std::vector<state_type> state_container_type;
        typedef std::vector<Real> time_container_type;
        // Constructor
        StateAndTimeBackInserter(state_container_type &states, time_container_type &times)
            : m_states(states), m_times(times)
        {;}

        void operator()(const state_type &x, double t)
        {
            m_states.push_back(x);
            m_times.push_back(t);
        }
        state_container_type &m_states;
        time_container_type &m_times;
    };

    ODESimulator_New(
        const boost::shared_ptr<Model>& model,
        const boost::shared_ptr<ODEWorld_New> &world,
        const ODESolverType solver_type = ROSENBROCK4_CONTROLLER)
        :base_type(model, world), dt_(inf), 
         abs_tol_(1e-6), rel_tol_(1e-6), solver_type_(solver_type)
    {;}
    ODESimulator_New(
        const boost::shared_ptr<ODEWorld_New>& world,
        const ODESolverType solver_type = ROSENBROCK4_CONTROLLER)
        : base_type(world), dt_(inf), abs_tol_(1e-6), rel_tol_(1e-6),
          solver_type_(solver_type)
    {
        initialize();
    }

    //ODESimulator_New(
    //    const boost::shared_ptr<Model>& model,
    //    const boost::shared_ptr<ODEWorld_New>& world,
    //    const ODESolverType solver_type = ROSENBROCK4_CONTROLLER)
    //    : base_type(boost::shared_ptr<NetworkModel>(new NetworkModel(model)), world),
    //      dt_(inf), abs_tol_(1e-6), rel_tol_(1e-6), solver_type_(solver_type)
    //{
    //    initialize();
    //}

    void initialize(void)
    {
        const std::vector<Species> species(model_->list_species());
        for(std::vector<Species>::const_iterator it = species.begin(); it != species.end(); it++)
        {
            if (world_->has_species(*it) != true) {
                world_->reserve_species(*it);
            }
        }
    }

    void step(void)
    {
        step(next_time());
        // FIXME Temporary, comment out
        //if ( this->model_->has_network_model() )
        //{
        //    this->model_->update_model();
        //}
    }
    bool step(const Real &upto);

    Real t(void) const
    {
        return world_->t();
    }
    void set_t(const Real &t)
    {
        world_->set_t(t);
    }
    Real dt(void) const
    {
        return this->dt_;
    }
    void set_dt(const Real &dt)
    {
        if (dt <= 0)
        {
            throw std::invalid_argument("The step size must be positive.");
        }
        dt_ = dt;
    }

    Real absolute_tolerance() const
    {
        return abs_tol_;
    }

    void set_absolute_tolerance(const Real abs_tol)
    {
        if (abs_tol < 0)
        {
            throw std::invalid_argument("A tolerance must be positive or zero.");
        }
        abs_tol_ = abs_tol;
    }

    Real relative_tolerance() const
    {
        return rel_tol_;
    }

    void set_relative_tolerance(const Real rel_tol)
    {
        if (rel_tol < 0)
        {
            throw std::invalid_argument("A tolerance must be positive or zero.");
        }
        rel_tol_ = rel_tol;
    }

    //Real evaluate(const ReactionRule& rr) const
    //{
    //    return evaluate(ODEReactionRule(rr));
    //}

    //Real evaluate(const ODEReactionRule& rr) const
    //{
    //    if (!rr.has_ratelaw())
    //    {
    //        // std::cout << "No ratelaw was bound. Return zero." << std::endl;
    //        return 0.0;
    //    }

    //    const ODEReactionRule::reactant_container_type reactants = rr.reactants();
    //    const ODEReactionRule::product_container_type products = rr.products();

    //    ODERatelaw::state_container_type::size_type cnt(0);

    //    ODERatelaw::state_container_type r(reactants.size());
    //    for(ODEReactionRule::reactant_container_type::const_iterator j(reactants.begin());
    //        j != reactants.end(); j++, cnt++)
    //    {
    //        r[cnt] = static_cast<double>(world_->get_value_exact(*j));
    //    }

    //    cnt = 0;

    //    ODERatelaw::state_container_type p(products.size());
    //    for(ODEReactionRule::reactant_container_type::const_iterator j(products.begin());
    //        j != products.end(); j++, cnt++)
    //    {
    //        p[cnt] = static_cast<double>(world_->get_value_exact(*j));
    //    }

    //    return rr.get_ratelaw()->deriv_func(r, p, world_->volume(), world_->t(), rr);
    //}

protected:
    //std::pair<deriv_func, jacobi_func> generate_system() const;
    deriv_func generate_system() const;

protected:
    Real dt_;
    Real abs_tol_, rel_tol_;
    ODESolverType solver_type_;
};

} // ode
} // ecell4

#endif  // ECELL4_ODE_ODE_SIMULATOR2_HPP
