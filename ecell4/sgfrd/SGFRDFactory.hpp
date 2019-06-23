#ifndef ECELL4_SGFRD_SGFRD_FACTORY_HPP
#define ECELL4_SGFRD_SGFRD_FACTORY_HPP
#include <ecell4/core/SimulatorFactory.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Polygon.hpp>
#include "SGFRDWorld.hpp"
#include "SGFRDSimulator.hpp"

namespace ecell4
{
namespace sgfrd
{

class SGFRDFactory :
    public SimulatorFactory<SGFRDWorld, SGFRDSimulator>
{
  public:

    typedef SimulatorFactory<SGFRDWorld, SGFRDSimulator> base_type;
    typedef base_type::world_type     world_type;
    typedef base_type::simulator_type simulator_type;
    typedef SGFRDFactory              this_type;

  public:

    SGFRDFactory(const Integer3& matrix_sizes   = default_matrix_sizes(),
                 Real bd_dt_factor              = default_bd_dt_factor(),
                 Real bd_reaction_length_factor = default_bd_reaction_length_factor())
        : base_type(), rng_(nullptr), polygon_("", STLFormat::Ascii),
          matrix_sizes_(matrix_sizes), bd_dt_factor_(bd_dt_factor),
          bd_reaction_length_factor_(bd_reaction_length_factor)
    {
        ; // do nothing
    }
    virtual ~SGFRDFactory() override = default;

    static inline Integer3 default_matrix_sizes()
    {
        return Integer3(3, 3, 3);
    }

    static inline Real default_bd_dt_factor()
    {
        return 0.01; // relative to the upper limit; for safety
    }

    static inline Real default_bd_reaction_length_factor()
    {
        return 0.1; // 0.05 ~ 0.1
    }

    this_type& rng(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        rng_ = rng;
        return (*this);
    }

    inline this_type* rng_ptr(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        return std::addressof(this->rng(rng));
    }

    this_type& polygon(const std::string& fname, const STLFormat fmt)
    {
        this->polygon_ = std::make_pair(fname, fmt);
        return (*this);
    }
    this_type* polygon_ptr(const std::string& fname, const STLFormat fmt)
    {
        return std::addressof(this->polygon(fname, fmt));
    }

  protected:

    virtual world_type* create_world(const Real3& edge_lengths) const
    {
        if (rng_)
        {
            if(this->polygon_.first.empty())
            {
                return new world_type(edge_lengths, matrix_sizes_, rng_);
            }
            else
            {
                return new world_type(edge_lengths, matrix_sizes_, rng_,
                        this->polygon_.first, this->polygon_.second);
            }
        }
        else
        {
            if(this->polygon_.first.empty())
            {
                return new world_type(edge_lengths, matrix_sizes_);
            }
            else
            {
                return new world_type(edge_lengths, matrix_sizes_,
                        this->polygon_.first, this->polygon_.second);
            }
        }
    }

    virtual simulator_type* create_simulator(
        const boost::shared_ptr<world_type>& w, const boost::shared_ptr<Model>& m) const
    {
        return new simulator_type(w, m, bd_dt_factor_, bd_reaction_length_factor_);
    }

  protected:

    boost::shared_ptr<RandomNumberGenerator> rng_;
    std::pair<std::string, STLFormat>        polygon_;
    Integer3 matrix_sizes_;
    Real     bd_dt_factor_;
    Real     bd_reaction_length_factor_;
};

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_SGFRD_FACTORY_HPP */
