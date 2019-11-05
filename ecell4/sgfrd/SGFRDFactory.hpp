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
        : base_type(), rng_(nullptr), polygon_(nullptr),
          polygon_file_("", STLFormat::Ascii),
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

    // -----------------------------------------------------------------------
    // get polygon from a list of triangles

    this_type& polygon(const Real3& el, const std::vector<Triangle>& ts)
    {
        polygon_file_.first = ""; // XXX clear polygon file

        this->polygon_ = boost::make_shared<Polygon>(el, ts);
        return (*this);
    }
    this_type* polygon_ptr(const Real3& el, const std::vector<Triangle>& ts)
    {
        return std::addressof(this->polygon(el, ts));
    }

    // -----------------------------------------------------------------------
    // read polygon from .STL file

    this_type& polygon(const std::string& fname, const STLFormat fmt)
    {
        polygon_ = nullptr; // XXX clear polygon structure

        this->polygon_file_ = std::make_pair(fname, fmt);
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
            if(this->polygon_)
            {
                return new world_type(edge_lengths, matrix_sizes_, rng_,
                                      this->polygon_);
            }
            else if(!this->polygon_file_.first.empty())
            {
                return new world_type(edge_lengths, matrix_sizes_, rng_,
                        this->polygon_file_.first, this->polygon_file_.second);
            }
            else
            {
                return new world_type(edge_lengths, matrix_sizes_, rng_);
            }
        }
        else
        {
            if(this->polygon_)
            {
                return new world_type(edge_lengths, matrix_sizes_,
                                      this->polygon_);
            }
            else if(!this->polygon_file_.first.empty())
            {
                return new world_type(edge_lengths, matrix_sizes_,
                        this->polygon_file_.first, this->polygon_file_.second);
            }
            else
            {
                return new world_type(edge_lengths, matrix_sizes_);
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
    boost::shared_ptr<Polygon>               polygon_;
    std::pair<std::string, STLFormat>        polygon_file_;
    Integer3 matrix_sizes_;
    Real     bd_dt_factor_;
    Real     bd_reaction_length_factor_;
};

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_SGFRD_FACTORY_HPP */
