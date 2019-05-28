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

    SGFRDFactory(const Integer3& matrix_sizes = default_matrix_sizes(),
                 Real bd_dt_factor = default_bd_dt_factor(),
                 Real bd_reaction_length_factor = default_bd_reaction_length_factor())
        : base_type(), rng_(nullptr), polygon_(nullptr), matrix_sizes_(matrix_sizes),
          bd_dt_factor_(bd_dt_factor), bd_reaction_length_factor_(bd_reaction_length_factor)
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

    this_type& polygon(const boost::shared_ptr<Polygon>& poly)
    {
        this->polygon_ = poly;
        return (*this);
    }
    this_type* polygon_ptr(const boost::shared_ptr<Polygon>& poly)
    {
        return std::addressof(this->polygon(poly));
    }

  protected:

    virtual world_type* create_world(const Real3& edge_lengths) const
    {
        if (!polygon_)
        {
            // generate default polygon: XY-plane at the middle of the Z axis.

            const std::size_t x_size = matrix_sizes_[0];
            const std::size_t y_size = matrix_sizes_[1];

            assert(x_size != 0);
            assert(y_size != 0);

            const Real dx = edge_lengths[0] / x_size;
            const Real dy = edge_lengths[1] / y_size;
            const Real z  = edge_lengths[2] / 2.0; // at the middle of Z-axis

            std::vector<Triangle> ts;
            ts.reserve(x_size * y_size * 2);

            for(std::size_t yi = 0; yi < y_size; ++yi)
            {
                for(std::size_t xi = 0; xi < x_size; ++xi)
                {
                    //   4___ 3
                    // y  | /|  upper left  = {1, 3, 4}
                    // ^  |/_|  lower right = {1, 2, 3}
                    // | 1    2
                    // |
                    // --> x

                    const Real3 v1(dx *  xi   , dy *  yi   , z);
                    const Real3 v2(dx * (xi+1), dy *  yi   , z);
                    const Real3 v3(dx * (xi+1), dy * (yi+1), z);
                    const Real3 v4(dx *  xi   , dy * (yi+1), z);

                    ts.push_back(Triangle(v1, v3, v4));
                    ts.push_back(Triangle(v1, v2, v3));
                }
            }
            this->polygon_ = boost::make_shared<Polygon>(edge_lengths, ts);
        }

        if (rng_)
        {
            return new world_type(edge_lengths, matrix_sizes_, polygon_, rng_);
        }
        else
        {
            return new world_type(edge_lengths, matrix_sizes_, polygon_);
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
    Integer3 matrix_sizes_;
    Real     bd_dt_factor_;
    Real     bd_reaction_length_factor_;
};

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_SGFRD_FACTORY_HPP */
