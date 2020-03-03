#ifndef ECELL4_BOUNDARY_CONDITION_HPP
#define ECELL4_BOUNDARY_CONDITION_HPP
#include "types.hpp"
#include "exceptions.hpp"
#include "Real3.hpp"
#include <cmath>

namespace ecell4
{

class Boundary
{
public:
    virtual ~Boundary() = default;
    virtual Real3 periodic_transpose(const Real3&, const Real3&) const noexcept = 0;
    virtual Real3 apply_boundary    (const Real3&) const noexcept = 0;
    virtual Real3 const& edge_lengths() const noexcept = 0;
    virtual void update(const Real3&) noexcept = 0;
};

class UnlimitedBoundary final : public Boundary
{
public:

    UnlimitedBoundary() noexcept :
        edge_lengths_(std::numeric_limits<Real>::infinity(),
                      std::numeric_limits<Real>::infinity(),
                      std::numeric_limits<Real>::infinity())
    {}
    explicit UnlimitedBoundary(const Real3& /*ignored*/) noexcept :
        edge_lengths_(std::numeric_limits<Real>::infinity(),
                      std::numeric_limits<Real>::infinity(),
                      std::numeric_limits<Real>::infinity())
    {}
    ~UnlimitedBoundary() override = default;

    Real3 periodic_transpose(const Real3& pos, const Real3&) const noexcept override
    {
        return pos;
    }
    Real3 apply_boundary(const Real3& pos) const noexcept override
    {
        return pos;
    }
    Real3 const& edge_lengths() const noexcept override
    {
        return this->edge_lengths_;
    }
    void update(const Real3&) noexcept override
    {
        return ;
    }

private:
    Real3 edge_lengths_;
};

class PeriodicBoundary final : public Boundary
{
public:

    PeriodicBoundary() noexcept
        : edge_lengths_(0.0, 0.0, 0.0), half_widths_(0.0, 0.0, 0.0)
    {}
    explicit PeriodicBoundary(const Real3& edges) noexcept
        : edge_lengths_(edges), half_widths_(edges * 0.5)
    {}
    ~PeriodicBoundary() override = default;

    Real3 periodic_transpose(const Real3& pos1, const Real3& pos2) const noexcept override
    {
        Real3 retval(pos1);
        for(Real3::size_type dim(0); dim < 3; ++dim)
        {
            const Real edge_length(edge_lengths_[dim]);
            const Real diff(pos2[dim] - pos1[dim]), half(half_widths_[dim]);

            if (diff > half)
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
    Real3 apply_boundary(const Real3& pos) const noexcept override
    {
        return modulo(pos, this->edge_lengths_);
    }
    Real3 const& edge_lengths() const noexcept override
    {
        return edge_lengths_;
    }
    void update(const Real3& edges) noexcept override
    {
        this->edge_lengths_ = edges;
        this->half_widths_  = edges * 0.5;
        return ;
    }

private:
    Real3 edge_lengths_;
    Real3 half_widths_;
};

} // ecell4
#endif//ECELL4_BOUNDARY_CONDITION_HPP
