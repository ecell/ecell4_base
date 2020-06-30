#ifndef ECELL4_EGFRD_VOLUME_CLEARER_HPP
#define ECELL4_EGFRD_VOLUME_CLEARER_HPP
#include <ecell4/core/Identifier.hpp>
#include <ecell4/core/Sphere.hpp>

namespace ecell4
{
namespace egfrd
{

class VolumeClearer
{
public:
    typedef ecell4::Sphere     particle_shape_type;
    typedef ecell4::ParticleID particle_id_type;

public:
    virtual ~VolumeClearer() {}

    virtual bool operator()(particle_shape_type const& shape, particle_id_type const& ignore) = 0;

    virtual bool operator()(particle_shape_type const& shape, particle_id_type const& ignore0, particle_id_type const& ignore1) = 0;
};

} // egfrd
} // ecell4
#endif /* VOLUME_CLEARER_HPP */
