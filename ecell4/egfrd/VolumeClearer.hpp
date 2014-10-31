#ifndef VOLUME_CLEARER_HPP
#define VOLUME_CLEARER_HPP

template<typename Tps_, typename Tpid_>
class VolumeClearer
{
public:
    typedef Tps_ particle_shape_type;
    typedef Tpid_ particle_id_type;

public:
    virtual ~VolumeClearer() {}

    virtual bool operator()(particle_shape_type const& shape, particle_id_type const& ignore) = 0;

    virtual bool operator()(particle_shape_type const& shape, particle_id_type const& ignore0, particle_id_type const& ignore1) = 0;
};

#endif /* VOLUME_CLEARER_HPP */
