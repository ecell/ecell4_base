#ifndef __ECELL4_SPATIOCYTE_COORDINATE_LOGGER
#define __ECELL4_SPATIOCYTE_COORDINATE_LOGGER

#include <fstream>
#include "SpatiocyteWorld.hpp"


namespace ecell4
{

namespace spatiocyte
{

class CoordinateLogger
{
public:

    typedef SpatiocyteWorld space_type;

protected:

    typedef std::vector<Species> species_container_type;

public:

    CoordinateLogger(
        boost::shared_ptr<space_type> space,
        const std::string& filename = "CoordinateLog.csv")
        : space_(space), filename_(filename)
    {
        ;
    }

    std::string filename() const
    {
        return filename_;
    }

    void set_filename(const std::string& filename)
    {
        filename_ = filename;
    }

    void initialize()
    {
        std::ofstream fout(filename().c_str(), std::ios::trunc);

        const Real& voxel_radius((*space_).voxel_radius());
        const Real voxel_size(2 * voxel_radius);
        const Position3& edge_lengths((*space_).edge_lengths());

        fout << "log interval=" // << dt_
             << ",world width=" << edge_lengths[2] / voxel_size
             << ",world height=" << edge_lengths[1] / voxel_size
             << ",world length=" << edge_lengths[0] / voxel_size
             << ",voxel radius=" << voxel_radius;

        for (species_container_type::const_iterator i(species_.begin());
             i != species_.end(); ++i)
        {
            fout << "," << (*i).serial()
                 << "=" << voxel_radius; // (*i).radius()
        }

        fout << std::endl;
        fout.flush();
        fout.close();
    }

    void log()
    {
        const Real voxel_size(2 * (*space_).voxel_radius());
        std::ofstream fout(filename().c_str(), std::ios::out | std::ios::app);

        for (species_container_type::const_iterator i(species_.begin());
             i != species_.end(); ++i)
        {
            fout << (*space_).t();

            std::vector<std::pair<ParticleID, Particle> >
                particles((*space_).list_particles(*i));
            for (std::vector<std::pair<ParticleID, Particle> >::const_iterator
                     i(particles.begin()); i != particles.end(); ++i)
            {
                const Position3 pos(divide((*i).second.position(), voxel_size));
                fout << "," << pos[0] << "," << pos[1] << "," << pos[2];
            }

            fout << std::endl;
        }

        fout.flush();
        fout.close();
    }

    bool has_species(const Species& sp)
    {
        for (species_container_type::const_iterator i(species_.begin());
             i != species_.end(); ++i)
        {
            if (*i == sp)
            {
                return true;
            }
        }
        return false;
    }

    void add_species(const Species& sp)
    {
        if (!has_species(sp))
        {
            species_.push_back(sp);
        }
    }

protected:

    boost::shared_ptr<space_type> space_;
    std::string filename_;

    species_container_type species_;
};

} // spatiocyte

} // ecell4

#endif /* __ECELL4_SPATIOCYTE_COORDINATE_LOGGER */
