#ifndef __ECELL4_SPATIOCYTE_SPATIOCYTE_VISUALIZATION_LOGGER
#define __ECELL4_SPATIOCYTE_SPATIOCYTE_VISUALIZATION_LOGGER

namespace ecell4
{

namespace spatiocyte
{

class SpatiocyteVisualizationLogger
{
public:

    typedef SpatiocyteWorld space_type;

protected:

    typedef std::vector<Species> species_container_type;
    // typedef space_type::particle_id_pair_ector_type
    // particle_id_pair_vctor_type;
    typedef std::vector<std::pair<ParticleID, Particle> >
    particle_id_pair_vector_type;

public:

    SpatiocyteVisualizationLogger(
        boost::shared_ptr<space_type> space,
        const std::string& filename = "VisualLog.dat")
        : space_(space), filename_(filename),
          log_polymer_(true), log_marker_(UINT_MAX)
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
        for (species_container_type::const_iterator i(species_.begin());
             i != species_.end(); ++i)
        {
            const space_type::particle_info_type
                info((*space_).get_particle_info(*i));
            if (!info.is_lattice)
            {
                offlattice_species_.push_back(*i);
            }
            else
            {
                lattice_species_.push_back(*i);
                if (info.is_polymer && log_polymer_)
                {
                    polymer_indices_.push_back(lattice_species_.size());
                    polymer_species_.push_back(*i);
                }
            }
        }

        log_settings();
        log_vacants();
    }

    void log_settings()
    {
        std::ofstream fout(
            filename().c_str(), std::ios::binary | std::ios::trunc);

        unsigned int lattice_type(HCP_LATTICE);
        fout.write((char*)(&lattice_type), sizeof(lattice_type));

        unsigned int mean_count(0);
        fout.write((char*)(&mean_count), sizeof(mean_count));

        unsigned int start_coord(0);
        fout.write((char*)(&start_coord), sizeof(start_coord));

        const boost::array<Integer, 3> sizes((*space_).lattice_size());
        unsigned int
            col_size(sizes[0]), layer_size(sizes[1]), row_size(sizes[2]);
        fout.write((char*)(&row_size), sizeof(row_size));
        fout.write((char*)(&layer_size), sizeof(layer_size));
        fout.write((char*)(&col_size), sizeof(col_size));

        Position3 edge_lengths(
            divide((*space_).edge_lengths(), 2 * (*space_).voxel_radius()));
        ::Point center = {
            edge_lengths[0] * 0.5, edge_lengths[0] * 0.5, edge_lengths[0] * 0.5};

        double real_row_size(center.z * 2), real_layer_size(center.y * 2),
            real_col_size(center.x * 2);
        fout.write((char*)(&real_row_size), sizeof(real_row_size));
        fout.write((char*)(&real_layer_size), sizeof(real_layer_size));
        fout.write((char*)(&real_col_size), sizeof(real_col_size));

        unsigned int lattice_species_size(lattice_species_.size()),
            polymer_species_size(polymer_species_.size()),
            offlattice_species_size(offlattice_species_.size());
        unsigned int reserved_size(0);
        fout.write(
            (char*)(&lattice_species_size), sizeof(lattice_species_size));
        fout.write(
            (char*)(&polymer_species_size), sizeof(polymer_species_size));
        fout.write((char*)(&reserved_size), sizeof(reserved_size));
        fout.write(
            (char*)(&offlattice_species_size), sizeof(offlattice_species_size));

        fout.write((char*)(&log_marker_), sizeof(log_marker_));

        double voxel_radius((*space_).voxel_radius());
        fout.write((char*)(&voxel_radius), sizeof(voxel_radius));

        for (species_container_type::const_iterator i(lattice_species_.begin());
             i != lattice_species_.end(); ++i)
        {
            // std::string fullid(
            //     theLatticeSpecies[i]->getVariable()->getFullID().asString());
            std::string fullid((*i).serial());
            unsigned int string_size(fullid.size());
            fout.write((char*)(&string_size), sizeof(string_size));
            fout.write(fullid.c_str(), string_size);

            // double radius(theLatticeSpecies[i]->getRadius());
            double radius((*space_).voxel_radius());
            fout.write((char*)(&radius), sizeof(radius));
        }

        {
            species_container_type::const_iterator i(polymer_species_.begin());
            std::vector<species_container_type::size_type>::const_iterator
                j(polymer_indices_.begin());
            while (i != polymer_species_.end())
            {
                unsigned int polymer_index(*j);
                fout.write((char*)(&polymer_index), sizeof(polymer_index));

                // double radius(thePolymerSpecies[i]->getRadius());
                double radius((*space_).voxel_radius());
                fout.write((char*)(&radius), sizeof(radius));

                ++i;
                ++j;
            }
        }

        for (species_container_type::const_iterator
                 i(offlattice_species_.begin());
             i != offlattice_species_.end(); ++i)
        {
            // std::string fullid(
            //     theLatticeSpecies[i]->getVariable()->getFullID().asString());
            std::string fullid((*i).serial());
            unsigned int string_size(fullid.size());
            fout.write((char*)(&string_size), sizeof(string_size));
            fout.write(fullid.c_str(), string_size);

            // double radius(theLatticeSpecies[i]->getRadius());
            double radius((*space_).voxel_radius());
            fout.write((char*)(&radius), sizeof(radius));
        }

        fout.flush();
        fout.close();
    }

    void log_vacants()
    {
        const Real voxel_size((*space_).voxel_radius() * 2);

        std::ofstream fout(
            filename().c_str(), std::ios::binary | std::ios::app);

        double t((*space_).t());
        fout.write((char*)(&t), sizeof(t));

        for (species_container_type::const_iterator i(lattice_species_.begin());
             i != lattice_species_.end(); ++i)
        {
            if ((*space_).get_particle_info(*i).is_vacant)
            {
                // The species index in the process:
                species_container_type::const_iterator
                    iter(lattice_species_.begin());
                unsigned int idx(distance(iter, i));
                fout.write((char*)(&idx), sizeof(idx));

                std::vector<unsigned int> coordinates((*space_).coordinates(*i));

                // The species molecule size:
                unsigned int vacant_size(coordinates.size());
                fout.write((char*)(&vacant_size), sizeof(vacant_size));

                for (std::vector<unsigned int>::const_iterator
                         j(coordinates.begin()); j != coordinates.end(); ++j)
                {
                    unsigned int coord(*j);
                    fout.write((char*)(&coord), sizeof(coord));
                }
            }
        }

        fout.write((char*)(&log_marker_), sizeof(log_marker_));

        for (species_container_type::const_iterator
                 i(offlattice_species_.begin());
             i != offlattice_species_.end(); ++i)
        {
            if ((*space_).get_particle_info(*i).is_vacant)
            {
                // The species index in the process:
                species_container_type::const_iterator
                    iter(offlattice_species_.begin());
                unsigned int idx(std::distance(iter, i));
                fout.write((char*)(&idx), sizeof(idx));

                particle_id_pair_vector_type particles(
                    (*space_).list_particles(*i));

                // The species molecule size:
                unsigned int vacant_size(particles.size());
                fout.write((char*)(&vacant_size), sizeof(vacant_size));

                for (particle_id_pair_vector_type::const_iterator
                         j(particles.begin()); j != particles.end(); ++j)
                {
                    const Position3 pos(divide((*j).second.position(), voxel_size));
                    Point point = {pos[0], pos[1], pos[2]};
                    fout.write((char*)(&point), sizeof(point));
                }
            }
        }

        fout.write((char*)(&log_marker_), sizeof(log_marker_));

        fout.flush();
        fout.close();
    }

    void log()
    {
        const Real voxel_size((*space_).voxel_radius() * 2);

        std::ofstream fout(
            filename().c_str(), std::ios::binary | std::ios::app);

        double t((*space_).t());
        fout.write((char*)(&t), sizeof(t));

        for (species_container_type::const_iterator i(lattice_species_.begin());
             i != lattice_species_.end(); ++i)
        {
            if ((*space_).get_particle_info(*i).is_vacant)
            {
                continue;
            }

            species_container_type::const_iterator
                iter(lattice_species_.begin());
            int idx(std::distance(iter, i)); // why not unsigned here?
            fout.write((char*)(&idx), sizeof(idx));

            std::vector<unsigned int> coordinates((*space_).coordinates(*i));

            // The species molecule size
            int species_size(coordinates.size()); // why not unsigned here?
            fout.write((char*)(&species_size), sizeof(species_size));

            for (std::vector<unsigned int>::const_iterator
                     j(coordinates.begin()); j != coordinates.end(); ++j)
            {
                unsigned int coord(*j);
                fout.write((char*)(&coord), sizeof(coord));
            }
        }

        // for (species_container_type::const_iterator i(polymer_species_.begin());
        //      i != polymer_species_.end(); ++i)
        // {
        //     log_source_molecules(*i);
        // }

        // for (species_container_type::const_iterator i(polymer_species_.begin());
        //      i != polymer_species_.end(); ++i)
        // {
        //     log_target_molecules(*i);
        // }

        // for (species_container_type::const_iterator i(polymer_species_.begin());
        //      i != polymer_species_.end(); ++i)
        // {
        //     log_shared_molecules(*i);
        // }

        // for (species_container_type::const_iterator i(reserved_species_.begin());
        //      i != reserved_species_.end(); ++i)
        // {
        //     log_reserved_molecules(*i);
        // }

        fout.write((char*)(&log_marker_), sizeof(log_marker_));

        // for (species_container_type::const_iterator i(polymer_species_.begin());
        //      i != polymer_species_.end(); ++i)
        // {
        //     log_polymers(*i);
        // }

        for (species_container_type::const_iterator
                 i(offlattice_species_.begin());
             i != offlattice_species_.end(); ++i)
        {
            const space_type::particle_info_type info(
                (*space_).get_particle_info(*i));
            if (info.is_vacant
                && !info.is_diffusive_vacant
                && !info.is_reactive_vacant)
            {
                continue;
            }

            // The species index in the process:
            species_container_type::const_iterator
                iter(offlattice_species_.begin());
            int idx(std::distance(iter, i)); // why not unsigned here?
            fout.write((char*)(&idx), sizeof(idx));

            particle_id_pair_vector_type particles(
                (*space_).list_particles(*i));

            // The species molecule size:
            int vacant_size(particles.size()); // why not unsigned here?
            fout.write((char*)(&vacant_size), sizeof(vacant_size));

            for (particle_id_pair_vector_type::const_iterator
                     j(particles.begin()); j != particles.end(); ++j)
            {
                const Position3 pos(divide((*j).second.position(), voxel_size));
                Point point = {pos[0], pos[1], pos[2]};
                fout.write((char*)(&point), sizeof(point));
            }
        }

        fout.write((char*)(&log_marker_), sizeof(log_marker_));
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
    bool log_polymer_;
    unsigned int log_marker_;

    species_container_type species_;

    species_container_type lattice_species_;
    species_container_type offlattice_species_;
    species_container_type polymer_species_;
    std::vector<species_container_type::size_type> polymer_indices_;
};

} // spatiocyte

} // ecell4

#endif /* __ECELL4_SPATIOCYTE_SPATIOCYTE_VISUALIZATION_LOGGER */
