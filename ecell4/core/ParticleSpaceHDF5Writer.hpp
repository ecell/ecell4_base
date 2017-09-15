#ifndef ECELL4_PARTICLE_SPACE_HDF5_WRITER_HPP
#define ECELL4_PARTICLE_SPACE_HDF5_WRITER_HPP

#include <cstring>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

#include <hdf5.h>
#include <H5Cpp.h>

#include "types.hpp"
#include "get_mapper_mf.hpp"
#include "Species.hpp"
#include "Particle.hpp"
#include "Space.hpp"


namespace ecell4
{

struct ParticleSpaceHDF5Traits
{
    typedef struct h5_species_struct {
        uint32_t id;
        char serial[32]; // species' serial may exceed the limit
    } h5_species_struct;

    typedef struct h5_particle_struct {
        int lot;
        int serial;
        uint32_t sid;
        double posx;
        double posy;
        double posz;
        double radius;
        double D;
    } h5_particle_struct;

    static H5::CompType get_particle_comp_type()
    {
        H5::CompType h5_particle_comp_type(sizeof(h5_particle_struct));
#define INSERT_MEMBER(member, type) \
        H5Tinsert(h5_particle_comp_type.getId(), #member,\
                HOFFSET(h5_particle_struct, member), type.getId())
        INSERT_MEMBER(lot, H5::PredType::NATIVE_INT);
        INSERT_MEMBER(serial, H5::PredType::NATIVE_INT);
        INSERT_MEMBER(sid, H5::PredType::STD_I32LE);
        INSERT_MEMBER(posx, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(posy, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(posz, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(radius, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(D, H5::PredType::NATIVE_DOUBLE);
#undef INSERT_MEMBER
        // h5_particle_comp_type.insertMember(
        //     std::string("lot"), HOFFSET(h5_particle_struct, lot),
        //     H5::PredType::NATIVE_INT);
        // h5_particle_comp_type.insertMember(
        //     std::string("serial"), HOFFSET(h5_particle_struct, serial),
        //     H5::PredType::NATIVE_INT);
        // h5_particle_comp_type.insertMember(
        //     std::string("sid"), HOFFSET(h5_particle_struct, sid),
        //     H5::PredType::STD_I32LE);
        // h5_particle_comp_type.insertMember(
        //     std::string("posx"), HOFFSET(h5_particle_struct, posx),
        //     H5::PredType::NATIVE_DOUBLE);
        // h5_particle_comp_type.insertMember(
        //     std::string("posy"), HOFFSET(h5_particle_struct, posy),
        //     H5::PredType::NATIVE_DOUBLE);
        // h5_particle_comp_type.insertMember(
        //     std::string("posz"), HOFFSET(h5_particle_struct, posz),
        //     H5::PredType::NATIVE_DOUBLE);
        // h5_particle_comp_type.insertMember(
        //     std::string("radius"), HOFFSET(h5_particle_struct, radius),
        //     H5::PredType::NATIVE_DOUBLE);
        // h5_particle_comp_type.insertMember(
        //     std::string("D"), HOFFSET(h5_particle_struct, D),
        //     H5::PredType::NATIVE_DOUBLE);
        return h5_particle_comp_type;
    }

    static H5::CompType get_species_comp_type()
    {
        H5::CompType h5_species_comp_type(sizeof(h5_species_struct));
#define INSERT_MEMBER(member, type) \
        H5Tinsert(h5_species_comp_type.getId(), #member,\
                HOFFSET(h5_species_struct, member), type.getId())
        INSERT_MEMBER(id, H5::PredType::STD_I32LE);
        INSERT_MEMBER(serial, H5::StrType(H5::PredType::C_S1, 32));
#undef INSERT_MEMBER
        // h5_species_comp_type.insertMember(
        //     std::string("id"), HOFFSET(h5_species_struct, id),
        //     H5::PredType::STD_I32LE);
        // h5_species_comp_type.insertMember(
        //     std::string("serial"), HOFFSET(h5_species_struct, serial),
        //     H5::StrType(H5::PredType::C_S1, 32));
        return h5_species_comp_type;
    }
};

template<typename Tspace_>
void save_particle_space(const Tspace_& space, H5::Group* root)
{
    typedef ParticleSpaceHDF5Traits traits_type;
    typedef typename traits_type::h5_species_struct h5_species_struct;
    typedef typename traits_type::h5_particle_struct h5_particle_struct;

    typedef std::vector<std::pair<ParticleID, Particle> >
        particle_container_type;
    const particle_container_type& particles(space.list_particles());
    const unsigned int num_particles(particles.size());

    std::vector<Species> species;
    typedef utils::get_mapper_mf<Species::serial_type, unsigned int>::type
        species_id_map_type;
    species_id_map_type species_id_map;

    boost::scoped_array<h5_particle_struct>
        h5_particle_table(new h5_particle_struct[num_particles]);
    for (unsigned int i(0); i < num_particles; ++i)
    {
        species_id_map_type::const_iterator
            it(species_id_map.find(particles[i].second.species_serial()));
        if (it == species_id_map.end())
        {
            species.push_back(particles[i].second.species());
            it = species_id_map.insert(
                std::make_pair(particles[i].second.species_serial(),
                               species.size())).first;
        }

        h5_particle_table[i].lot = particles[i].first.lot();
        h5_particle_table[i].serial = particles[i].first.serial();
        h5_particle_table[i].sid = (*it).second;
        h5_particle_table[i].posx = particles[i].second.position()[0];
        h5_particle_table[i].posy = particles[i].second.position()[1];
        h5_particle_table[i].posz = particles[i].second.position()[2];
        h5_particle_table[i].radius = particles[i].second.radius();
        h5_particle_table[i].D = particles[i].second.D();
    }

    boost::scoped_array<h5_species_struct>
        h5_species_table(new h5_species_struct[species.size()]);
    for (unsigned int i(0); i < species.size(); ++i)
    {
        h5_species_table[i].id = i + 1;
        std::strcpy(h5_species_table[i].serial,
                    species[i].serial().c_str());
    }

    const int RANK = 1;
    hsize_t dim1[] = {num_particles};
    H5::DataSpace dataspace1(RANK, dim1);
    boost::scoped_ptr<H5::DataSet> dataset1(new H5::DataSet(
        root->createDataSet(
            "particles", traits_type::get_particle_comp_type(), dataspace1)));

    hsize_t dim2[] = {species.size()};
    H5::DataSpace dataspace2(RANK, dim2);
    boost::scoped_ptr<H5::DataSet> dataset2(new H5::DataSet(
        root->createDataSet(
            "species", traits_type::get_species_comp_type(), dataspace2)));

    dataset1->write(h5_particle_table.get(), dataset1->getDataType());
    dataset2->write(h5_species_table.get(), dataset2->getDataType());

    const uint32_t space_type = static_cast<uint32_t>(Space::PARTICLE);
    H5::Attribute attr_space_type(
        root->createAttribute(
            "type", H5::PredType::STD_I32LE, H5::DataSpace(H5S_SCALAR)));
    attr_space_type.write(H5::PredType::STD_I32LE, &space_type);

    const double t = space.t();
    H5::Attribute attr_t(
        root->createAttribute(
            "t", H5::PredType::IEEE_F64LE, H5::DataSpace(H5S_SCALAR)));
    attr_t.write(H5::PredType::IEEE_F64LE, &t);

    const Real3 edge_lengths = space.edge_lengths();
    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    H5::Attribute attr_lengths(
        root->createAttribute(
            "edge_lengths", lengths_type, H5::DataSpace(H5S_SCALAR)));
    double lengths[] = {edge_lengths[0], edge_lengths[1], edge_lengths[2]};
    attr_lengths.write(lengths_type, lengths);
}

template<typename Tspace_>
void load_particle_space(const H5::Group& root, Tspace_* space)
{
    typedef ParticleSpaceHDF5Traits traits_type;
    typedef typename traits_type::h5_species_struct h5_species_struct;
    typedef typename traits_type::h5_particle_struct h5_particle_struct;

    Real3 edge_lengths;
    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    root.openAttribute("edge_lengths").read(lengths_type, &edge_lengths);
    space->reset(edge_lengths);

    double t;
    root.openAttribute("t").read(H5::PredType::IEEE_F64LE, &t);
    space->set_t(t);

    {
        H5::DataSet species_dset(root.openDataSet("species"));
        const unsigned int num_species(
            species_dset.getSpace().getSimpleExtentNpoints());
        boost::scoped_array<h5_species_struct> h5_species_table(
            new h5_species_struct[num_species]);
        species_dset.read(
            h5_species_table.get(), traits_type::get_species_comp_type());
        species_dset.close();

        H5::DataSet particle_dset(root.openDataSet("particles"));
        const unsigned int num_particles(
            particle_dset.getSpace().getSimpleExtentNpoints());
        boost::scoped_array<h5_particle_struct> h5_particle_table(
            new h5_particle_struct[num_particles]);
        particle_dset.read(
            h5_particle_table.get(), traits_type::get_particle_comp_type());
        particle_dset.close();

        typedef utils::get_mapper_mf<unsigned int, Species::serial_type>::type
            species_id_map_type;
        species_id_map_type species_id_map;
        for (unsigned int i(0); i < num_species; ++i)
        {
            species_id_map[h5_species_table[i].id] = h5_species_table[i].serial;
        }

        for (unsigned int i(0); i < num_particles; ++i)
        {
            space->update_particle(ParticleID(std::make_pair(h5_particle_table[i].lot, h5_particle_table[i].serial)), Particle(Species(species_id_map[h5_particle_table[i].sid]), Real3(h5_particle_table[i].posx, h5_particle_table[i].posy, h5_particle_table[i].posz), h5_particle_table[i].radius, h5_particle_table[i].D));
        }

        // boost::scoped_array<h5_particle_struct>
        //     h5_particle_table(new h5_particle_struct[num_particles]);
        // for (unsigned int i(0); i < num_particles; ++i)
        // {
        //     species_id_map_type::const_iterator
        //         it(species_id_map.find(particles[i].second.species_serial()));
        //     if (it == species_id_map.end())
        //     {
        //         species.push_back(particles[i].second.species());
        //         it = species_id_map.insert(
        //             std::make_pair(particles[i].second.species_serial(),
        //                            species.size())).first;
        //     }

        //     h5_particle_table[i].lot = particles[i].first.lot();
        //     h5_particle_table[i].serial = particles[i].first.serial();
        //     h5_particle_table[i].sid = (*it).second;
        //     h5_particle_table[i].posx = particles[i].second.position()[0];
        //     h5_particle_table[i].posy = particles[i].second.position()[1];
        //     h5_particle_table[i].posz = particles[i].second.position()[2];
        //     h5_particle_table[i].radius = particles[i].second.radius();
        //     h5_particle_table[i].D = particles[i].second.D();
        // }

        // boost::scoped_array<h5_species_struct>
        //     h5_species_table(new h5_species_struct[species.size()]);
        // for (unsigned int i(0); i < species.size(); ++i)
        // {
        //     h5_species_table[i].id = i + 1;
        //     std::strcpy(h5_species_table[i].serial,
        //                 species[i].serial().c_str());
        // }
    }
}

} // ecell4

#endif /*  ECELL4_PARTICLE_SPACE_HDF5_WRITER_HPP */
