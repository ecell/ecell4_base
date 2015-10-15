#ifndef __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP
#define __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP

#include <cstring>
#include <iostream>
#include <sstream>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/lexical_cast.hpp>

#include <hdf5.h>
#include <H5Cpp.h>

#include "types.hpp"
#include "Species.hpp"
#include "Voxel.hpp"


namespace ecell4
{

struct LatticeSpaceHDF5Traits
{
    typedef struct h5_species_struct {
        uint32_t id;
        char serial[32]; // species' serial may exceed the limit
        double radius;
        double D;
        uint32_t location;
    } h5_species_struct;

    typedef struct h5_voxel_struct {
        int lot;
        int serial;
        uint32_t sid;
        int64_t coord;
    } h5_voxel_struct;

    static H5::CompType get_species_comp_type()
    {
        H5::CompType h5_species_comp_type(sizeof(h5_species_struct));
        h5_species_comp_type.insertMember(
            std::string("id"), HOFFSET(h5_species_struct, id),
            H5::PredType::STD_I32LE);
        h5_species_comp_type.insertMember(
            std::string("serial"), HOFFSET(h5_species_struct, serial),
            H5::StrType(H5::PredType::C_S1, 32));
        h5_species_comp_type.insertMember(
            std::string("radius"), HOFFSET(h5_species_struct, radius),
            H5::PredType::NATIVE_DOUBLE);
        h5_species_comp_type.insertMember(
            std::string("D"), HOFFSET(h5_species_struct, D),
            H5::PredType::NATIVE_DOUBLE);
        h5_species_comp_type.insertMember(
            std::string("location"), HOFFSET(h5_species_struct, location),
            H5::PredType::STD_I32LE);
        return h5_species_comp_type;
    }

    static H5::CompType get_voxel_comp_type()
    {
        H5::CompType h5_voxel_comp_type(sizeof(h5_voxel_struct));
        h5_voxel_comp_type.insertMember(
            std::string("lot"), HOFFSET(h5_voxel_struct, lot),
            H5::PredType::NATIVE_INT);
        h5_voxel_comp_type.insertMember(
            std::string("serial"), HOFFSET(h5_voxel_struct, serial),
            H5::PredType::NATIVE_INT);
        h5_voxel_comp_type.insertMember(
            std::string("sid"), HOFFSET(h5_voxel_struct, sid),
            H5::PredType::STD_I32LE);
        h5_voxel_comp_type.insertMember(
            std::string("coordinate"), HOFFSET(h5_voxel_struct, coord),
            H5::PredType::STD_I64LE);
        return h5_voxel_comp_type;
    }

    static void update_h5_species(h5_species_struct& h5_species,
        unsigned int sid, const Species& species)
    {
        h5_species.id = sid;
        h5_species.radius = atof(species.get_attribute("radius").c_str());
        h5_species.D = atof(species.get_attribute("D").c_str());
        std::strcpy(h5_species.serial, species.serial().c_str());
        h5_species.location = 0; // Default location is Vacant (sid=0).
    }

    static void update_h5_voxel(h5_voxel_struct& h5_voxel,
        unsigned int sid, const std::pair<ParticleID, Voxel>& id_voxel_pair)
    {
        h5_voxel.lot = id_voxel_pair.first.lot();
        h5_voxel.serial = id_voxel_pair.first.serial();
        h5_voxel.sid = sid;
        h5_voxel.coord = id_voxel_pair.second.coordinate();
    }

};

template<typename Tspace_>
void save_lattice_space(const Tspace_& space, H5::Group* root)
{
    typedef LatticeSpaceHDF5Traits traits_type;
    typedef typename traits_type::h5_species_struct h5_species_struct;
    typedef typename traits_type::h5_voxel_struct h5_voxel_struct;

    typedef std::vector<std::pair<ParticleID, Voxel> >
        voxel_container_type;

    const unsigned int num_voxels(space.num_voxels());
    const std::vector<Species> species(space.list_species());

    boost::scoped_array<h5_voxel_struct>
        h5_voxel_table(new h5_voxel_struct[num_voxels]);
    boost::scoped_array<h5_species_struct>
        h5_species_table(new h5_species_struct[species.size()]);

    typedef utils::get_mapper_mf<Species::serial_type, uint32_t>::type
        species_id_map_type;
    species_id_map_type species_id_map;

    // Convert species and voxels
    unsigned int vidx(0);
    for (unsigned int i(0); i < species.size(); ++i)
    {
        const unsigned int sid(i + 1);
        traits_type::update_h5_species(h5_species_table[i], sid, species[i]);
        species_id_map[species[i].serial()] = sid;

        const voxel_container_type& voxels(space.list_voxels_exact(species[i]));
        for (unsigned int j(0); j < voxels.size(); ++j)
        {
            traits_type::update_h5_voxel(h5_voxel_table[vidx], sid, voxels[j]);
            ++vidx;
        }
    }

    // Update locations
    for (unsigned int i(0); i < species.size(); ++i)
    {
        species_id_map_type::iterator itr(
                species_id_map.find(species[i].get_attribute("location")));
        if (itr != species_id_map.end())
            h5_species_table[i].location = (*itr).second;
    }

    const int RANK = 1;
    hsize_t dim1[] = {num_voxels};
    H5::DataSpace dataspace1(RANK, dim1);
    boost::scoped_ptr<H5::DataSet> dataset1(new H5::DataSet(
        root->createDataSet(
            "voxels", traits_type::get_voxel_comp_type(), dataspace1)));

    hsize_t dim2[] = {species.size()};
    H5::DataSpace dataspace2(RANK, dim2);
    boost::scoped_ptr<H5::DataSet> dataset2(new H5::DataSet(
        root->createDataSet(
            "species", traits_type::get_species_comp_type(), dataspace2)));

    dataset1->write(h5_voxel_table.get(), dataset1->getDataType());
    dataset2->write(h5_species_table.get(), dataset2->getDataType());

    const uint32_t space_type(static_cast<uint32_t>(Space::LATTICE));
    H5::Attribute attr_space_type(
        root->createAttribute(
            "type", H5::PredType::STD_I32LE, H5::DataSpace(H5S_SCALAR)));
    attr_space_type.write(H5::PredType::STD_I32LE, &space_type);

    const double t(space.t());
    H5::Attribute attr_t(
        root->createAttribute(
            "t", H5::PredType::IEEE_F64LE, H5::DataSpace(H5S_SCALAR)));
    attr_t.write(H5::PredType::IEEE_F64LE, &t);

    const double voxel_radius(space.voxel_radius());
    H5::Attribute attr_voxel_radius(
        root->createAttribute(
            "voxel_radius", H5::PredType::IEEE_F64LE, H5::DataSpace(H5S_SCALAR)));
    attr_voxel_radius.write(H5::PredType::IEEE_F64LE, &voxel_radius);

    const uint32_t is_periodic(space.is_periodic()? 1 : 0);
    H5::Attribute attr_is_periodic(
        root->createAttribute(
            "is_periodic", H5::PredType::STD_I32LE, H5::DataSpace(H5S_SCALAR)));
    attr_is_periodic.write(H5::PredType::STD_I32LE, &is_periodic);


    const Real3 edge_lengths(space.edge_lengths());
    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    H5::Attribute attr_lengths(
        root->createAttribute(
            "edge_lengths", lengths_type, H5::DataSpace(H5S_SCALAR)));
    double lengths[] = {edge_lengths[0], edge_lengths[1], edge_lengths[2]};
    attr_lengths.write(lengths_type, lengths);
}

template<typename Tspace_>
void load_lattice_space(const H5::Group& root, Tspace_* space)
{
    typedef LatticeSpaceHDF5Traits traits_type;
    typedef typename traits_type::h5_species_struct h5_species_struct;
    typedef typename traits_type::h5_voxel_struct h5_voxel_struct;

    double t;
    root.openAttribute("t").read(H5::PredType::IEEE_F64LE, &t);
    space->set_t(t);

    double voxel_radius;
    root.openAttribute("voxel_radius").read(H5::PredType::IEEE_F64LE, &voxel_radius);

    Real3 edge_lengths;
    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    root.openAttribute("edge_lengths").read(lengths_type, &edge_lengths);

    uint32_t is_periodic;
    root.openAttribute("is_periodic").read(
        H5::PredType::STD_I32LE, &is_periodic);

    space->reset(edge_lengths, voxel_radius, (is_periodic != 0));

    {
        H5::DataSet species_dset(root.openDataSet("species"));
        const unsigned int num_species(
            species_dset.getSpace().getSimpleExtentNpoints());
        boost::scoped_array<h5_species_struct> h5_species_table(
            new h5_species_struct[num_species]);
        species_dset.read(
            h5_species_table.get(), traits_type::get_species_comp_type());
        species_dset.close();

        H5::DataSet voxel_dset(root.openDataSet("voxels"));
        const unsigned int num_voxels(
            voxel_dset.getSpace().getSimpleExtentNpoints());
        boost::scoped_array<h5_voxel_struct> h5_voxel_table(
            new h5_voxel_struct[num_voxels]);
        voxel_dset.read(
            h5_voxel_table.get(), traits_type::get_voxel_comp_type());
        voxel_dset.close();

        typedef utils::get_mapper_mf<uint32_t, Species>::type
            species_map_type;
        species_map_type species_map;
        for (unsigned int i(0); i < num_species; ++i)
        {
            std::ostringstream radius, D;
            radius << h5_species_table[i].radius;
            D << h5_species_table[i].D;
            species_map[h5_species_table[i].id] = Species(
                    h5_species_table[i].serial,
                    radius.str(), D.str());
        }

        // Update locations
        for (unsigned int i(0); i < num_species; ++i)
        {
            const unsigned int location(h5_species_table[i].location);
            if (location != 0) {
                species_map[h5_species_table[i].id].set_attribute(
                    "location", species_map[location].serial());
            }
        }

        for (unsigned int i(0); i < num_voxels; ++i)
        {
            const Species& species(species_map[h5_voxel_table[i].sid]);
            space->update_voxel(
                ParticleID(std::make_pair(h5_voxel_table[i].lot, h5_voxel_table[i].serial)),
                Voxel(species, h5_voxel_table[i].coord,
                    atof(species.get_attribute("radius").c_str()),
                    atof(species.get_attribute("D").c_str()),
                    species.get_attribute("location")));
        }
    }
}

} // ecell4

#endif /*  __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP */
