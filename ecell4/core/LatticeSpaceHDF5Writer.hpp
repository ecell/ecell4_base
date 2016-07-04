#ifndef __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP
#define __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP

#include <cstring>
#include <iostream>
#include <sstream>
#include <map>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/lexical_cast.hpp>

#include <hdf5.h>
#include <H5Cpp.h>

#include "types.hpp"
#include "Species.hpp"
#include "Voxel.hpp"
#include "MolecularTypeBase.hpp"
#include "StructureType.hpp"
#include "VacantType.hpp"


namespace ecell4
{

struct LatticeSpaceHDF5Traits
{
    struct h5_species_struct {
        double radius;
        double D;
        char location[32];
        uint32_t is_structure;
        uint32_t dimension;
    };

    struct h5_voxel_struct {
        uint32_t lot;
        uint32_t serial;
        uint64_t coordinate;
    };

    static H5::CompType get_property_comp()
    {
        H5::CompType property_comp_type(sizeof(h5_species_struct));
// #define INSERT_MEMBER(member, type) \
//         property_comp_type.insertMember(#member,\
//                 HOFFSET(h5_species_struct, member), type)
#define INSERT_MEMBER(member, type) \
        H5Tinsert(property_comp_type.getId(), #member,\
                HOFFSET(h5_species_struct, member), type.getId())
        INSERT_MEMBER(radius, H5::PredType::IEEE_F64LE);
        INSERT_MEMBER(D, H5::PredType::IEEE_F64LE);
        INSERT_MEMBER(location, H5::StrType(H5::PredType::C_S1, 32));
        INSERT_MEMBER(is_structure, H5::PredType::STD_I32LE);
        INSERT_MEMBER(dimension, H5::PredType::STD_I32LE);
#undef INSERT_MEMBER
        return property_comp_type;
    }

    static H5::CompType get_voxel_comp()
    {
        H5::CompType voxel_comp_type(sizeof(h5_voxel_struct));
// #define INSERT_MEMBER(member, type) \
//         voxel_comp_type.insertMember(std::string(#member),\
//             HOFFSET(h5_voxel_struct, member), type)
#define INSERT_MEMBER(member, type) \
        H5Tinsert(voxel_comp_type.getId(), #member,\
                HOFFSET(h5_voxel_struct, member), type.getId())
        INSERT_MEMBER(lot, H5::PredType::NATIVE_INT);
        INSERT_MEMBER(serial, H5::PredType::NATIVE_INT);
        INSERT_MEMBER(coordinate, H5::PredType::STD_I64LE);
#undef INSERT_MEMBER
        return voxel_comp_type;
    }

    static void save_molecular_type(const VoxelPool* mtb,
            std::vector<std::pair<ParticleID, Voxel> > voxels, H5::Group* group)
    {
        const Species species(mtb->species());
        boost::scoped_ptr<H5::Group> mtgroup(
                new H5::Group(group->createGroup(species.serial().c_str())));

        h5_species_struct property;
        property.radius = mtb->radius();
        property.D = mtb->D();
        const VoxelPool* loc(mtb->location());
        if (loc->is_vacant())
            std::strcpy(property.location, "");
        else
            std::strcpy(property.location, loc->species().serial().c_str());
        property.is_structure = mtb->is_structure() ? 1 : 0;
        property.dimension = mtb->get_dimension();

        H5::CompType property_comp_type(get_property_comp());
        mtgroup->createAttribute("property", property_comp_type,
                H5::DataSpace(H5S_SCALAR)).write(property_comp_type, &property);

        // Save voxels
        const Integer num_voxels(voxels.size());
        std::size_t vidx(0);
        boost::scoped_array<h5_voxel_struct> h5_voxel_array(new h5_voxel_struct[num_voxels]);
        for (std::vector<std::pair<ParticleID, Voxel> >::const_iterator itr(voxels.begin());
                itr != voxels.end(); ++itr)
        {
            h5_voxel_array[vidx].lot = (*itr).first.lot();
            h5_voxel_array[vidx].serial = (*itr).first.serial();
            h5_voxel_array[vidx].coordinate = (*itr).second.coordinate();
            ++vidx;
        }

        H5::CompType voxel_comp_type(get_voxel_comp());
        hsize_t dims[] = {num_voxels};
        H5::DataSpace dspace(/* RANK= */1, dims);
        boost::scoped_ptr<H5::DataSet> dset(new H5::DataSet(
            mtgroup->createDataSet("voxels", voxel_comp_type, dspace)));
        dset->write(h5_voxel_array.get(), dset->getDataType());
    }

    template<typename Tspace_>
    static void save_molecular_type_recursively(const Species& location,
            std::multimap<Species, const VoxelPool*>& location_map,
            Tspace_& space, H5::Group* root)
    {
        std::multimap<Species, const VoxelPool*>::iterator itr;
        while ((itr = location_map.find(location)) != location_map.end())
        {
            const VoxelPool* mtb((*itr).second);
            const Species species(mtb->species());
            save_molecular_type(mtb, space.list_voxels_exact(species), root);
            save_molecular_type_recursively(species, location_map, space, root);
            location_map.erase(itr);
        }
    }

    static void sort_by_location(std::multimap<std::string, Species> location_map,
            std::vector<Species>& species, const std::string location="")
    {
        std::multimap<std::string, Species>::iterator itr;
        while ((itr = location_map.find(location)) != location_map.end())
        {
            const Species sp((*itr).second);
            species.push_back(sp);
            sort_by_location(location_map, species, sp.serial());
            location_map.erase(itr);
        }
    }

};

template<typename Tspace_>
void save_lattice_space(const Tspace_& space, H5::Group* root)
{
    typedef LatticeSpaceHDF5Traits traits_type;

    boost::scoped_ptr<H5::Group> spgroup(new H5::Group(root->createGroup("species")));

    const std::vector<Species> species(space.list_species());
    std::multimap<Species, const VoxelPool*> location_map;
    for (std::vector<Species>::const_iterator itr(species.begin());
            itr != species.end(); ++itr)
    {
        const VoxelPool *mtb(space.find_molecular_type(*itr));
        Species location(mtb->location()->species());
        location_map.insert(std::make_pair(location, mtb));
    }
    traits_type::save_molecular_type_recursively(VacantType::getInstance().species(),
            location_map, space, spgroup.get());

    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    const Real3 lengths(space.edge_lengths());

    const uint32_t space_type(static_cast<uint32_t>(Space::LATTICE));
    const double t(space.t());
    const double voxel_radius(space.voxel_radius());
    const uint32_t is_periodic(space.is_periodic()? 1 : 0);
    double edge_lengths[] = {lengths[0], lengths[1], lengths[2]};

#define CREATE_ATTRIBUTE(attribute, type) \
    root->createAttribute(#attribute, type,\
        H5::DataSpace(H5S_SCALAR)).write(type, &attribute)

    CREATE_ATTRIBUTE(space_type, H5::PredType::STD_I32LE);
    CREATE_ATTRIBUTE(t, H5::PredType::IEEE_F64LE);
    CREATE_ATTRIBUTE(voxel_radius, H5::PredType::IEEE_F64LE);
    CREATE_ATTRIBUTE(is_periodic, H5::PredType::STD_I32LE);
    CREATE_ATTRIBUTE(edge_lengths, lengths_type);

#undef CREATE_ATTRIBUTE

}

template<typename Tspace_>
void load_lattice_space(const H5::Group& root, Tspace_* space)
{
    typedef LatticeSpaceHDF5Traits traits_type;

    uint32_t space_type; // not use
    double t;
    double voxel_radius;
    Real3 edge_lengths;
    const hsize_t dims[] = {3};
    uint32_t is_periodic;

#define OPEN_ATTRIBUTE(attribute, type) \
    root.openAttribute(#attribute).read(type, &attribute)

    OPEN_ATTRIBUTE(space_type, H5::PredType::STD_I32LE);
    OPEN_ATTRIBUTE(t, H5::PredType::IEEE_F64LE);
    OPEN_ATTRIBUTE(voxel_radius, H5::PredType::IEEE_F64LE);
    OPEN_ATTRIBUTE(edge_lengths, H5::ArrayType(H5::PredType::NATIVE_DOUBLE, 1, dims));
    OPEN_ATTRIBUTE(is_periodic, H5::PredType::STD_I32LE);

#undef OPEN_ATTRIBUTE

    space->set_t(t);
    space->reset(edge_lengths, voxel_radius, (is_periodic != 0));

    std::map<Species, traits_type::h5_species_struct> struct_map;
    std::map<Species, std::vector<std::pair<ParticleID, Integer> > > voxels_map;
    std::multimap<std::string, Species> location_map;

    std::map<Species, std::pair<traits_type::h5_species_struct,
            std::vector<std::pair<ParticleID, Integer> > > > tmp_map;

    H5::Group spgroup(root.openGroup("species"));
    char name_C[32 + 1];
    for (hsize_t idx(0); idx < spgroup.getNumObjs(); ++idx)
    {
        memset(name_C, 0, 32 + 1);  // clear buffer
        const ssize_t name_len = H5Lget_name_by_idx(spgroup.getLocId(), ".", H5_INDEX_NAME, H5_ITER_INC, idx, name_C, 32, H5P_DEFAULT);
        H5::Group group(spgroup.openGroup(name_C));
        const std::string name_S(name_C);
        Species species(name_S);

        // const H5std_string serial = spgroup.getObjnameByIdx(idx);
        // H5::Group group(spgroup.openGroup(serial.c_str()));
        // Species species(std::string(serial.c_str()));

        traits_type::h5_species_struct property;
        group.openAttribute("property").read(
                traits_type::get_property_comp(), &property);
        struct_map.insert(std::make_pair(species, property));
        location_map.insert(std::make_pair(property.location, species));

        H5::DataSet voxel_dset(group.openDataSet("voxels"));
        const unsigned int num_voxels(
            voxel_dset.getSpace().getSimpleExtentNpoints());
        boost::scoped_array<traits_type::h5_voxel_struct> h5_voxel_array(
                new traits_type::h5_voxel_struct[num_voxels]);
        voxel_dset.read(
                h5_voxel_array.get(), traits_type::get_voxel_comp());
        voxel_dset.close();
        group.close();

        std::vector<std::pair<ParticleID, Integer> > voxels;
        for (unsigned int idx(0); idx < num_voxels; ++idx)
        {
            voxels.push_back(std::make_pair(
                        ParticleID(std::make_pair(h5_voxel_array[idx].lot, h5_voxel_array[idx].serial)),
                        h5_voxel_array[idx].coordinate));
        }
        voxels_map.insert(std::make_pair(species, voxels));
    }
    spgroup.close();

    std::vector<Species> sp_list;
    traits_type::sort_by_location(location_map, sp_list);
    for (std::vector<Species>::iterator itr(sp_list.begin());
            itr != sp_list.end(); ++itr)
    {
        Species species(*itr);
        traits_type::h5_species_struct property((*struct_map.find(species)).second);
        std::vector<std::pair<ParticleID, Integer> > voxels((*voxels_map.find(species)).second);
        if (property.is_structure == 0)
            space->make_molecular_type(species, property.radius, property.D, property.location);
        else
            space->make_structure_type(species, static_cast<Shape::dimension_kind>(property.dimension), property.location);
        space->add_voxels(species, voxels);
    }
}

} // ecell4

#endif /*  __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP */
