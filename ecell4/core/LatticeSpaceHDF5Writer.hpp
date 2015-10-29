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
#include "MolecularTypeBase.hpp"


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
#define INSERT_MEMBER(member, type) \
        property_comp_type.insertMember(#member,\
                HOFFSET(h5_species_struct, member), type)
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
#define INSERT_MEMBER(member, type) \
        voxel_comp_type.insertMember(std::string(#member),\
            HOFFSET(h5_voxel_struct, member), type)
        INSERT_MEMBER(lot, H5::PredType::NATIVE_INT);
        INSERT_MEMBER(serial, H5::PredType::NATIVE_INT);
        INSERT_MEMBER(coordinate, H5::PredType::STD_I64LE);
#undef INSERT_MEMBER
        return voxel_comp_type;
    }

    static void save_molecular_type(const MolecularTypeBase* mtb,
            std::vector<std::pair<ParticleID, Voxel> > voxels, H5::Group* group)
    {
        const Species species(mtb->species());
        boost::scoped_ptr<H5::Group> mtgroup(
                new H5::Group(group->createGroup(species.serial())));

        h5_species_struct property;
        property.radius = mtb->radius();
        property.D = mtb->D();
        const MolecularTypeBase* loc(mtb->location());
        if (loc->is_vacant())
            std::strcpy(property.location, "");
        else
            std::strcpy(property.location, loc->species().serial().c_str());
        property.is_structure = mtb->is_structure() ? 1 : 0;
        property.dimension = loc->is_vacant() ? mtb->get_dimension() : Shape::UNDEF;

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

};

template<typename Tspace_>
void save_lattice_space(const Tspace_& space, H5::Group* root)
{
    typedef LatticeSpaceHDF5Traits traits_type;

    boost::scoped_ptr<H5::Group> spgroup(new H5::Group(root->createGroup("species")));

    const std::vector<Species> species(space.list_species());
    /*
     * TODO sort species with location
     */
    for (std::vector<Species>::const_iterator itr(species.begin());
            itr != species.end(); ++itr)
    {
        const MolecularTypeBase *mtb(space.find_molecular_type(*itr));
        traits_type::save_molecular_type(mtb, space.list_voxels(*itr), spgroup.get());
    }

#define CREATE_ATTRIBUTE(attribute, type) \
    root->createAttribute(#attribute, type,\
        H5::DataSpace(H5S_SCALAR)).write(type, &attribute)

    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    const Real3 lengths(space.edge_lengths());

    const uint32_t space_type(static_cast<uint32_t>(Space::LATTICE));
    const double t(space.t());
    const double voxel_radius(space.voxel_radius());
    const uint32_t is_periodic(space.is_periodic()? 1 : 0);
    double edge_lengths[] = {lengths[0], lengths[1], lengths[2]};

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

    H5::Group spgroup(root.openGroup("species"));
    for (hsize_t idx(0); idx < spgroup.getNumObjs(); ++idx)
    {
        std::string serial(spgroup.getObjnameByIdx(idx));
        H5::Group group(spgroup.openGroup(serial));
        traits_type::h5_species_struct property;
        group.openAttribute("property").read(
                traits_type::get_property_comp(), &property);
        /*
        Species species(serial, property.radius, property.D,
                std::string(property.location));
                */

        H5::DataSet voxel_dset(group.openDataSet("voxels"));
        const unsigned int num_voxels(
            voxel_dset.getSpace().getSimpleExtentNpoints());
        boost::scoped_array<traits_type::h5_voxel_struct> h5_voxel_array(
                new traits_type::h5_voxel_struct[num_voxels]);
        voxel_dset.read(
                h5_voxel_array.get(), traits_type::get_voxel_comp());
        voxel_dset.close();

        std::vector<std::pair<ParticleID, Integer> > voxels;
        for (unsigned int idx(0); idx < num_voxels; ++idx)
        {
            voxels.push_back(std::make_pair(
                        ParticleID(std::make_pair(h5_voxel_array[idx].lot, h5_voxel_array[idx].serial)),
                        h5_voxel_array[idx].coordinate));
        }

        group.close();
    }

}

} // ecell4

#endif /*  __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP */
