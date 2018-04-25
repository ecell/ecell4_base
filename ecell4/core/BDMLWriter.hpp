#ifndef ECELL4_BDML_WRITER_HPP
#define ECELL4_BDML_WRITER_HPP

#include <cstring>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <sstream>

#include "types.hpp"
#include "get_mapper_mf.hpp"
#include "Species.hpp"
#include "Particle.hpp"
#include "Space.hpp"
#include "exceptions.hpp"

#ifdef WITH_HDF5

#include <hdf5.h>
#include <H5Cpp.h>

namespace ecell4
{

struct BDMLTraits
{
    typedef struct bdml_objectDef_struct {
        uint32_t oid;
        char name[32];
    } bdml_objectDef_struct;

    static H5::CompType get_bdml_objectDef_comp_type()
    {
        H5::CompType comp_type(sizeof(bdml_objectDef_struct));
#define INSERT_MEMBER(member, type) \
        H5Tinsert(comp_type.getId(), #member,\
                HOFFSET(bdml_objectDef_struct, member), type.getId())
        INSERT_MEMBER(oid, H5::PredType::STD_I32LE);
        INSERT_MEMBER(name, H5::StrType(H5::PredType::C_S1, 32));
#undef INSERT_MEMBER
        return comp_type;
    }

    typedef struct bdml_scaleUnit_struct {
        char dimension[32];
        double xscale;
        double yscale;
        double zscale;
        char sunit[32];
        double tscale;
        char tunit[32];
    } bdml_scaleUnit_struct;

    static H5::CompType get_bdml_scaleUnit_comp_type()
    {
        H5::CompType comp_type(sizeof(bdml_scaleUnit_struct));
#define INSERT_MEMBER(member, type) \
        H5Tinsert(comp_type.getId(), #member,\
                HOFFSET(bdml_scaleUnit_struct, member), type.getId())
        INSERT_MEMBER(dimension, H5::StrType(H5::PredType::C_S1, 32));
        INSERT_MEMBER(xscale, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(yscale, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(zscale, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(sunit, H5::StrType(H5::PredType::C_S1, 32));
        INSERT_MEMBER(tscale, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(tunit, H5::StrType(H5::PredType::C_S1, 32));
#undef INSERT_MEMBER
        return comp_type;
    }

    typedef struct bdml_sphere_struct {
        char id[32];
        double t;
        char entity[32];
        double x;
        double y;
        double z;
        double radius;
        char optional[32];
    } bdml_sphere_struct;

    static H5::CompType get_bdml_sphere_comp_type()
    {
        H5::CompType comp_type(sizeof(bdml_sphere_struct));
#define INSERT_MEMBER(member, type) \
        H5Tinsert(comp_type.getId(), #member,\
                HOFFSET(bdml_sphere_struct, member), type.getId())
        INSERT_MEMBER(id, H5::StrType(H5::PredType::C_S1, 32));
        INSERT_MEMBER(t, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(entity, H5::StrType(H5::PredType::C_S1, 32));
        INSERT_MEMBER(x, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(y, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(z, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(radius, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(optional, H5::StrType(H5::PredType::C_S1, 32));
#undef INSERT_MEMBER
        return comp_type;
    }
};

// template <typename Tspace_>
void save_bd5(
    const Space& space, const std::string& filename,
    const int group_index,
    const std::string& object_name,
    const std::string& spatial_unit,
    const std::string& time_unit,
    const bool trunc
    )
{
    //XXX: group_index = 0
    //XXX: object_name = "molecule"
    //XXX: spatial_unit = "meter"
    //XXX: time_unit = "second"

    typedef BDMLTraits traits_type;

    assert(group_index >= 0);
    const std::string group_name
        = static_cast<std::ostringstream&>(std::ostringstream() << std::dec << group_index).str();

    H5E_auto2_t func;
    void* client_data;
    H5::Exception::getAutoPrint(func, &client_data);

    boost::scoped_ptr<H5::H5File> fout;
    bool file_exists = false;
    if (trunc)
    {
        boost::scoped_ptr<H5::H5File> tmp(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
        fout.swap(tmp);
    }
    else
    {
        H5::Exception::dontPrint();
        try
        {
            boost::scoped_ptr<H5::H5File> tmp(new H5::H5File(filename.c_str(), H5F_ACC_RDWR));
            fout.swap(tmp);
            H5::Exception::setAutoPrint(func, client_data);
            file_exists = true;
        }
        catch (H5::FileIException &file_exists_err)
        {
            H5::Exception::setAutoPrint(func, client_data);
            boost::scoped_ptr<H5::H5File> tmp(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
            fout.swap(tmp);
        }
    }

    boost::scoped_ptr<H5::Group> data;

    if (!file_exists)
    {
        boost::scoped_ptr<H5::Group>
            tmp(new H5::Group(fout->createGroup("data")));
        data.swap(tmp);

        {
            boost::scoped_array<traits_type::bdml_objectDef_struct>
                objectDef_table(new traits_type::bdml_objectDef_struct[1]);
            objectDef_table[0].oid = 0;
            std::strcpy(objectDef_table[0].name, object_name.c_str());

            const int RANK = 1;
            hsize_t dim[] = {1};
            H5::DataSpace dataspace(RANK, dim);
            boost::scoped_ptr<H5::DataSet> objectDef(
                new H5::DataSet(data->createDataSet("objectDef", traits_type::get_bdml_objectDef_comp_type(), dataspace)));
            objectDef->write(objectDef_table.get(), objectDef->getDataType());
        }

        {
            boost::scoped_array<traits_type::bdml_scaleUnit_struct>
                scaleUnit_table(new traits_type::bdml_scaleUnit_struct[1]);
            std::strcpy(scaleUnit_table[0].dimension, "3D+T");
            scaleUnit_table[0].xscale = 1.0;
            scaleUnit_table[0].yscale = 1.0;
            scaleUnit_table[0].zscale = 1.0;
            std::strcpy(scaleUnit_table[0].sunit, spatial_unit.c_str());
            scaleUnit_table[0].tscale = 1.0;
            std::strcpy(scaleUnit_table[0].tunit, time_unit.c_str());

            const int RANK = 1;
            hsize_t dim[] = {1};
            H5::DataSpace dataspace(RANK, dim);
            boost::scoped_ptr<H5::DataSet> scaleUnit(
                new H5::DataSet(data->createDataSet("scaleUnit", traits_type::get_bdml_scaleUnit_comp_type(), dataspace)));
            scaleUnit->write(scaleUnit_table.get(), scaleUnit->getDataType());
        }
    }
    else
    {
        boost::scoped_ptr<H5::Group>
            tmp(new H5::Group(fout->openGroup("data")));
        data.swap(tmp);
    }

    H5::Exception::dontPrint();
    try
    {
        data->openGroup(group_name.c_str());
        H5::Exception::setAutoPrint(func, client_data);
        std::stringstream ss;
        ss << "Group [" << group_name << "] already exists. Do nothing";
        throw AlreadyExists(ss.str());
    }
    catch (H5::Exception &err)
    {
        H5::Exception::setAutoPrint(func, client_data);

        boost::scoped_ptr<H5::Group>
            data_zero(new H5::Group(data->createGroup(group_name.c_str())));
        boost::scoped_ptr<H5::Group>
            data_zero_object(new H5::Group(data_zero->createGroup("object")));

        typedef std::vector<std::pair<ParticleID, Particle> >
            particle_container_type;
        const particle_container_type particles(space.list_particles());
        // const particle_container_type& particles(space.list_particles());
        const unsigned int NUM_MOL = space.num_particles();

        boost::scoped_array<traits_type::bdml_sphere_struct>
            data_table(new traits_type::bdml_sphere_struct[NUM_MOL]);
        for (unsigned int i(0); i < NUM_MOL; ++i)
        {
            const Particle& p(particles[i].second);

            std::strcpy(data_table[i].id, "1");
            data_table[i].t = space.t();
            std::strcpy(data_table[i].entity, "sphere");
            data_table[i].x = p.position()[0];
            data_table[i].y = p.position()[1];
            data_table[i].z = p.position()[2];
            data_table[i].radius = p.radius();
            std::strcpy(data_table[i].optional, "none");
        }

        const int RANK = 1;
        hsize_t dim[] = {NUM_MOL};
        H5::DataSpace dataspace(RANK, dim);
        boost::scoped_ptr<H5::DataSet> data_zero_object_zero(
            new H5::DataSet(data_zero_object->createDataSet("0", traits_type::get_bdml_sphere_comp_type(), dataspace)));
        data_zero_object_zero->write(data_table.get(), data_zero_object_zero->getDataType());
    }
}

}; // ecell4

#else // WITH_HDF5

namespace ecell4
{

void save_bdml(
    const Space& space, const std::string& filename,
    const std::string& group_name,
    const std::string& object_name,
    const std::string& spatial_unit,
    const std::string& time_unit,
    const bool trunc
    )
{
    throw NotSupported(
        "This method requires HDF5. The HDF5 support is turned off.");
}

}; // ecell4

#endif // WITH_HDF5

#endif /* ECELL4_BDML_WRITER_HPP */
