#include "ODEWorld.hpp"


namespace ecell4
{

namespace ode
{

void ODEWorld::save(const std::string& filename) const
{
    boost::scoped_ptr<H5::H5File>
        fout(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
    boost::scoped_ptr<H5::Group>
        group(new H5::Group(fout->createGroup("CompartmentSpace")));
    save_compartment_space<ODEWorldHDF5Traits<ODEWorld> >(*this, group.get());

    const uint32_t space_type = static_cast<uint32_t>(Space::ELSE);
    group->openAttribute("type").write(H5::PredType::STD_I32LE, &space_type);
}

void ODEWorld::load(const std::string& filename)
{
    boost::scoped_ptr<H5::H5File>
        fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));
    const H5::Group group(fin->openGroup("CompartmentSpace"));
    load_compartment_space<ODEWorldHDF5Traits<ODEWorld> >(group, this);
}

} // ode

} // ecell4
