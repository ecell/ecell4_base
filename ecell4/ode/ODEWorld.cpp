#include "ODEWorld.hpp"
#include <sstream>
#include <ecell4/core/extras.hpp>
#include <ecell4/core/exceptions.hpp>

namespace ecell4
{

namespace ode
{

void ODEWorld::bind_to(std::shared_ptr<Model> model)
{
    if (std::shared_ptr<Model> bound_model = lock_model())
    {
        if (bound_model.get() != model.get())
        {
            std::cerr << "Warning: Model already bound to BDWorld"
                << std::endl;
        }
    }

    model_ = model;
}

void ODEWorld::save(const std::string& filename) const
{
#ifdef WITH_HDF5
    std::unique_ptr<H5::H5File>
        fout(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
    std::unique_ptr<H5::Group>
        group(new H5::Group(fout->createGroup("CompartmentSpace")));
    save_compartment_space<ODEWorldHDF5Traits<ODEWorld> >(*this, group.get());

    const uint32_t space_type = static_cast<uint32_t>(Space::ELSE);
    group->openAttribute("type").write(H5::PredType::STD_I32LE, &space_type);

    extras::save_version_information(fout.get(), std::string("ecell4-ode-") + std::string(VERSION_INFO));
#else
    throw NotSupported(
        "This method requires HDF5. The HDF5 support is turned off.");
#endif
}

void ODEWorld::load(const std::string& filename)
{
#ifdef WITH_HDF5
    std::unique_ptr<H5::H5File>
        fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));

    const std::string required = "ecell4-ode-0.0";
    try
    {
        const std::string version = extras::load_version_information(*fin);
        if (!extras::check_version_information(version, required))
        {
            std::stringstream ss;
            ss << "The version of the given file [" << version
                << "] is too old. [" << required << "] or later is required.";
            throw NotSupported(ss.str());
        }
    }
    catch(H5::GroupIException not_found_error)
    {
        throw NotFound("No version information was found.");
    }

    const H5::Group group(fin->openGroup("CompartmentSpace"));
    load_compartment_space<ODEWorldHDF5Traits<ODEWorld> >(group, this);
#else
    throw NotSupported(
        "This method requires HDF5. The HDF5 support is turned off.");
#endif
}

} // ode

} // ecell4
