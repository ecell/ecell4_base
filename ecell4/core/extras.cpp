#include "extras.hpp"


namespace ecell4
{

namespace extras
{

#ifdef WITH_HDF5
void save_version_information(H5::CommonFG* root, const std::string& version)
{
    using namespace H5;
    boost::scoped_ptr<DataSet> dataset(
        new DataSet(root->createDataSet("version", H5::StrType(H5::PredType::C_S1, 32), H5::DataSpace(H5S_SCALAR))));
    dataset->write(version.c_str(), dataset->getDataType());
}

std::string load_version_information(const H5::CommonFG& root)
{
    using namespace H5;
    char buf[32];
    const DataSet dataset(DataSet(root.openDataSet("version")));
    dataset.read(buf, dataset.getDataType());
    return std::string(buf);
}
#endif

std::string load_version_information(const std::string& filename)
{
#ifdef WITH_HDF5
    using namespace H5;
    boost::scoped_ptr<H5::H5File>
        fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));
    return load_version_information(*fin);
#else
    return "";
#endif
}

} // extras

} // ecell4
