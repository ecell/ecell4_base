#include <boost/scoped_ptr.hpp>
#include <gsl/gsl_rng.h>
#include <sstream>

#include <hdf5.h>
#include <H5Cpp.h>

#include "RandomNumberGenerator.hpp"


namespace ecell4
{

void GSLRandomNumberGenerator::save(H5::CommonFG* root) const
{
    using namespace H5;

    boost::scoped_ptr<DataType> optype(new DataType(H5T_OPAQUE, 1));
    hsize_t bufsize(gsl_rng_size(rng_.get()));
    DataSpace dataspace(1, &bufsize);
    optype->setTag("GSLRandomNumberGenerator state type");
    boost::scoped_ptr<DataSet> dataset(
        new DataSet(root->createDataSet("rng", *optype, dataspace)));
    dataset->write((unsigned char*)(gsl_rng_state(rng_.get())), *optype);
}

void GSLRandomNumberGenerator::load(const H5::CommonFG& root)
{
    using namespace H5;

    const DataSet dataset(DataSet(root.openDataSet("rng")));
    // size_t bufsize(gsl_rng_size(rng_.get()));
    boost::scoped_ptr<DataType> optype(new DataType(H5T_OPAQUE, 1));
    optype->setTag("GSLRandomNumberGenerator state type");
    unsigned char* state = (unsigned char*)(gsl_rng_state(rng_.get()));
    dataset.read(state, *optype);
}

} // ecell4
