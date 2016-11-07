#include <boost/scoped_ptr.hpp>
#include <gsl/gsl_rng.h>
#include <sstream>

#include "RandomNumberGenerator.hpp"

#ifdef WITH_HDF5
#include "extras.hpp"
#endif

namespace ecell4
{

#ifdef WITH_HDF5
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

void GSLRandomNumberGenerator::save(const std::string& filename) const
{
    boost::scoped_ptr<H5::H5File>
        fout(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
    this->save(fout.get());
    extras::save_version_information(fout.get(), std::string("ecell4-gsl_number_generator-") + std::string(ECELL4_VERSION));
}

void GSLRandomNumberGenerator::load(const std::string& filename)
{
    boost::scoped_ptr<H5::H5File>
        fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));
    this->load(*fin);
}
#endif

Real GSLRandomNumberGenerator::random()
{
    return gsl_rng_uniform(rng_.get());
}

Real GSLRandomNumberGenerator::uniform(Real min, Real max)
{
    return gsl_rng_uniform(rng_.get()) * (max - min) + min;
}

Integer GSLRandomNumberGenerator::uniform_int(Integer min, Integer max)
{
    if (max < min)
    {
        throw std::invalid_argument(
            "the max value must be larger than the min value.");
    }

    const unsigned long int n(max - min + 1);
    const unsigned long int range(rng_->type->max - rng_->type->min);

    if (n <= range)
    {
        return gsl_rng_uniform_int(rng_.get(), n) + min;
    }
    else
    {
        const Integer m((max - min) / range);
        Integer k;
        do
        {
            k = min + gsl_rng_uniform_int(rng_.get(), range)
                + range * gsl_rng_uniform_int(rng_.get(), m + 1);
        } while (k > max);
        return k;
    }
}

Real GSLRandomNumberGenerator::gaussian(Real sigma, Real mean)
{
    return gsl_ran_gaussian(rng_.get(), sigma) + mean;
}

Integer GSLRandomNumberGenerator::binomial(Real p, Integer n)
{
    return gsl_ran_binomial(rng_.get(), p, n);
}

Real3 GSLRandomNumberGenerator::direction3d(Real length)
{
    double x, y, z;
    gsl_ran_dir_3d(rng_.get(), &x, &y, &z);
    return Real3(x * length, y * length, z * length);
}

void GSLRandomNumberGenerator::seed(Integer val)
{
    gsl_rng_set(rng_.get(), val);
}

void GSLRandomNumberGenerator::seed()
{
    gsl_rng_set(rng_.get(), unsigned(std::time(0)));
}

} // ecell4
