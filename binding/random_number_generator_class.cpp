#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "RandomNumberGenerator.hpp"
#include "../gsl_rng_base.hpp"
#include "binding_common.hpp"

namespace binding {

struct static_gsl_rng: public gsl_rng_base<static_gsl_rng>
{
    typedef unsigned int result_type;

    static const char name[];
    static const unsigned long int min = 0;
    static const unsigned long int max = static_cast<result_type>(-1);

    void set(unsigned long int seed)
    {
        idx_ = std::min(seed, static_cast<unsigned long int>(PyObject_Size(seq_.ptr())));
    }

    unsigned long int get()
    {
        Py_ssize_t nelems(PyObject_Size(seq_.ptr()));
        if (idx_ >= nelems)
        {
            return min; 
        }
        boost::python::handle<> i(
                boost::python::allow_null(
                    PySequence_GetItem(seq_.ptr(), idx_)));
        if (!i)
        {
            return min;
        }
        ++idx_;
        return PyLong_AsUnsignedLong(i.get());
    }

    double get_double()
    {
        return static_cast<double>(get()) / (static_cast<double>(max) + 1);
    }

    static_gsl_rng(boost::python::object seq)
        : seq_(seq), idx_(0) {}

private:
    boost::python::object seq_;
    Py_ssize_t idx_;
};

const char static_gsl_rng::name[] = "static_gsl_rng";

template<gsl_rng_type const*& Prng_>
static GSLRandomNumberGenerator create_gsl_rng()
{
    return GSLRandomNumberGenerator(gsl_rng_alloc(Prng_));
}

static GSLRandomNumberGenerator create_static_gsl_rng(boost::python::object seq)
{
    return GSLRandomNumberGenerator(
            GSLRandomNumberGenerator::rng_handle(new static_gsl_rng(seq)));
}

void register_random_number_generator_class()
{
    using namespace boost::python;
    register_random_number_generator_class<GSLRandomNumberGenerator>("RandomNumberGenerator");
    def("create_gsl_rng", &create_gsl_rng<gsl_rng_mt19937>);
    def("create_static_gsl_rng", &create_static_gsl_rng);
}

} // namespace binding
