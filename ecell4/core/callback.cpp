#include <ecell4/core/callback.hpp>

#include <vector>

namespace ecell4
{
CallbackWrapper::CallbackWrapper(stepladder_type stepladder, pyfunc_type pyfunc):
    pyfunc_(pyfunc), stepladder_(stepladder)
{
    ;
}

CallbackWrapper::~CallbackWrapper()
{
    ;
}

void CallbackWrapper::call()
{
    std::vector<double> a;
    a.push_back( double(0.5) );
    a.push_back( double(0.6) );
    a.push_back( double(0.7) );
    if (this->is_available()) {
        this->stepladder_(this->pyfunc_, a);
    }
    ;
}

}
