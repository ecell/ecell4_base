#ifndef __ECELL4_CALLBACK_HPP
#define __ECELL4_CALLBACK_HPP

#include <vector>
#include <ecell4/core/types.hpp>
#include <ecell4/core/exceptions.hpp>

namespace ecell4
{


//template <typename T>
class CallbackWrapper
{
public:
    typedef void *pyfunc_type;
    typedef int (*stepladder_type)(pyfunc_type, std::vector<double>);

    CallbackWrapper(stepladder_type stepladder, pyfunc_type pyfunc);
    ~CallbackWrapper();
    bool is_available()
    {
        return (this->pyfunc_ != NULL && this->stepladder_ != NULL);
    }
    void call();
protected:
    pyfunc_type pyfunc_;
    stepladder_type stepladder_;
};


}
#endif /* __ECELL4_HPP */

