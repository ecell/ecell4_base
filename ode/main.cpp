#include "MyFunction0.hpp"
#include "ODESolver.hpp"
#include "Defs.hpp"

#include <iostream>
#include <string>

void test_run(ODESolver& s, Time an_end_time)
{

    Time the_time(s.get_current_time());

    while (the_time < an_end_time)
    {
        s.integrate(the_time);

        Integer status_code = s.step();
        if (status_code != 0)break;
      
        std::cout << s.get_current_time() <<  " " 
                  << s.get_value(0) << std::endl;

        the_time = s.reschedule();
    }
    
}

int main()
{
    const Integer a_variable_array_size=1;

    try {
        if (a_variable_array_size <= 0)
        {
            throw 1;
        }
    }
    catch (int e)
    {
        std::cout << 
            "The number of variables (or equations) 'N' must be positive." 
        << std::endl;
        exit(1);
    }

    ODESolver s;

    Real va[a_variable_array_size] = {1.0};

    MyFunction0 f0;
    s.register_function(&f0);
  
    try {
        if (s.get_number_equations() != a_variable_array_size)
        {
            throw 1;
        }
    }
    catch (int e)
    {
        std::cout << 
            "The number of equations must be equal to \
            the number of variables 'N'." 
        << std::endl;
        exit(1);
    }

    const Time a_duration = 10.;
    try {
        if (a_duration <= 0.)
        {
            throw 1;
        }
    }
    catch (int e)
    {
        std::cout << 
            "The time to end this simulation 'aDuration' must be positive." 
        << std::endl;
        exit(1);
    }

    try {
        s.initialize(va, a_variable_array_size);
    
        StatusEvent status_event0 = {0, 0.5, 3};
        s.register_status_event(status_event0);

        test_run(s, a_duration);
    }
    catch (std::string str)
    {
        std::cout << str << std::endl;
    }
    
    return 0;
}

