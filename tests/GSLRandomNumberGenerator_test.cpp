#include <gtest/gtest.h>

#include "../RandomNumberGenerator.hpp"

using namespace ecell4;


TEST(GSLRandomNumberGeneratorTest, Seed)
{
    GSLRandomNumberGenerator rng;
    rng.seed(0);
}
