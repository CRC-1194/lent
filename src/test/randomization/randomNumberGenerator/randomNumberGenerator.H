// TODO: add header
//
// description:
// Random number generator as describes in Stroustrups'
// "The C++ Programming language", section 5.6.3.
// This class behaves as a function object and returns a
// random value in the interval [-1.0,1.0]
// Uses fixed uniform distribution.
// TODO: make the distribution type a template parameter

#ifndef randomNumberGenerator_H
#define randomNumberGenerator_H

#include "scalar.H"

#include <random>

namespace Foam {
namespace FrontTracking {

class randomNumberGenerator
{

private:
    std::random_device seed_;
    std::default_random_engine randomEngine_;
    std::uniform_real_distribution<> distribution_;

public:

    // Constructor
    randomNumberGenerator();
    randomNumberGenerator(scalar low, scalar high);

    // Member operators
    scalar operator()()
    {
        return distribution_(randomEngine_);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
