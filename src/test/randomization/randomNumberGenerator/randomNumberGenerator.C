#include "randomNumberGenerator.H"

namespace Foam {
namespace FrontTracking {


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
randomNumberGenerator::randomNumberGenerator()
:
    seed_{},
    randomEngine_{seed_()},
    distribution_{-1.0,1.0}
{}

randomNumberGenerator::randomNumberGenerator(scalar low, scalar high)
:
    seed_{},
    randomEngine_{seed_()},
    distribution_{low,high}
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
