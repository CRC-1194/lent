CorrectPhi
(
    U,
    phi,
    p_rgh,
    dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1),
    geometricZeroField(),
    lentSC
);

#include "continuityErrs.H"
