#include "auxFunctions.H"

// Compute first parameter of the parametrization of an ellipsoid
scalar computeU(vector v)
{
    scalar u = 0;

    vector vtilde = v;
    vector ex(1,0,0);

    // Project into xy-plane
    vtilde[2] = 0;

    scalar utilde = Foam::acos(vtilde & ex / mag(vtilde));

    // Since u element of [0,2pi], above expression is not unique;
    // therefore distiction of cases using sign of y-entry
    if (vtilde[1] >= 0)
    {
        u = utilde;
    }
    else
    {
        u = 2*pi - utilde;
    }

    return u;
}
