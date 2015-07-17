#include "fvCFD.H"
#include "triSurfaceFields.H"

#include <fstream>

#include "auxFunctions.H"

// Compute first parameter of parametrization of an ellipsoid
scalar computeU(vector& p)
{
    scalar u = 0;

    vector ptilde = p;
    vector ex(1,0,0);

    // Project into xy-plane
    ptilde[2] = 0;

    scalar utilde = Foam::acos(ptilde & ex / mag(ptilde));

    // Since u element of [0,2pi], above expression is not unique;
    // therefore distiction of cases using sign of y-entry
    if (ptilde[1] >= 0)
    {
        u = utilde;
    }
    else
    {
        u = 2*pi - utilde;
    }

    return u;
}

// Compute second parameter of parametrization of an ellipsoid
scalar computeV(vector& p)
{
    scalar v = 0;

    vector ptilde = p;

    // Project into xy-plane
    ptilde[2] = 0;

    v = Foam::acos(p & ptilde / (mag(p)*mag(ptilde)));

    return v;
}

// Compute deviation of front vertices from the analytical ellipsoid
// normalized with the analytical value
scalar ellipsoidDeviation(scalar u, scalar v, vector& axis, vector& p,
                          vector& center)
{
    scalar deviation = 0;

    vector pointAnalytic(axis[0]*Foam::cos(u)*Foam::sin(v),
                         axis[1]*Foam::sin(u)*Foam::sin(v),
                         axis[2]*Foam::cos(v));

    scalar distanceNumeric = mag(p - center);
    scalar distanceAnalytic = mag(pointAnalytic - center);

    deviation = mag(distanceNumeric - distanceAnalytic) / distanceAnalytic;

    return deviation;
}

// Compute normal at point on ellipsoid given the parameters (u,v)
// Vector axis contains length of axis
vector ellipsoidNormal(scalar u, scalar v, vector& axis)
{
    vector normal(0.0, 0.0, 0.0);

    vector tangentU(-axis[0]*Foam::sin(u)*Foam::sin(v), 
                    -axis[1]*Foam::cos(u)*Foam::sin(v),
                    0.0);
    vector tangentV(axis[0]*Foam::cos(u)*Foam::cos(v),
                    axis[1]*Foam::sin(u)*Foam::cos(v), 
                    -axis[2]*Foam::sin(v));
    
    normal = tangentU ^ tangentV;
    normal = normal / mag(normal);
   
    return normal;
}

// Compute twice the mean curvature according to wolfram alpha
scalar ellipsoidCurvature(scalar u, scalar v, vector& axis)
{
    scalar curvature = 0;

    scalar a = axis[0];
    scalar b = axis[1];
    scalar c = axis[2];

    scalar denominator = a*a*b*b*Foam::cos(v)*Foam::cos(v) +
                         c*c*(b*b*Foam::cos(u)*Foam::cos(u) + 
                         a*a*Foam::sin(u)*Foam::sin(u))
                         *Foam::sin(v)*Foam::sin(v);
    denominator = denominator*denominator*denominator;
    denominator = 8.0*Foam::sqrt(denominator);


    curvature = a*b*c*( 3*(a*a + b*b) + 2*c*c + (a*a+b*b-2*c*c)*Foam::cos(2*v)
                - 2*(a*a-b*b)*Foam::cos(2*u)*Foam::sin(v)*Foam::sin(v) )
                / denominator;

    return 2*curvature;
}

// Compare exact curvature with numerical curvature at each vertex
void checkEllipsoidCurvature(const triSurfacePointScalarField& cn,
                             const triSurface& front, vector& center,
                             vector& axis, std::fstream& errorFile)
{
    const pointField& localPoints = front.localPoints();

    scalar minCurvature = 1.0e15;
    scalar maxCurvature = 0.0;

    scalar linearDeviation = 0.0;
    scalar maxDeviation = 0.0;

    scalar counter = 0.0;

    scalar numCurvature = 0.0;
    scalar anaCurvature = 0.0;
    scalar deviation = 0.0;

    point p(0,0,0);

    forAll(cn, V)
    {
        p = localPoints[V];

        scalar u = computeU(p);
        scalar v = computeV(p);

        numCurvature = mag(cn[V]);
        anaCurvature = ellipsoidCurvature(u, v, axis);

        // Normalize directly since curvature varies
        deviation = mag(numCurvature - anaCurvature) / anaCurvature;

        // Set extremal values
        if (numCurvature < minCurvature)
        {
            minCurvature = numCurvature;
        }

        if (numCurvature > maxCurvature)
        {
            maxCurvature = numCurvature;
        }

        if (deviation > maxDeviation)
        {
            maxDeviation = deviation;
        }

        // Increment deviations
        linearDeviation += deviation;

        counter++;
    }

    linearDeviation = linearDeviation / counter;

    // Print results
    Info << "\n=== Curvature results ===" << endl;
    Info << "L_inf = " << maxDeviation << endl;
    Info << "Linear deviation = " << linearDeviation << endl;
    Info << "Minimum curvature = " << minCurvature << endl;
    Info << "Maximum curvature = " << maxCurvature << endl;

    // Write to file
    errorFile << "# Curvature header" << std::endl;
    errorFile << "# n_points\tn_trias\tL_inf\t"
              << "linDev\tcurv(min)\tcurv(max)" << std::endl;
    errorFile << front.meshPoints().size() << "\t"
              << front.localFaces().size() << "\t"
              << maxDeviation << "\t"
              << linearDeviation << "\t"
              << minCurvature << "\t"
              << maxCurvature << "\t" << std::endl;
}

// Compare exact front normal with numerical normal
void checkEllipsoidNormal(const triSurfaceVectorField& cn, const triSurface& front,
                          vector axis, std::fstream& errorFile)
{
    const pointField& localPoints = front.localPoints();

    scalar maxDeviation = 0.0;
    scalar devAngle = 0.0;
    scalar linearDeviation = 0.0;

    scalar counter = 0.0;

    vector numNormal(0.0, 0.0, 0.0);
    vector anaNormal(0.0, 0.0, 0.0);
    point p(0.0, 0.0, 0.0);
    scalar angle = 0.0;
    scalar deviation = 0.0;
    
    forAll(localPoints, P)
    {
        p = localPoints[P];

        scalar u = computeU(p);
        scalar v = computeV(p);

        anaNormal = ellipsoidNormal(u, v, axis);
        numNormal = cn[P] / mag(cn[P]);

        // Do not compute actual angle; instead use dot product
        // only as measure for angular deviation
        angle = anaNormal & numNormal;

        deviation = mag(mag(angle) - 1.0);

        if (deviation > maxDeviation)
        {
            maxDeviation = deviation;

            // Compute deviation angle in degree
            devAngle = Foam::acos(angle) * 180.0/pi;
        }

        linearDeviation += deviation;

        counter++;
    }

    // Average
    linearDeviation = linearDeviation / counter;
    //
    // Print results
    Info << "\n=== Normal vector results ===" << endl;
    Info << "Linear deviation = " << linearDeviation << endl;
    Info << "Maximum deviation angle = " << devAngle << endl;

    // Write to file
    errorFile << "# Normal vector header" << std::endl;
    errorFile << "# n points\tn trias\tmax. deviation (degree)\t"
              << "linDev" << std::endl;
    errorFile << front.meshPoints().size() << "\t"
              << front.localFaces().size() << "\t"
              << devAngle << "\t\t"
              << linearDeviation << std::endl;
}
