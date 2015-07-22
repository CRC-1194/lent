#include "fvCFD.H"
#include "triSurfaceFields.H"

#include <fstream>

#include "auxFunctions.H"

// Compute first parameter of parametrization of an ellipsoid
scalar computeU(vector p)
{
    scalar u = 0.0;

    vector ptilde = p;
    vector ex(1.0, 0.0, 0.0);

    // Project into xy-plane
    ptilde[2] = 0.0;

    scalar utilde = 0.0;
    // Avoid divide-by-zero for vectors which point solely in
    // z-direction. In this case the projection of p results
    // in the zero vector.
    // However, in this case an arbitray value of the interval
    // [0, 2pi] can be chosen for u, thus to leave it zero is fine.
    if (mag(ptilde) > SMALL)
    {
        utilde = Foam::acos(ptilde & ex / mag(ptilde));
    }

    // Since u element of [0,2pi], above expression is not unique;
    // therefore distiction of cases using sign of y-entry
    if (ptilde[1] >= 0.0)
    {
        u = utilde;
    }
    else
    {
        u = 2.0*pi - utilde;
    }

    return u;
}

// Compute second parameter of parametrization of an ellipsoid
scalar computeV(vector p)
{
    scalar v = 0.0;

    vector ez(0.0, 0.0, 1.0);

    v = Foam::acos((p & ez) / mag(p));

    return v;
}

// Analytical function describing the surface of an ellipsoid
point ellipsoidPoint(scalar u, scalar v, vector& axis)
{
    point ePoint(axis[0]*Foam::cos(u)*Foam::sin(v),
                         axis[1]*Foam::sin(u)*Foam::sin(v),
                         axis[2]*Foam::cos(v));

    return ePoint;
}

// Compute deviation of front vertices from the analytical ellipsoid
// normalized with the analytical value
scalar ellipsoidDeviationNormalized(scalar u, scalar v, vector& axis, vector p,
                          vector center)
{
    scalar deviation = 0.0;

    vector pointAnalytic = ellipsoidPoint(u, v, axis);

    // Center of ellipsoid has already been moved to origin
    scalar distanceNumeric = mag(p);
    scalar distanceAnalytic = mag(pointAnalytic);

    deviation = mag(distanceNumeric - distanceAnalytic) / distanceAnalytic;

    return deviation;
}

// Compute normal at point on ellipsoid given the parameters (u,v)
// Vector axis contains length of axis
vector ellipsoidNormal(scalar u, scalar v, vector axis)
{
    vector normal(0.0, 0.0, 0.0);

    vector tangentU(axis[0]*Foam::sin(u)*Foam::sin(v), 
                    -axis[1]*Foam::cos(u)*Foam::sin(v),
                    0.0);
    vector tangentV(axis[0]*Foam::cos(u)*Foam::cos(v),
                    axis[1]*Foam::sin(u)*Foam::cos(v), 
                    -axis[2]*Foam::sin(v));
    
    normal = tangentU ^ tangentV;

    // Tangent approach fails for points coinciding with the
    // z-axis. However, for those points, the normal is simply the
    // basis vector e_z
    if (mag(normal) > SMALL)
    {
        normal = normal / mag(normal);
    }
    else
    {
        normal[0] = 0.0;
        normal[1] = 0.0;

        if (v < pi/2.0)
        {
            normal[2] = 1.0;
        }
        else
        {
            normal[2] = -1.0;
        }

    }
   
    return normal;
}

// Compute twice the mean curvature according to wolfram alpha
scalar ellipsoidCurvature(scalar u, scalar v, vector axis)
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
void checkEllipsoidCurvature(const triSurfacePointVectorField& cn,
                             const triSurface& front, vector center,
                             vector axis, std::fstream& errorFile)
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

    // Metrics to find location of maximum error
    scalar uMax = 0.0;
    scalar vMax = 0.0;
    point pMax(0.0, 0.0, 0.0);
    vector normal(0.0, 0.0, 0.0);
    scalar curvatureNum = 0.0;
    scalar curvatureAna = 0.0;

    forAll(cn, V)
    {
        p = localPoints[V];

        // Move center of ellipsoid to origin
        p = p - center;

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

            // Save max error metrics
            uMax = u;
            vMax = v;
            pMax = localPoints[V];
            curvatureNum = numCurvature;
            curvatureAna = anaCurvature;
            normal = cn[V]/mag(cn[V]);
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

    Info << "Location of maximum error:\n"
         << "u = " << uMax << "; v = " << vMax
         << "\n point = " << pMax
         << "\n Curvature numeric = " << curvatureNum
         << "; Curvature analytic = " << curvatureAna
         << "\n normal vector = " << normal
         << endl;

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
void checkEllipsoidNormal(const triSurfacePointVectorField& cn,
                          const triSurface& front,
                          vector center,
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

        // Move center of ellipsoid to origin
        p = p - center;

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
            // Factor in that curvature normal should point inwards
            // whereas the analytical normal points outwards
            devAngle = (Foam::acos(angle)) * 180.0/pi;
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

// Compute front vertices' deviation from a ellipsoid ssurface
void ellipsoidDeviation(const triSurface& front, vector center, vector axis,
                        std::fstream& errorFile)
{
    scalar deviation = 0.0;
    scalar linearDeviation = 0.0;
    scalar maxDeviation = 0.0;
    scalar counter = 0.0;

    const labelList& vertices = front.meshPoints();
    const pointField& localPoints = front.localPoints();

    scalar u = 0.0;
    scalar v = 0.0;
    point p(0.0, 0.0, 0.0);

    forAll(vertices, V)
    {
        p = localPoints[V];

        // Move center of ellipsoid to origin
        p = p - center;

        u = computeU(p);
        v = computeV(p);

        deviation = ellipsoidDeviationNormalized(u, v, axis, p, center);

        if (deviation > maxDeviation)
        {  
            maxDeviation = deviation;
        }

        linearDeviation += deviation;

        counter++;
    }

    linearDeviation = linearDeviation / counter;

    // Print results
    Info << "\n=== Front deviation from ellipsoid ===" << endl;
    Info << "Linear deviation = " << linearDeviation << endl;
    Info << "Maximum deviation = " << maxDeviation << endl;

    // Write to file
    errorFile << "# Ellipsoid deviation header" << std::endl;
    errorFile << "# n_points\tn_trias\tlinDev\tmaxDev" << std::endl;
    errorFile << front.meshPoints().size() << "\t"
              << front.localFaces().size() << "\t"
              << linearDeviation << "\t"
              << maxDeviation << std::endl;
}

// Experiment: correct the vertices of the gmsh front
void correctFront(triSurface& front, vector center, vector axis)
{
    pointField& frontPoints = const_cast<pointField&> (front.points());

    point p(0.0 ,0.0, 0.0);

    forAll(frontPoints, V)
    {
        p = frontPoints[V];

        p = p - center;

        scalar u = computeU(p);
        scalar v = computeV(p);

        // Compute correct location from parameters and translate
        // the result by center
        frontPoints[V] = ellipsoidPoint(u, v, axis) + center;
    }
}
