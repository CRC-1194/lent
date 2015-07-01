#include "fvCFD.H"
#include "triSurfaceFields.H"
#include "interfaceProperties.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"

#include "lentMethod.H"

#include <fstream>

#include "auxFunctions.H"
/*****************************************************************************\
 *
 *      Error metrics and checks
 *
\*****************************************************************************/

// Compare exact curvature with numerical curvature
void checkCurvature(const triSurfacePointVectorField& cn, const triSurface& front,
                    scalar radius, std::fstream& errorFile)
{
    scalar curvatureExact = 2.0 / radius;

    scalar minCurvature = curvatureExact;
    scalar maxCurvature = 0.0;

    scalar averageCurvature = 0.0;
    scalar linearDeviation = 0.0;
    scalar quadDeviation = 0.0;
    scalar maxDeviation = 0.0; // L_inf norm as in Francois' paper

    scalar counter = 0;

    forAll(cn, V)
    {
        scalar numCurvature = mag(cn[V]);
        scalar deviation = mag(numCurvature - curvatureExact);

        // Set extremal values
        if (numCurvature < minCurvature)
        {
            minCurvature = numCurvature;
        }

        if (numCurvature > maxCurvature)
        {
            maxCurvature = numCurvature;
        }

        // Contribution to average value
        averageCurvature += numCurvature;

        // Increment deviations
        linearDeviation += deviation;
        quadDeviation += deviation*deviation;

        counter++;
    }

    averageCurvature = averageCurvature / counter;

    if (mag(curvatureExact - minCurvature) > mag(maxCurvature - curvatureExact))
    {
        maxDeviation = mag(curvatureExact - minCurvature);
    }
    else
    {
        maxDeviation = mag(maxCurvature - curvatureExact);
    }

    // Normalize with exact curvature for comparability
    linearDeviation = linearDeviation / (counter*curvatureExact);
    quadDeviation = Foam::sqrt(quadDeviation / (counter*curvatureExact));
    maxDeviation = maxDeviation / curvatureExact;

    // Print results
    Info << "\n=== Curvature results ===" << endl;
    Info << "Exact curvature = " << curvatureExact << endl;
    Info << "L_inf = " << maxDeviation << endl;
    Info << "Average curvature = " << averageCurvature << endl;
    Info << "Linear deviation = " << linearDeviation << endl;
    Info << "Quadratic deviation = " << quadDeviation << endl;
    Info << "Minimum curvature = " << minCurvature << endl;
    Info << "Maximum curvature = " << maxCurvature << endl;

    // Write to file
    errorFile << "# Curvature header" << std::endl;
    errorFile << "# n_points\tn_trias\tcurv(exact)\tL_inf\tcurv(avg.)\t"
              << "linDev\tquadDev\tcurv(min)\tcurv(max)" << std::endl;
    errorFile << front.meshPoints().size() << "\t"
              << front.localFaces().size() << "\t"
              << curvatureExact << "\t"
              << maxDeviation << "\t"
              << averageCurvature << "\t"
              << linearDeviation << "\t"
              << quadDeviation << "\t"
              << minCurvature << "\t"
              << maxCurvature << "\t" << std::endl;
}

// Compare exact front normal with numerical normal
void checkNormal(const triSurfacePointVectorField& cn, const triSurface& front,
                 vector center, std::fstream& errorFile)
{
    scalar maxDeviation = 0.0;
    scalar devAngle = 0.0;
    scalar linearDeviation = 0.0;
    scalar quadDeviation = 0.0;

    scalar counter = 0;

    const pointField& localPoints = front.localPoints();

    forAll(localPoints, P)
    {
        vector normalExact = localPoints[P] - center;
        vector normalNumerical = cn[P];

        // Normalize vectors
        normalExact = normalExact / mag(normalExact);
        normalNumerical = normalNumerical / mag(normalNumerical);

        // Do not compute actual angle; instead use dot product
        // only as measure for angular deviation
        scalar angle = normalExact & normalNumerical;

        scalar deviation = mag(mag(angle) - 1.0);

        if (deviation > maxDeviation)
        {
            maxDeviation = deviation;

            // Compute deviation angle in degree
            devAngle = Foam::acos(angle) * 57.296;
        }

        linearDeviation += deviation;
        quadDeviation += deviation*deviation;

        counter++;
    }

    // Average
    linearDeviation = linearDeviation / counter;
    quadDeviation = Foam::sqrt(quadDeviation / counter);

    // Print results
    Info << "\n=== Normal vector results ===" << endl;
    Info << "Linear deviation = " << linearDeviation << endl;
    Info << "Quadratic deviation = " << quadDeviation << endl;
    Info << "Maximum deviation angle = " << devAngle << endl;

    // Write to file
    errorFile << "# Normal vector header" << std::endl;
    errorFile << "# n points\tn trias\tmax. deviation (degree)\t"
              << "linDev\tquadDev" << std::endl;
    errorFile << front.meshPoints().size() << "\t"
              << front.localFaces().size() << "\t"
              << devAngle << "\t\t"
              << linearDeviation << "\t"
              << quadDeviation << std::endl;
}

// Compute front vertices' deviation from a sphere surface
void sphereDeviation(const triSurface& front, scalar radius, vector center,
                     std::fstream& errorFile)
{
    scalar averageDeviation = 0.0;
    scalar maxDeviation = 0.0;
    scalar counter = 0;

    const labelList& vertices = front.meshPoints();
    const pointField& localPoints = front.localPoints();

    forAll(vertices, V)
    {
        scalar deviation = mag(radius - mag(center - localPoints[V]));

        if (deviation > maxDeviation)
        {
            maxDeviation = deviation;
        }

        averageDeviation += deviation;

        counter++;
    }

    // Normalize results
    averageDeviation = averageDeviation / (counter*radius);
    maxDeviation = maxDeviation / radius;

    // Print results
    Info << "\n=== Front deviation from sphere ===" << endl;
    Info << "Linear deviation = " << averageDeviation << endl;
    Info << "Maximum deviation = " << maxDeviation << endl;

    // Write to file
    errorFile << "# Sphere deviation header" << std::endl;
    errorFile << "# n_points\tn_trias\tlinDev\tmaxDev" << std::endl;
    errorFile << front.meshPoints().size() << "\t"
              << front.localFaces().size() << "\t"
              << averageDeviation << "\t"
              << maxDeviation << std::endl;
}


// Determine mesh quality
void meshQuality(const triSurface& front, std::fstream& file)
{
    // Get references to surface mesh
    const pointField& localPoints = front.localPoints();
    const List<labelledTri>& localFaces = front.localFaces();

    // Extrema of mesh properties
    scalar minAngle = 2.0*pi;  
    scalar maxAngle = 0.0;   

    scalar maxRatio = 1.0; 

    scalar minArea = 1000.0; 
    scalar maxArea = 0.0; 

    forAll(localFaces, T)
    {
        labelledTri tri = localFaces[T];

        label lx1 = tri[0];
        label lx2 = tri[1];
        label lx3 = tri[2];

        // points
        point x1 = localPoints[lx1];
        point x2 = localPoints[lx2];
        point x3 = localPoints[lx3];

        // edges
        vector a = x1 - x2;
        vector b = x1 - x3;
        vector c = x2 - x3;

        // angles
        scalar alpha = getAngle(b, c);
        scalar gamma = getAngle(a, b);
        scalar beta = pi - (alpha + gamma);

        // area
        scalar area = areaHeron(a, b, c);

        // angle check
        minAngle = alpha < minAngle ? alpha : minAngle;
        minAngle = beta < minAngle ? beta : minAngle;
        minAngle = gamma < minAngle ? gamma : minAngle;

        maxAngle = alpha > maxAngle ? alpha : maxAngle;
        maxAngle = beta > maxAngle ? beta : maxAngle;
        maxAngle = gamma > maxAngle ? gamma : maxAngle;

        // edge ratio check
        scalar minlength = mag(a);
        scalar maxlength = mag(b);

        minlength = mag(c) < minlength ? mag(c) : minlength;
        minlength = mag(b) < minlength ? mag(b) : minlength;
        maxlength = mag(c) > maxlength ? mag(c) : maxlength;
        maxlength = mag(a) > maxlength ? mag(a) : maxlength;

        // Edge ratio check
        if (maxlength/minlength > maxRatio)
        {
            maxRatio = maxlength/minlength;
        }

        // area check
        minArea = area < minArea ? area : minArea;
        maxArea = area > maxArea ? area : maxArea;
    }

    // Write results
    file << "# Mesh metrics:\n"
         << "# Minimum angle = " << 57.3*minAngle << "\n"
         << "# Maximum angle = " << 57.3*maxAngle << "\n"
         << "# Max edge ratio = " << maxRatio << "\n"
         << "# Min area = " << minArea << "\n"
         << "# Max Area = " << maxArea << std::endl;
}
