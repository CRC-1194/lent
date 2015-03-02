/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2.2.x                               
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Test application for the interface reconstruction algorithm of the LENT method.

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "triSurfaceFields.H"
#include "interfaceProperties.H"
#include "incompressibleTwoPhaseMixture.H"
#include "lentMethod.H"

#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Aux functions

// Calculate angle between two vectors
scalar getAngle(vector& a, vector& b)
{
    return Foam::acos((a & b) / (mag(a) * mag(b)));
}

// Cotangent function to resemble algorihm notion from paper
scalar cot(scalar angle)
{
    return 1.0/Foam::tan(angle);
}

// Triangle area functions
// Use heron's formula
scalar areaHeron(vector& A, vector& B, vector& C)
{
    scalar a = mag(A);
    scalar b = mag(B);
    scalar c = mag(C);
    scalar s = 0.5 * (a + b + c);
                
    return Foam::sqrt(s*(s-a)*(s-b)*(s-c));
}

// Use cross product property
scalar areaCross(vector& a, vector& b)
{
    return 0.5*mag(a ^ b);
}

// Quality metric in order to get to the core of the problem
// of completely wrong curvatures
bool badQuality(point& A, point& B, point& C)
{
    bool isBad = false;

    scalar pi = 3.141592653589793238; 

    // Edges
    vector a = C - B;
    vector b = C - A;
    vector c = B - A;

    // Angles
    scalar alpha = getAngle(b, c);
    scalar gamma = getAngle(a, b);
    scalar beta = pi - (alpha + gamma);

    // Area
    scalar area = areaHeron(a, b, c);

    // Actual quality check, hardcoded for now
    scalar minAngle = 0.17453;  // 10 degree
    scalar maxAngle = 2.4435;   // 140 degree

    scalar maxRatio = 5.0; // maximum ratio of edge lengths

    scalar minArea = 1e-5; // Should be smaller than the value computed by checkSTL
    scalar maxArea = 1e-1; // Should be bigger than value computed by checkSTL

    // Angle check
    if (alpha < minAngle || beta < minAngle || gamma < minAngle)
    {
        isBad = true;
    }

    if (alpha > maxAngle || beta > maxAngle || gamma > maxAngle)
    {
        isBad = true;
    }

    // Edge ratio check
    scalar minLength = mag(a);
    scalar maxLength = mag(b);

    minLength = mag(c) < minLength ? mag(c) : minLength;
    minLength = mag(b) < minLength ? mag(b) : minLength;
    maxLength = mag(c) > maxLength ? mag(c) : maxLength;
    maxLength = mag(a) > maxLength ? mag(a) : maxLength;

    isBad = maxRatio > maxLength/minLength ? false : true;

    // Area check
    if (area < minArea || area > maxArea) {isBad = true;}

    return isBad;
}

void curvatureNormals(triSurfacePointVectorField& cn, const triSurface& front)
{
    scalar pi = 3.141592653589793238; 

    cn = dimensionedVector("zero",
                            dimless/dimLength,
                            vector(0,0,0)
                          );

    // TODO: modify so that mutliple patches, e.g. in the case of
    // bubble breakup, are supported

    // All of the following functions are taken from
    // "PrimitivePatch.H"
    // Get label list of all mesh vertices
    const labelList& vertices = front.meshPoints();

    // Get assignment point --> faces
    const labelListList& adjacentFaces = front.pointFaces();

    // List of local point references for current patch
    const pointField& localPoints = front.localPoints();
    
    // Needed: List of local face references
    const List<labelledTri>& localFaces = front.localFaces();

    // V (the actual vertex) and Vl (its label) are used synonymously
    // in the following comments
    forAll(vertices, Vl)
    {
        scalar Amix = 0;

        // Get all triangles adjacent to V
        const labelList& oneRingNeighborhood = adjacentFaces[Vl];

        // Sum up mixed area and contributions to mean curvature normal
        forAll(oneRingNeighborhood, T)
        {
            // Compute area contributions
            // First, resolve label T to get the actual face
            label triLabel = oneRingNeighborhood[T];
            labelledTri currentTri = localFaces[triLabel];

            // Read point labels from triangle and determine which
            // one matches Vl to resolve point coordinates correctly
            label tri0 = currentTri[0];
            label tri1 = currentTri[1];
            label tri2 = currentTri[2];

            point V;
            point Q;
            point R;

            if (tri0 == Vl)
            {
                V = localPoints[Vl];
                Q = localPoints[tri1];
                R = localPoints[tri2];
            }
            else if (tri1 == Vl)
            {
                V = localPoints[Vl];
                Q = localPoints[tri0];
                R = localPoints[tri2];
            }
            else
            {
                V = localPoints[Vl];
                Q = localPoints[tri0];
                R = localPoints[tri1];
            }

            // Edge vectors
            vector VQ = Q - V;
            vector VR = R - V;
            vector QR = R - Q;

            // Get angles
            scalar Va = getAngle(VQ, VR);
            scalar Ra = getAngle(VR, QR);
            scalar Qa = pi - (Va + Ra);

            /*
            if (badQuality(V, Q, R))
            {
                Info << "Bad triangle detected:\n"
                     << "Label: l = " << oneRingNeighborhood[T] << "\n"
                     << "Points:\nV" << V << "\nQ" << Q << "\nR" << R
                     << " \n"
                     << "Angles: v = " << Va*57.3 << ", q = " << Qa*57.3
                     << ", r = " << Ra*57.3 << "\n" << endl;
            }
            */

            // Check if non-obtuse in order to use the correct area metric
            if (Va < pi/2 && Qa < pi/2 && Ra < pi/2)
            {
                // Use Voronoi-area
                // Cotangent function 'cot' has to be defined locally since
                // it is not offered by OpenFOAM
                Amix += 0.125*(magSqr(VR)*cot(Qa) + magSqr(VQ)*cot(Ra));    
            }
            else if (Va >= pi/2)
            {
                // Obtuse angle at V, use half area
                // Use Heron's formula for now until problem
                // with vector approach is solved
                Amix += 0.5 * areaHeron(VQ, VR, QR);
            }
            else
            {
                // Obtuse angle not at V, use quarter area
                // Use Heron's formula for now until problem
                // with vector approach is solved
                Amix += 0.25 * areaHeron(VQ, VR, QR);
            }

            // Compute mean curvature normal contributions
            cn[Vl] += cot(Qa)*VR + cot(Ra)*VQ;
        }
        cn[Vl] = cn[Vl] / (2.0 * Amix);
    }
}

// Compare exact curvature with numerical curvature
void checkCurvature(const triSurfacePointVectorField& cn, const triSurface& front,
                    scalar radius, std::fstream& errorFile)
{
    scalar curvatureExact = 2.0 / radius;

    scalar minCurvature = curvatureExact;
    scalar maxCurvature = 0;

    scalar averageCurvature = 0;
    scalar linearDeviation = 0;
    scalar quadDeviation = 0;
    scalar maxDeviation = 0; // L_inf norm as in Francois' paper

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
    scalar maxDeviation = 0;
    scalar devAngle = 0;
    scalar linearDeviation = 0;
    scalar quadDeviation = 0;

    scalar counter = 0;

    const pointField& localPoints = front.localPoints();

    forAll(localPoints, P)
    {
        vector normalExact = center - localPoints[P];
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
    scalar averageDeviation = 0;
    scalar maxDeviation = 0;
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

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "errorFile",
        "Path to file where results will be appended"
    );

    argList::addOption
    (
        "radius",
        "Radius of the sphere used for curvature calculation"
    );

    argList::addOption
    (
        "center",
        "Center of the sphere."    
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"
    
    if (!args.optionFound("errorFile"))
    {
        FatalErrorIn("main()")
            << "Please use option '-errorFile' to name the output file"
            << endl << exit(FatalError);
    }

    if (!args.optionFound("radius"))
    {
        FatalErrorIn("main()")
            << "Please use option '-radius' to specify the sphere's radius"
            << endl << exit(FatalError);
    }

    if (!args.optionFound("center"))
    {
        FatalErrorIn("main()")
            << "Please use option '-center' to specify the sphere's center"
            << endl << exit(FatalError);
    }

    // Read command line arguments
    const scalar radius = args.optionRead<scalar>("radius");
    const vector center = args.optionRead<vector>("center");
    const std::string errorFileNameBase = args.optionRead<fileName>("errorFile");

    // Define three file names for different errors
    const std::string errorFileNameCurvature = errorFileNameBase + ".curvature";
    const std::string errorFileNameNormalVector = errorFileNameBase + ".normvec";
    const std::string errorFileNameSphereDev = errorFileNameBase + ".spheredev";

    // Open file to write results to
    const char* errorFileNameCurvaturePtr = errorFileNameCurvature.c_str();
    const char* errorFileNameNormalVectorPtr = errorFileNameNormalVector.c_str();
    const char* errorFileNameSphereDevPtr = errorFileNameSphereDev.c_str();

    std::fstream errorFileCurvature;
    std::fstream errorFileNormalVector;
    std::fstream errorFileSphereDev;

    errorFileCurvature.open(errorFileNameCurvaturePtr, std::ios_base::app);
    errorFileNormalVector.open(errorFileNameNormalVectorPtr, std::ios_base::app);
    errorFileSphereDev.open(errorFileNameSphereDevPtr, std::ios_base::app);

    // Read and intialize front 
    triSurface front(
        "front/front.stl"
    );

    // Initialize field for curvature normals
    triSurfacePointVectorField cn
    (
        IOobject
        (
            "curvatureNormals",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        front,
        dimensionedVector
        (
            "zero",
            pow(dimLength, -1),
            vector(0,0,0)
        )
    );

    // Print number of mesh points and faces
    Info << "Number of front mesh points: " << front.meshPoints().size() << endl;
    Info << "Number of front mesh triangles: " << front.localFaces().size() << endl;

    // Finally call the function
    curvatureNormals(cn, front);

    // Check deviation from sphere
    sphereDeviation(front, radius, center, errorFileSphereDev);

    // Check curvature
    checkCurvature(cn, front, radius, errorFileCurvature);

    // Check normals
    checkNormal(cn, front, center, errorFileNormalVector);

    errorFileCurvature.close();
    errorFileNormalVector.close();
    errorFileSphereDev.close();

    return 0;
}

// ************************************************************************* //
