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

// Determine mesh quality
void meshQuality(const triSurface& front, std::fstream& file)
{
    scalar pi = 3.141592653589793238; 

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

// Aux function to get a reference point for consistent normal vectors
// Does only work for convex, closed surfaces
point getRefPoint(const triSurface& front)
{
    point refPoint = vector(0,0,0);

    // Approach: get extremal points for each axis (x,y,z), thereby putting
    // the front in a box. The center of this cuboid should also be located
    // inside the front (assuming the surface is convex...)
    const labelList& vertices = front.meshPoints();
    const pointField& localPoints = front.localPoints();

    point initValue = localPoints[0];

    scalar xmin = initValue[0];
    scalar xmax = initValue[0];
    scalar ymin = initValue[1];
    scalar ymax = initValue[1];
    scalar zmin = initValue[2];
    scalar zmax = initValue[2];

    forAll(vertices, Vl)
    {
        point V = localPoints[Vl];

        xmin = V[0] < xmin ? V[0] : xmin;
        xmax = V[0] > xmax ? V[0] : xmax;

        ymin = V[1] < ymin ? V[1] : ymin;
        ymax = V[1] > ymax ? V[1] : ymax;

        zmin = V[2] < zmin ? V[2] : zmin;
        zmax = V[2] > zmax ? V[2] : zmax;
    }

    // Assemble reference point
    refPoint[0] = 0.5 * (xmin+xmax);
    refPoint[1] = 0.5 * (ymin+ymax);
    refPoint[2] = 0.5 * (zmin+zmax);

    return refPoint;
}

// Method taken from tryggvason book "Direct numerical simulations..."
// Despite 
void noCurvature(triSurfacePointVectorField& cn, const triSurface& front)
{
    cn = dimensionedVector("zero",
                            dimless/dimLength,
                            vector(0,0,0)
                           );

    // Get necessary references
    const labelList& vertices = front.meshPoints();
    const labelListList& adjacentFaces = front.pointFaces();
    const pointField& localPoints = front.localPoints();
    const List<labelledTri>& localFaces = front.localFaces();

    // Get reference point inside surface
    point refPoint = getRefPoint(front);

    // Iterate over all vertices instead of triangles
    // --> easier consistency check and normalization
    // back to K*n
    forAll(vertices, Vl)
    {
        // Get all triangles adjacent to V
        const labelList& oneRingNeighborhood = adjacentFaces[Vl];

        forAll(oneRingNeighborhood, T)
        {
            labelledTri tri = localFaces[T];

            label tri0 = tri[0];
            label tri1 = tri[1];
            label tri2 = tri[2];

            // Triangle vertices
            point x1;
            point x2;
            point x3;

            if (tri0 == Vl)
            {
                x1 = localPoints[Vl];
                x2 = localPoints[tri1];
                x3 = localPoints[tri2];
            }
            else if (tri1 == Vl)
            {
                x1 = localPoints[Vl];
                x2 = localPoints[tri0];
                x3 = localPoints[tri2];
            }
            else
            {
                x1 = localPoints[Vl];
                x2 = localPoints[tri0];
                x3 = localPoints[tri1];
            }

            // Triangle edges
            vector x13 = x3 - x1;
            vector x21 = x1 - x2;
            vector x32 = x2 - x3;

            // Midpoint of edge opposite to V
            // Needed to ensure that the "pull" points outwards
            point midpoint = 0.5*(x2+x3);
            vector direction = x1 - midpoint;

            // Get consistent normal vector of face
            // The folllowing method only works for completlely convex surfaces
            // which are closed
            vector faceNormal = (x13^x21) / mag(x13^x21);
            
            // Consistency check
            scalar angle = (x1-refPoint) & faceNormal;
            if (angle < 0.0)
            {
                faceNormal = -faceNormal;
            }

            // Curvature normal contribution
            angle = (faceNormal ^ x32) & direction;
            if (angle > 0.0)
            {
                cn[Vl] += 0.5*(faceNormal ^ x32);
            }
            else
            {
                cn[Vl] += 0.5*(x32 ^ faceNormal);
            }
        }

        // Divide by area to finally get curvature normal
        cn[Vl] = cn[Vl];
    }
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
    noCurvature(cn, front);

    // Check deviation from sphere
    meshQuality(front, errorFileSphereDev);
    sphereDeviation(front, radius, center, errorFileSphereDev);

    // Check curvature
    meshQuality(front, errorFileCurvature);
    checkCurvature(cn, front, radius, errorFileCurvature);

    // Check normals
    meshQuality(front, errorFileNormalVector);
    checkNormal(cn, front, center, errorFileNormalVector);

    // Write curvature normals field for visual inspection
    // Additionally, write scalar curvature error field
    triSurfacePointScalarField curvatureError
    (
        IOobject
        (
            "curvatureErrors",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        front,
        dimensionedScalar
        (
            "zero",
            pow(dimLength, -1),
            0
        )
    );

    dimensionedScalar dradius = dimensionedScalar("zero", dimLength, radius);
    curvatureError = mag(mag(cn) - 2.0/dradius);

    cn.write();
    curvatureError.write();

    errorFileCurvature.close();
    errorFileNormalVector.close();
    errorFileSphereDev.close();

    return 0;
}

// ************************************************************************* //
