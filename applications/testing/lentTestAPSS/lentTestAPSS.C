/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    lentTestAPSS

Description
    For now: store normals and curvature in separate fields

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "triSurfaceFields.H"
#include "interfaceProperties.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"

#include "lentMethod.H"

#include <fstream>

// Include Eigen library for basic linear algebra operations
#include "Eigen/Core"
#include "Eigen/Dense"

#include "auxFunctions/auxFunctions.H"
#include "auxFunctions/errorMetrics.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Typedefs for more concise use of Eigen
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;
typedef Eigen::Matrix<double, 5, 5> Matrix5d;
typedef Eigen::Matrix<double, 5, 1> Vector5d;

// Weighting function from APSS; remember that the argument x is the 
// distance squared divided by h squared
scalar weight(scalar x2byh2)
{
    if (x2byh2 < 1.0)
    {
        return (1.0 - x2byh2)*(1.0 - x2byh2)*(1.0 - x2byh2)*(1.0 - x2byh2);
    }
    else
    {
        return 0.0;
    }
}

// Computation of the feature support size h; the most viable approach is not
// clear yet. h should be chosen so sufficient points are included to
// achieve a little smoothing/averaging without erasing features.
//
// h should be related to the local size of the Eulerian mesh, not
// to front as its elements are of highly different sizes 
// For now, hard code the value for the reconstruction 32*32*32
// test case: h = 2.0*1/32
scalar supportSize()
{
    scalar h = 0.06;
    return h;
}


// Return true if label check is not yet in list elements
// Used do detect new unique points
bool isNotInList(const DynamicList<label,1,1,1>& elements, label check)
{
    bool notInList = true;

    forAll(elements, N)
    {
        if (check == elements[N]) {notInList = false;}
    }

    return notInList;
}

// Get unique new neighbors for points listed in center;
// points listed in known points are ignored
DynamicList<label,1,1,1> getNeighborPoints(const DynamicList<label,1,1,1>&
                                                 centers,
                                           const DynamicList<label,1,1,1>&
                                                 knownPoints,
                                           const labelListList& adjacentEdges,
                                           const edgeList& edges)
{
    DynamicList<label,1,1,1> neighborPoints(1);

    forAll(centers, Cl)
    {
        const label curPoint = centers[Cl];
        const labelList& edgeStar = adjacentEdges[curPoint];

        forAll(edgeStar, El)
        {
            edge E = edges[edgeStar[El]];

            label newVertex = E.otherVertex(curPoint);

            if (isNotInList(knownPoints, newVertex) &&
                isNotInList(neighborPoints, newVertex))
            {
                neighborPoints.append(newVertex);
            }
        }
    }

    return neighborPoints;
}

// Collect topological neighbors of a given reference point;
// only the topological distance in terms of number of edges is
// considered, but no geometrical information. Geometric filtering
// is performed by the weight function
DynamicList<label,1,1,1> findSupportPoints(const label refPoint,
                                           const triSurface& front,
                                           label numRings)
{
    // Note: numRings describes up to which neighborhood potential support
    // points are included.
    // 0 means only the reference point itself is included, 1 returns the
    // 1-ring neighborhood and so forth

    DynamicList<label,1,1,1> knownPoints(1);
    DynamicList<label,1,1,1> neighborhoodCenters(1);
    neighborhoodCenters.append(refPoint);
    
    // Get point-edge assignment
    const labelListList& adjacentEdges = front.pointEdges();

    // Get list of edges from which point labels are retrieved
    const edgeList& edges = front.edges();

    for (int i = 0; i <= numRings; i++)
    {
        neighborhoodCenters = getNeighborPoints(neighborhoodCenters,
                                                knownPoints,
                                                adjacentEdges,
                                                edges);

        knownPoints.append(neighborhoodCenters);
    }

    return knownPoints;
}

// Compute normal at a front vertex using area weighted triangle
// normals
void computeFrontVertexNormals(triSurfacePointVectorField& normals,
                               const triSurface& front)
{
    normals = vector(0.0, 0.0, 0.0);
                
    // Get label list of all mesh vertices
    const labelList& vertices = front.meshPoints();

    // Get assignment point -> faces
    const labelListList& adjacentFaces = front.pointFaces();

    // List of local point references for current patch
    const pointField& localPoints = front.localPoints();

    // List of local faces references
    const List<labelledTri>& localFaces = front.localFaces();

    forAll(vertices, Vl)
    {
        const labelList& oneRingNeighborhood = adjacentFaces[Vl];

        forAll(oneRingNeighborhood, T)
        {
            label triLabel = oneRingNeighborhood[T];
            labelledTri currentTri = localFaces[triLabel];

            // Resemble point notion from lentFoam paper
            point VT0 = localPoints[currentTri[0]];
            point VT1 = localPoints[currentTri[1]];
            point VT2 = localPoints[currentTri[2]];

            // Compute area normal according to lentFoam paper
            vector areaNormal = (VT1 - VT0) ^ (VT2 - VT0);

            normals[Vl] += areaNormal;
        }

        normals[Vl] = normals[Vl] / mag(normals[Vl]);
    }  
}

// Perform matrix and vector setup in functions for better readabilty of the
// method itself
void setUpW(MatrixXd& W, const DynamicList<label,1,1,1>& neighborPoints,
            const pointField& localPoints, const point& refPoint)
{
    W.setZero();

    label numPoints = neighborPoints.size();
    scalar h = supportSize();
        
    forAll(neighborPoints, Pl)
    {
        point P = localPoints[neighborPoints[Pl]];
        W(Pl, Pl) = weight(magSqr(P - refPoint)/(h*h));
    }
    // NOTE: changing the definition of the feature support size h may
    // require definition of beta in a separate function
    scalar beta = 1e6*h*h;

    for (label counter = 0; counter < numPoints; counter++)
    {
        W(numPoints+counter, numPoints+counter) = beta*W(counter, counter);
        W(2*numPoints+counter, 2*numPoints+counter) = beta*W(counter, counter);
        W(3*numPoints+counter, 3*numPoints+counter) = beta*W(counter, counter);
    }
}

void setUpD(MatrixXd& D, const DynamicList<label,1,1,1>& neighborPoints,
            const pointField& localPoints)
{
        D.setZero();

        label numPoints = neighborPoints.size();

        vector e0(1, 0, 0); 
        vector e1(0, 1, 0); 
        vector e2(0, 0, 1); 
        
        forAll(neighborPoints, Pl)
        {
            point P = localPoints[Pl];

            D(Pl, 0) = 1;
            D(Pl, 1) = P[0];
            D(Pl, 2) = P[1];
            D(Pl, 3) = P[2];
            D(Pl, 4) = P & P;

            D(numPoints+Pl, 1) = e0[0];
            D(numPoints+Pl, 2) = e0[1];
            D(numPoints+Pl, 3) = e0[2];
            D(numPoints+Pl, 4) = 2*e0 & P;

            D(2*numPoints+Pl, 1) = e1[0];
            D(2*numPoints+Pl, 2) = e1[1];
            D(2*numPoints+Pl, 3) = e1[2];
            D(2*numPoints+Pl, 4) = 2*e1 & P;

            D(3*numPoints+Pl, 1) = e2[0];
            D(3*numPoints+Pl, 2) = e2[1];
            D(3*numPoints+Pl, 3) = e2[2];
            D(3*numPoints+Pl, 4) = 2*e2 & P;
        }        
}

void setUpb(VectorXd& b, const DynamicList<label,1,1,1>& neighborPoints,
            const triSurfacePointVectorField& normals)
{
        b.setZero();

        label numPoints = neighborPoints.size();

        vector e0(1, 0, 0); 
        vector e1(0, 1, 0); 
        vector e2(0, 0, 1); 

        forAll(neighborPoints, Pl)
        {
            vector n = normals[Pl];

            b(numPoints+Pl) = e0 & n;
            b(2*numPoints+Pl) = e1 & n;
            b(3*numPoints+Pl) = e2 & n;
        }
}

/******************************************************************************\
 Actual method from paper "Algebraic point set surfaces"
\******************************************************************************/
void computeCurvature(triSurfacePointScalarField& curvature,
                      const triSurface& front,
                      triSurfacePointVectorField& normals)
{
    curvature = dimensionedScalar("zero",
                                  dimless/dimLength,
                                  0.0
                                 );
                
    // Get label list of all mesh vertices
    const labelList& vertices = front.meshPoints();

    // List of local point references for current patch
    const pointField& localPoints = front.localPoints();

    forAll(vertices, Vl)
    {
        point V = localPoints[Vl];

        DynamicList<label,1,1,1> neighborPoints = findSupportPoints(Vl, front, 3);
        label numPoints = neighborPoints.size();

        // Allocate matrices and vectors
        MatrixXd W(4*numPoints, 4*numPoints);
        MatrixXd D(4*numPoints, 5);
        VectorXd b(4*numPoints);

        // Set up matrices and right hand side vector
        setUpW(W, neighborPoints, localPoints, V);
        setUpD(D, neighborPoints, localPoints);
        setUpb(b, neighborPoints, normals);

        // Obtain solution vector u
        Matrix5d A = D.transpose() * W * D;
        Vector5d bhat = D.transpose() * W * b;
        Vector5d u = A.colPivHouseholderQr().solve(bhat);

        vector c(u(1), u(2), u(3));
        c = c/(-2*u(4));

        scalar t = u(0)/u(4);

        curvature[Vl] = 2 / (Foam::sqrt(magSqr(c) - t));
    }

}


// Main program:
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

    argList::addOption
    (
        "reconTimes",
        "Number of times the front is to be reconstructed before curvature calculation."
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

    if (!args.optionFound("reconTimes"))
    {
        FatalErrorIn("main")
            << "Please use option 'reconTimes' to specify the number of "
            << "reconstructions."
            << endl << exit(FatalError);
    }

    // Read command line arguments
    const scalar radius = args.optionRead<scalar>("radius");
    const vector center = args.optionRead<vector>("center");
    const std::string errorFileNameBase = args.optionRead<fileName>("errorFile");
    const label reconTimes = args.optionRead<label>("reconTimes");

    // Respect number of time steps defined in lent reconstruction trest cases
    // for now
    if (reconTimes < 0 || reconTimes > 3)
    {
        FatalErrorIn("main")
            << "Option -reconTimes is out of range. Please use n=0...4"
            << endl << exit(FatalError);
    }

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
    // Get correct file name first
    std::string frontFileName = "";
    switch (reconTimes)
    {
        case 0:
            frontFileName = "front/front.stl";
            break;
        case 1:
            frontFileName = "front/front-00000000.vtk";
            break;
        case 2:
            frontFileName = "front/front-00000001.vtk";
            break;
        case 3:
            frontFileName = "front/front-00000002.vtk";
            break;
        case 4:
            frontFileName = "front/front-00000003.vtk";
            break;
    }

    Info << "Loading front" << endl;
    triSurface front(
        frontFileName
    );

    Info << "Front loaded\n" << endl;

    // Initialize field for normals
    triSurfacePointVectorField normals
    (
        IOobject
        (
            "normals",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        front,
        vector(0.0,0.0,0.0)
    );

    // Initialize field for curvature
    triSurfacePointScalarField curvature
    (
        IOobject
        (
            "curvature",
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

    // Initialize field for curvatureNormals
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
            vector(0.0, 0.0, 0.0)
        )
    );

    // Print number of mesh points and faces
    Info << "Number of front mesh points: " << front.meshPoints().size() << endl;
    Info << "Number of front mesh triangles: " << front.localFaces().size() << endl;

    // Finally call the functions
    computeFrontVertexNormals(normals, front);
    computeCurvature(curvature, front, normals);

    // Fuse normals and curvature into single field
    cn = normals * curvature;

    // Check deviation from sphere
    meshQuality(front, errorFileSphereDev);
    sphereDeviation(front, radius, center, errorFileSphereDev);

    // Check curvature
    meshQuality(front, errorFileCurvature);
    checkCurvature(cn, front, radius, errorFileCurvature);

    // Check normals
    meshQuality(front, errorFileNormalVector);
    checkNormal(normals, front, center, errorFileNormalVector);

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
            dimless,
            0.0
        )
    );

    dimensionedScalar dradius = dimensionedScalar("zero", dimLength, radius);

    // Compute difference to exact curvature and relate the difference to
    // the exact curvature
    curvatureError = mag(mag(cn) - 2.0/dradius)*dradius/2.0;

    cn.write();
    curvatureError.write();

    errorFileCurvature.close();
    errorFileNormalVector.close();
    errorFileSphereDev.close();
    

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
