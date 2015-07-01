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
// TODO: include error metrics from curvatureNormal

// Typedefs for more concise use of Eigen
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;
typedef Eigen::Matrix<double, 5, 5> Matrix5d;
typedef Eigen::Matrix<double, 5, 1> Vector5d;

// Aux functions
// Weighting function from APSS; remember that the argument x is the 
// distance squared divided by h squared
scalar weight(scalar x2byh2)
{
    scalar result = (1.0 - x2byh2)*(1.0 - x2byh2)*(1.0 - x2byh2)*(1.0 - x2byh2);

    if (result < 0) {result = 0;}

    return result;
}

// Computation of the feature support size h; the most viable approach is not
// clear yet. h should be chosen so sufficient points are included to
// achieve a little smoothing/averaging without erasing features
// A relation to the cell size of the Eulerian grid or the one/two ring
// neighborhood seem reasonable starting points
//
// For now: compute average distance between point at hand and its
// one ring neighbors; multiply this with an arbitrarily chosen factor (e.g. 1.5)
// so that at least the direct neighbors are assigned a weight larger than zero.
scalar supportSize(point center, const labelList& neighbors,
                   const pointField& localPoints)
{
    scalar h = 0;
    scalar scaleNeighborhoodSize = 1.5;

    forAll(neighbors, N)
    {
        point neighbor = localPoints[neighbors[N]];

        h += mag(center - neighbor);
    }

    // Exclude center point from number of points
    h = (h / (neighbors.size() - 1)) * scaleNeighborhoodSize;

    return h;
}

// Return the neighbors of a given point 'center' as a label list
// It is useful to include the reference point itself for a more
// consistent set up of the matrices in the course of algebraic
// sphere fitting
labelList getNeighborPoints(const label center,
                            const labelList& oneRingNeighborhood,
                            const List<labelledTri>& localFaces)
{
    // Number of trinagles in one ring neighborhood matches the number
    // of unique points different from the center point
    label size = oneRingNeighborhood.size() + 1;
    labelList neighbors(size);

    // Set to negative value so the check for double entries starts
    // work after the first triangle has been looked upon
    neighbors[0] = center;

    // Only inspect every second triangle --> less effort to find unique
    // point labels
    label even = 0;
    label index = 1;

    forAll(oneRingNeighborhood, T)
    {
        if (even%2 == 0)
        {
            label triLabel = oneRingNeighborhood[T];
            labelledTri currentTri = localFaces[triLabel];

            forAll(currentTri, C)
            {
                label currentLabel = currentTri[C];

                if (currentLabel != center && currentLabel != neighbors[1])
                {
                    neighbors[index] = currentLabel;
                    index++;
                }
            }   
        }

        even++;
    }

    return neighbors;
}


// Compute normal at a front vertex using area weighted triangle
// normals
// NOTE: seems that this is already implemented in OpenFOAM:
// --> const Field<point>& normals =  triSurface.pointNormals()
// --> Check this
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

    // Get assignment point -> faces
    const labelListList& adjacentFaces = front.pointFaces();

    // List of local point references for current patch
    const pointField& localPoints = front.localPoints();

    // List of local faces references
    const List<labelledTri>& localFaces = front.localFaces();

    forAll(vertices, Vl)
    {
        const labelList& oneRingNeighborhood = adjacentFaces[Vl];
        point V = localPoints[Vl];

        const labelList neighborPoints = getNeighborPoints(Vl,
                                                oneRingNeighborhood,
                                                localFaces);

        // + 1 since reference point is not included in neighborPoints
        const label numPoints = neighborPoints.size() + 1;

        scalar h = supportSize(V, neighborPoints, localPoints);
        
        // For first tests, skip search for further points which may lay in
        // the feature sphere

        // Set up matrices and vectors
        MatrixXd W(4*numPoints, 4*numPoints);
        MatrixXd D(4*numPoints, 5);
        VectorXd b(4*numPoints);

        W.setZero();
        D.setZero();
        b.setZero();

        // Set up W
        label counter = 0;
        forAll(neighborPoints, Pl)
        {
            point P = localPoints[Pl];
            W(counter, counter) = weight(magSqr(P - V)/(h*h));
            counter++;
        }

        // NOTE: changing the definition of the feature support size h may
        // require definition of beta in a separate function
        scalar beta = 1e6*h*h;

        for (counter = 0; counter < numPoints; counter++)
        {
            W(numPoints+counter, numPoints+counter) = beta*W(counter, counter);
            W(2*numPoints+counter, 2*numPoints+counter) = beta*W(counter, counter);
            W(3*numPoints+counter, 3*numPoints+counter) = beta*W(counter, counter);
        }

        // Set up D
        vector e0(1, 0, 0); 
        vector e1(0, 1, 0); 
        vector e2(0, 0, 1); 
        
        counter = 0;
        forAll(neighborPoints, Pl)
        {
            point P = localPoints[Pl];

            D(counter, 0) = 1;
            D(counter, 1) = P[0];
            D(counter, 2) = P[1];
            D(counter, 3) = P[2];
            D(counter, 4) = P & P;

            D(numPoints+counter, 1) = e0[0];
            D(numPoints+counter, 2) = e0[1];
            D(numPoints+counter, 3) = e0[2];
            D(numPoints+counter, 4) = 2*e0 & P;

            D(2*numPoints+counter, 1) = e1[0];
            D(2*numPoints+counter, 2) = e1[1];
            D(2*numPoints+counter, 3) = e1[2];
            D(2*numPoints+counter, 4) = 2*e1 & P;

            D(3*numPoints+counter, 1) = e2[0];
            D(3*numPoints+counter, 2) = e2[1];
            D(3*numPoints+counter, 3) = e2[2];
            D(3*numPoints+counter, 4) = 2*e2 & P;

            counter++;
        }        

        // Set up b
        counter = 0;
        forAll(neighborPoints, Pl)
        {
            vector n = normals[Pl];

            b(numPoints+counter) = e0 & n;
            b(2*numPoints+counter) = e1 & n;
            b(3*numPoints+counter) = e2 & n;

            counter++;
        }

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

    // Print number of mesh points and faces
    Info << "Number of front mesh points: " << front.meshPoints().size() << endl;
    Info << "Number of front mesh triangles: " << front.localFaces().size() << endl;

    // Finally call the functions
    computeFrontVertexNormals(normals, front);
    computeCurvature(curvature, front, normals);

    // Fuse normals and curvature into single field
    normals = normals * curvature;

    // Check deviation from sphere
    meshQuality(front, errorFileSphereDev);
    sphereDeviation(front, radius, center, errorFileSphereDev);

    // Check curvature
    meshQuality(front, errorFileCurvature);
    checkCurvature(normals, front, radius, errorFileCurvature);

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
    curvatureError = mag(mag(normals) - 2.0/dradius)*dradius/2.0;

    normals.write();
    curvatureError.write();

    errorFileCurvature.close();
    errorFileNormalVector.close();
    errorFileSphereDev.close();
    

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
