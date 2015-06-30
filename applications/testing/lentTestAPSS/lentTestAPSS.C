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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// TODO: include error metrics from curvatureNormal

// Typedefs for more concise use of Eigen
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;

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

    h = (h / neighbors.size()) * scaleNeighborhoodSize;

    return h;
}

// Return the neighbors of a given point 'center' as a label list
labelList getNeighborPoints(const label center,
                            const labelList& oneRingNeighborhood,
                            const List<labelledTri>& localFaces)
{
    // Number of trinagles in one ring neighborhood matches the number
    // of unique points different from the center point
    label size = oneRingNeighborhood.size();
    labelList neighbors(size);

    // Set to negative value so the check for double entries starts
    // work after the first triangle has been looked upon
    neighbors[0] = -1;

    // Only inspect every second triangle --> less effort to find unique
    // point labels
    label even = 0;
    label index = 0;

    forAll(oneRingNeighborhood, T)
    {
        if (even%2 == 0)
        {
            label triLabel = oneRingNeighborhood[T];
            labelledTri currentTri = localFaces[triLabel];

            forAll(currentTri, C)
            {
                label currentLabel = currentTri[C];

                if (currentLabel != center && currentLabel != neighbors[0])
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
        scalar area = 0.0;

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
            area += mag(areaNormal);
        }

        normals[Vl] = normals[Vl] / area;
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

        const label numPoints = neighborPoints.size() + 1;

        scalar h = supportSize(V, neighborPoints, localPoints);
        
        // For first tests, skip search for further points which may lay in
        // the feature sphere

        // Set up matrices and vectors for 
        MatrixXd W(4*numPoints, 4*numPoints);
        MatrixXd D(4*numPoints, 5);
        VectorXd b(4*numPoints);

        W.setZero();
        D.setZero();
        b.setZero();

        // Set up W
        W(0,0) = 1;

        label counter = 1;
        forAll(neighborPoints, Pl)
        {
            point P = localPoints[Pl];
            W(counter, counter) = weight(magSqr(P - V)/(h*h));
            counter++;
        }

        // NOTE: changing the definition of the feature support size h may
        // require definition of beta in a separate function
        scalar beta = 1e6*h*h;

        // 1. Compute feature size h
        // 2. Determine points located in sphere with center Vl and radius h
        // 3. Set up matrices W, D and b
    }

}


// Main program:
int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
