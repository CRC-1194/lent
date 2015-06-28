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
//#include "frontCurvatureMeyer.H"

#include <fstream>

#include "auxFunctions/auxFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// TODO: include error metrics from curvatureNormal

// Aux functions


// Compute normal at a front vertex using barycentric area weighted triangle
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
                                  dimless/dimlength,
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

    }

}

// Weighting function from APSS; remember that the argument x is the 
// distance squared divided by h squared
scalar weight(scalar x2byh2)
{
    scalar result = (1.0 - x2byh2)*(1.0 - x2byh2)*(1.0 - x2byh2)*(1.0 - x2byh2);

    if (result < 0) {result = 0};

    return result;
}

// Computation of the feature support size h; the most viable approach is not
// clear yet. h should be chosen so sufficient points are included to
// achieve a little smoothing/averaging without erasing features
// A relation to the cell size of the Eulerian grid or the one/two ring
// neighborhood seem reasonable starting points
scalar supportSize(...)
{
    scalar h = 0;

    return h;
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
