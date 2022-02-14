/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
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

Class
    Foam::meyerCurvatureModel

SourceFiles
    meyerCurvatureModel.C

Authors:
    Tobias Tolle   (tolle@mma.tu-darmstadt.de)
    Tomislav Maric (tolle@mma.tu-darmstadt.de)
 
Affiliations:
    Mathematical Modeling and Analysis Institute, Mathematics Department, 
    TU Darmstadt, Germany

Funding:
    German Research Foundation (DFG) - Project-ID 265191195 - SFB 1194

    German Research Foundation (DFG) - Project-ID MA 8465/1-1, 
    Initiation of International Collaboration 
    "Hybrid Level Set / Front Tracking methods for simulating 
    multiphase flows in geometrically complex systems"

\*---------------------------------------------------------------------------*/

#include "meyerCurvatureModel.H"
#include "addToRunTimeSelectionTable.H"
#include "triSurfaceFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(meyerCurvatureModel, 0);
    addToRunTimeSelectionTable(curvatureModel, meyerCurvatureModel, Dictionary);
    
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void meyerCurvatureModel::computeCurvature(const fvMesh& mesh, const triSurfaceFront& frontMesh) const
{
    triSurfaceFrontPointVectorField cn{
                    IOobject(
                        "curvature_normal", 
                        mesh.time().timeName(), 
                        frontMesh,
                        IOobject::NO_READ, 
                        IOobject::AUTO_WRITE
                    ), 
                    frontMesh, 
                    dimensionedVector(
                        "zero", 
                        dimless/dimLength, 
                        vector(0.0,0.0,0.0)
                    )
                };

    // TODO: modify so that mutliple patches, e.g. in the case of
    // bubble breakup, are supported

    // All of the following functions are taken from
    // "PrimitivePatch.H"

    // Get assignment point --> faces
    const labelListList& adjacentFaces = frontMesh.pointFaces();

    // List of local point references for current patch
    const pointField& vertices = frontMesh.localPoints();
    
    // Needed: List of local face references
    const List<labelledTri>& localFaces = frontMesh.localFaces();

    const auto& faceAreas = frontMesh.magSf();

    // V (the actual vertex) and Vl (its label) are used synonymously
    // in the following comments
    forAll(vertices, Vl)
    {
        scalar Amix = 0.0;

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
            labelList triVertices = orderVertices(currentTri, Vl);

            point V = vertices[triVertices[0]];
            point Q = vertices[triVertices[1]];
            point R = vertices[triVertices[2]];

            // Edge vectors
            vector VQ = Q - V;
            vector VR = R - V;
            vector QR = R - Q;

            // Get angles
            scalar Va = getAngle(VQ, VR);
            scalar Ra = getAngle(VR, QR);
            scalar Qa = M_PI - (Va + Ra);

            // Check if non-obtuse in order to use the correct area metric
            if (Va < 0.5*M_PI && Qa < 0.5*M_PI && Ra < 0.5*M_PI)
            {
                // Use Voronoi-area
                // Cotangent function 'cot' has to be defined locally since
                // it is not offered by OpenFOAM
                Amix += 0.125*(magSqr(VR)*cot(-VQ, QR) + magSqr(VQ)*cot(VR, QR));    
            }
            else if (Va >= 0.5*M_PI)
            {
                // Obtuse angle at V, use half area
                Amix += 0.5 * faceAreas[triLabel];
            }
            else
            {
                // Obtuse angle not at V, use quarter area
                Amix += 0.25 * faceAreas[triLabel];
            }

            // Compute mean curvature normal contributions
            cn[Vl] += cot(-VQ, QR)*VR + cot(VR, QR)*VQ;
        }
        cn[Vl] = cn[Vl] / (2.0 * Amix);
    }

    // ---------------------------------------------------------------
    // FIXME: currently only values located at triangles - not at vertices -
    // can be transfered to the Eulerian mesh. Thus, average the vertex values
    // for each triangle and then perform the transfer
    auto& cnTria = *curvatureBuffer(frontMesh);
    forAll(cnTria, I)
    {
        const auto& aFace = localFaces[I];

        for (const auto& vertexID : aFace)
        {
            cnTria[I] += cn[vertexID];
        }
    }

    cnTria /= 3.0;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meyerCurvatureModel::meyerCurvatureModel(const dictionary& configDict) 
    :
        frontBasedCurvatureModel{configDict}
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

labelList meyerCurvatureModel::orderVertices(labelledTri& tri, label V) const
{
    labelList vertices(3);

    if (tri[0] == V)
    {
        vertices[0] = tri[0];
        vertices[1] = tri[1];
        vertices[2] = tri[2];
    }
    else if (tri[1] == V)
    {
        vertices[0] = tri[1];
        vertices[1] = tri[0];
        vertices[2] = tri[2];
    }
    else
    {
        vertices[0] = tri[2];
        vertices[1] = tri[0];
        vertices[2] = tri[1];
    }

    return vertices;
}


scalar meyerCurvatureModel::getAngle(vector& a, vector& b) const
{
    return Foam::acos((a & b) / (mag(a) * mag(b)));
}


// Cotangent function to resemble algorithm notion from paper
// Removed trigonometric functions to prevent precision loss. TT
scalar meyerCurvatureModel::cot(vector a, vector b) const
{
    return a & b / mag(a ^ b);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
