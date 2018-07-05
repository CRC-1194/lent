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

Class
    Foam::frontCurvatureMeyer

SourceFiles
    frontCurvatureMeyer.C

Author
    Tobias Tolle bt@lefou-familie.org

Contributors
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Interface for the front curvature models. 

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

#include "frontCurvatureMeyer.H"
#include "addToRunTimeSelectionTable.H"
#include "triSurfaceFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontCurvatureMeyer, 0);
    addToRunTimeSelectionTable(frontCurvatureModel, frontCurvatureMeyer, Dictionary);
    
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void frontCurvatureMeyer::initializeCurvatureNormal(const fvMesh& mesh, const triSurfaceFront& front) const
{
    const Time& runTime = mesh.time();  

    if (curvatureNormalTmp_.empty())
    {
        curvatureNormalTmp_ =
            tmp<triSurfaceFrontPointVectorField>
            {
                new triSurfaceFrontPointVectorField
                (
                    IOobject(
                        "curvature_normal", 
                        runTime.timeName(), 
                        front,
                        IOobject::NO_READ, 
                        IOobject::AUTO_WRITE
                    ), 
                    front, 
                    dimensionedVector(
                        "zero", 
                        dimless/dimLength, 
                        vector(0.0,0.0,0.0)
                    )
                )
            };

        return;
    }
    else if (curvatureNormalTmp_->size() != front.localPoints().size())
    {
        curvatureNormalTmp_->resize(front.localPoints().size());
    }

    for (auto& cn : curvatureNormalTmp_.ref())
    {
        cn = vector{0,0,0};
    }
}

void frontCurvatureMeyer::computeCurvature(const fvMesh& mesh, const triSurfaceFront& frontMesh) const
{
    initializeCurvatureNormal(mesh, frontMesh);

    auto& cn = curvatureNormalTmp_.ref();

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
    triSurfaceFrontVectorField cnTria
    {
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
                        dimless, 
                        vector(0.0,0.0,0.0)
                    )
    };

    // FIXME: currently only values located at triangles - not at vertices -
    // can be transfered to the Eulerian mesh. Thus, average the vertex values
    // for each triangle and then perform the transfer
    forAll(cnTria, I)
    {
        const auto& aFace = localFaces[I];

        for (const auto& vertexID : aFace)
        {
            cnTria[I] += cn[vertexID];
        }
    }

    cnTria /= 3.0;

    const auto& frontToMesh = frontToMeshTmp_.ref();
    frontToMesh.transferCurvature(cnTria, frontMesh, mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontCurvatureMeyer::frontCurvatureMeyer(const dictionary& configDict) 
    :
        frontCurvatureModel{configDict},
        curvatureNormalTmp_{},
        frontToMeshTmp_{
            frontToMeshTransferModel::New(configDict.subDict("frontToMeshTransfer"))
        }
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

labelList frontCurvatureMeyer::orderVertices(labelledTri& tri, label V) const
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


scalar frontCurvatureMeyer::getAngle(vector& a, vector& b) const
{
    return Foam::acos((a & b) / (mag(a) * mag(b)));
}


// Cotangent function to resemble algorithm notion from paper
// Removed trigonometric functions to prevent precision loss. TT
scalar frontCurvatureMeyer::cot(vector a, vector b) const
{
    return a & b / mag(a ^ b);
}

tmp<volScalarField> frontCurvatureMeyer::cellCurvature(
    const fvMesh& mesh,
    const triSurfaceFront& frontMesh
) const
{
    if (curvatureNeedsUpdate(mesh))
    {
        computeCurvature(mesh, frontMesh);
        curvatureUpdated(mesh); 
    }

    const auto& frontToMesh = frontToMeshTmp_.ref();
    return frontToMesh.cellCurvature();
}

tmp<surfaceScalarField> frontCurvatureMeyer::faceCurvature(
    const fvMesh& mesh, 
    const triSurfaceFront& frontMesh
) const
{
    if (curvatureNeedsUpdate(mesh))
    {
        computeCurvature(mesh, frontMesh);
        curvatureUpdated(mesh); 
    }

    const auto& frontToMesh = frontToMeshTmp_.ref();
    return frontToMesh.faceCurvature();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
