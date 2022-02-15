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
    Foam::singleCellHarmonicMarkerFieldModel

SourceFiles
    singleCellHarmonicMarkerFieldModel.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description

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


#include "singleCellHarmonicMarkerFieldModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(singleCellHarmonicMarkerFieldModel, 0);
    addToRunTimeSelectionTable(markerFieldModel, singleCellHarmonicMarkerFieldModel, Dictionary);

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //
singleCellHarmonicMarkerFieldModel::singleCellHarmonicMarkerFieldModel(const dictionary& configDict)
    :
        sharpMarkerFieldModel(configDict), 
        pointDistFieldName_(configDict.get<word>("pointDistance"))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void singleCellHarmonicMarkerFieldModel::calcMarkerField(volScalarField& markerField) const
{
    const scalar pi = constant::mathematical::pi;
    const fvMesh& mesh = markerField.mesh(); 
    const edgeList& meshEdges = mesh.edges(); 
    const labelListList& meshEdgeCells = mesh.edgeCells();  

    const volScalarField& searchDistanceSqr = 
        mesh.lookupObject<volScalarField>(sqrSearchDistFieldName()); 

    const pointScalarField& pointDistance = 
        mesh.lookupObject<pointScalarField>(pointDistFieldName()); 

    const volScalarField& signedDistance = 
        mesh.lookupObject<volScalarField>(cellDistFieldName()); 

    // Compute the sharp marker field [0,1] based on the distance sign.
    sharpMarkerFieldModel::calcMarkerField(markerField); 

    // Find an edge of a cell where the point signed distance switches in sign.
    // Find the cells that share that edge and set the marker field using the
    // harmonic model. 
    forAll(meshEdges, edgeI)
    {
        const edge& meshEdge = meshEdges[edgeI]; 
        auto firstPointDist = pointDistance[meshEdge[0]];
        auto secondPointDist = pointDistance[meshEdge[1]];

        // If the distance switches sign for an edge of the cell cellI. 
        if ((firstPointDist * secondPointDist ) < 0)
        {
            const labelList& edgeCells = meshEdgeCells[edgeI]; 
            forAll(edgeCells, edgeJ)
            {
                const label cellLabel = edgeCells[edgeJ]; 
                scalar searchDistance = sqrt(searchDistanceSqr[cellLabel]); 

                markerField[cellLabel] = 0.5 * (
                    1 + signedDistance[cellLabel] / searchDistance + 1/pi *
                    sin((pi * signedDistance[cellLabel]) / searchDistance)
                );
            }
        }
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

