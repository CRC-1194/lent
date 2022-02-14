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
    Foam::twoCellHarmonicMarkerFieldModel

SourceFiles
    twoCellHarmonicMarkerFieldModel.C

Author
    Tomislav Maric (maric@mma.tu-darmstadt.de)

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


#include "twoCellHarmonicMarkerFieldModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(twoCellHarmonicMarkerFieldModel, 0);
    addToRunTimeSelectionTable(markerFieldModel, twoCellHarmonicMarkerFieldModel, Dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void twoCellHarmonicMarkerFieldModel::calcMarkerField(volScalarField& markerField) const
{
    const scalar pi = constant::mathematical::pi;
    const fvMesh& mesh = markerField.mesh(); 
    const cellList& cells = mesh.cells(); 
    const labelList& own = mesh.owner(); 
    const labelList& nei = mesh.neighbour(); 

    const volScalarField& signedDistance = 
        mesh.lookupObject<volScalarField>(cellDistFieldName()); 

    const volScalarField& searchDistanceSqr = 
        mesh.lookupObject<volScalarField>(sqrSearchDistFieldName()); 

    // For all cell centered values of signedDistance. 
    forAll(signedDistance, cellI)
    {
        // If the cell centered signed distance switches signs across at least
        // one cell face. 
        bool signedDistanceSwitches = false; 
        const cell& meshCell = cells[cellI];
        forAll(meshCell, faceI)
        {
            const label ownCell = own[meshCell[faceI]]; 
            const label neiCell = nei[meshCell[faceI]]; 

            if ((sign(signedDistance[ownCell]) * sign(signedDistance[neiCell])) < 0)
            {
                signedDistanceSwitches = true; 
                break; 
            }
        }
        
        if (signedDistanceSwitches)
        {
            scalar searchDistance = sqrt(searchDistanceSqr[cellI]);

            markerField[cellI] = 0.5 * (
                1 + signedDistance[cellI] / searchDistance + 1/pi *
                sin((pi * signedDistance[cellI]) / searchDistance)
            );
        }
        else if (signedDistance[cellI] > 0)
        {
            markerField[cellI] = 1;
        }
        else if (signedDistance[cellI] < 0)
        {
            markerField[cellI] = 0;
        }

    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

