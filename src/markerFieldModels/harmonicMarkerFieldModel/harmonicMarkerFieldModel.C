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
    Foam::harmonicMarkerFieldModel

SourceFiles
    harmonicMarkerFieldModel.C

Authors
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


#include "harmonicMarkerFieldModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(harmonicMarkerFieldModel, 0);
    addToRunTimeSelectionTable(markerFieldModel, harmonicMarkerFieldModel, Dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void harmonicMarkerFieldModel::calcMarkerField(volScalarField& markerField) const
{
    const fvMesh& mesh = markerField.mesh(); 

    const volScalarField& signedDistance = 
        mesh.lookupObject<volScalarField>(cellDistFieldName()); 

    const volScalarField& searchDistanceSqr = 
        mesh.lookupObject<volScalarField>(sqrSearchDistFieldName()); 

    scalar pi = constant::mathematical::pi;

    forAll (markerField, cellI)
    {
        scalar searchDistance = sqrt(searchDistanceSqr[cellI]);

        if (mag(signedDistance[cellI]) < searchDistance)
        {
            markerField[cellI] = 0.5 * (
                1 + signedDistance[cellI] / searchDistance + 1/pi *
                sin((pi * signedDistance[cellI]) / searchDistance)
            );
        }
        else
        {
            if (signedDistance[cellI] > 0)
            {
                markerField[cellI] = 1;
            }
            if (signedDistance[cellI] < 0)
            {
                markerField[cellI] = 0;
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

