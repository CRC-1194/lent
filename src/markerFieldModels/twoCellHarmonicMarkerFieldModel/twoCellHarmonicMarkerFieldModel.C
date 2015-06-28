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
    Foam::twoCellHarmonicMarkerFieldModel

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Abstract base class for the markerField function calculation from a signed
    distance field.

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

