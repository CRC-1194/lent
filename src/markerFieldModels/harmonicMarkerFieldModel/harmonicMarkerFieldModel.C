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
    Foam::harmonicMarkerFieldModel

SourceFiles
    diffuseInterfaceProperties.C

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


#include "harmonicMarkerFieldModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(harmonicMarkerFieldModel, 0);
    addToRunTimeSelectionTable(markerFieldModel, harmonicMarkerFieldModel, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

harmonicMarkerFieldModel::harmonicMarkerFieldModel(const dictionary& configDict)
:
    markerFieldModel(configDict), 
    cellDistFieldName_(configDict.lookup("cellDistanceField")),
    pointDistFieldName_(configDict.lookup("pointDistanceField"))
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

harmonicMarkerFieldModel::~harmonicMarkerFieldModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void harmonicMarkerFieldModel::calcMarkerField(const fvMesh& mesh) const
{
    volScalarField& markerField = 
        mesh.lookupObject<volScalarField>(markerFieldName()); 

    const volScalarField& signedDistance = 
        mesh.lookupObject<volScalarField>(cellDistFieldName_); 

    const volScalarField& searchDistanceSqr = 
        mesh.lookupObject<volScalarField>(markerFieldName()); 

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

