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
    Foam::maxNormalAngleFrontReconstructionModel

SourceFiles
     maxNormalAngleFrontReconstructionModel.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Abstract base class for the heaviside function calculation from a signed
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


#include "maxNormalAngleFrontReconstructionModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(maxNormalAngleFrontReconstructionModel, 0);
    addToRunTimeSelectionTable(frontReconstructionModel, maxNormalAngleFrontReconstructionModel, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

maxNormalAngleFrontReconstructionModel::maxNormalAngleFrontReconstructionModel(const dictionary& configDict)
:
    frontReconstructionModel(configDict),
    maxAngle_(readScalar(configDict.lookup("maxAngle")) * M_PI / 180.0),
    minAngleCos_(Foam::cos(maxAngle_))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool maxNormalAngleFrontReconstructionModel::reconstructionRequired(
    const triSurfaceFront& front,
    const volScalarField& signedDistance
) const
{
    if (signedDistance.time().timeIndex() <= 1)
        return true; 

    const auto& edges = front.edges(); 
    const auto& allEdgeFaces = front.edgeFaces(); 
    const auto& faceNormals = front.faceNormals(); 

    forAll(edges, I)
    {
        const auto& edgeFaces = allEdgeFaces[I];
        const auto& n0 = faceNormals[edgeFaces[0]]; 

        for(label J = 1; J < edgeFaces.size(); ++J)
        {
            const auto& n = faceNormals[edgeFaces[J]]; 

            if ((n0 & n) < minAngleCos_)
                return true; 
        }
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
