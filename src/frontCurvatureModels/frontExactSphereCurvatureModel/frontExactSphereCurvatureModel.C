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
    Foam::frontExactSphereCurvatureModel

SourceFiles
    frontExactSphereCurvatureModel.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#include "frontExactSphereCurvatureModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontExactSphereCurvatureModel, 0);
    addToRunTimeSelectionTable(frontCurvatureModel, frontExactSphereCurvatureModel, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontExactSphereCurvatureModel::frontExactSphereCurvatureModel(const dictionary& configDict)
    :
        frontExactCircleCurvatureModel(configDict)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar frontExactSphereCurvatureModel::curvatureAtPoint(const point& p) const
{
    scalar curvature = 0;

    scalar curvRadius = curvatureRadius(p); 

    if (curvRadius > SMALL) 
    {
        curvature = 2 / curvRadius; 
    }
    else
    {
        curvature = 0; 
    }

    return curvature;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
