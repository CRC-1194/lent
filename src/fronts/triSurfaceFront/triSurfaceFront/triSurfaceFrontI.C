/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Authors
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

\*---------------------------------------------------------------------------*/

#include "triSurfaceFront.H"
#include "volPointInterpolation.H"
#include "pointFields.H"

// * * * * * * * * * * * * * Member Function Templates * * * * * * * * * * * //

template<class Displacement>
void Foam::FrontTracking::triSurfaceFront::executeMovePoints(
    const Displacement& d
)
{
    //FIXME: remove the freaking copy of points, hardcoded to triSurface. 
    pointField newPoints(points()); 
    newPoints += d; 
    movePoints(newPoints);
}

// ************************************************************************* //
