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

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    tomislav.maric@gmx.com
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#include "levelSetFront.H"
#include "volPointInterpolation.H"
#include "pointFields.H"

// * * * * * * * * * * * * * Member Function Templates * * * * * * * * * * * //

template <class T>
T Foam::levelSetFrontTracking::levelSetFront::readWriteInterval (
    const dictionary& controlDict
)
{
    T writeInterval = controlDict.lookupOrDefault<T> (
        "writeInterval",
        -1
    );

    if (writeInterval == -1)
    {
        FatalErrorIn (
            "levelSetFront::write(const Time& runTime)"
        ) << "Wrong writeInterval value in the controlDict dictionary."
            << abort(FatalError);
    }

    return writeInterval;
}

template<class Displacement>
void Foam::levelSetFrontTracking::levelSetFront::executeMovePoints (
    const Displacement& d
)
{
    pointField newPoints(points()); 

    newPoints += d; 

    movePoints(newPoints);

    setMoving(true);

}

// ************************************************************************* //
