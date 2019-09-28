/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 AUTHOR,AFFILIATION
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

Class
    Foam::hypersurfaceData
    Foam::weingartenData

Description

    Data for the hypersurface information used by the higher-order initialization.

Authors
    Tomislav Maric, maric@mma.tu-darmstadt.de
    Johannes Kromer, kromer@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis Group, TU Darmstadt

SourceFiles
    hypersurfaceData.C
    hypersurfaceDataI.H

\*---------------------------------------------------------------------------*/

#include "hypersurfaceData.H"

namespace Foam 
{
namespace FrontTracking
{

/*---------------------------------------------------------------------------*\
                         Class weingartenData Definition 
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

weingartenData::weingartenData()
:
    tangent1_(VGREAT,VGREAT,VGREAT),
    tangent2_(VGREAT,VGREAT,VGREAT),
    kappa1_(VGREAT),
    kappa2_(VGREAT),
    normal_(VGREAT,VGREAT,VGREAT) 
{}

/*---------------------------------------------------------------------------*\
                         Class hypersurfaceData Definition 
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hypersurfaceData::hypersurfaceData()
:
    basePoint_(VGREAT, VGREAT, VGREAT), 
    weingartenData_()
{}

// ************************************************************************* //

} // End namespace FrontTracking 
} // End namespace Foam
