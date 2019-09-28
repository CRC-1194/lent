/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Tomislav Maric, Johannes Kromer
     \\/     M anipulation  |                    TU Darmstadt
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
    Foam::ellipsoidHypersurface

Description
    Ellipsoid hyper-surface class used for higher-order initialization.

Authors 
    Tomislav Maric, maric@mma.tu-darmstadt.de
    Johannes Kromer, kromer@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis Group, TU Darmstadt

SourceFiles
    ellipsoidHypersurfaceI.H
    ellipsoidHypersurface.C

\*---------------------------------------------------------------------------*/

#include "ellipsoidHypersurface.H"

namespace Foam { namespace FrontTracking {

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ellipsoidHypersurface::ellipsoidHypersurface(scalar a, scalar b, scalar c)
:
    a_(a), 
    b_(b), 
    c_(c), 
    aSqr_(a * a), 
    bSqr_(b * b), 
    cSqr_(c * c)
{}

}} // End namespace Foam::FrontTracking 

// ************************************************************************* //
