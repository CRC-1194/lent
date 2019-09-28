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
    Ellipsoid hyper-surface class used for higher-order VoF initialization.

Authors 
    Tomislav Maric, maric@mma.tu-darmstadt.de
    Johannes Kromer, kromer@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis Group, TU Darmstadt

SourceFiles
    ellipsoidHypersurfaceI.H
    ellipsoidHypersurface.C

\*---------------------------------------------------------------------------*/

namespace Foam { namespace FrontTracking {

scalar ellipsoidHypersurface::value(vector x) const
{
    return (
        ((x[0] * x[0]) / aSqr_) + 
        ((x[1] * x[1]) / bSqr_) + 
        ((x[2] * x[2]) / cSqr_) - 1.0
    );
}

vector ellipsoidHypersurface::gradient(vector x) const
{
    return vector
    ( 
        (2 * x[0]) / aSqr_,  
        (2 * x[1]) / bSqr_, 
        (2 * x[2]) / cSqr_
    );
}

}} // End namespace Foam::FrontTracking 

// ************************************************************************* //
