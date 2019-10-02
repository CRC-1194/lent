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

scalar ellipsoidHypersurface::value(point p) const
{
    double result {-1.};  

    for (char i = 0; i < 3; ++i)
        result += Foam::sqr((p[i] - center_[i]) / halfAxes_[i]); 

    return result;
}

vector ellipsoidHypersurface::gradient(point p) const
{
    return vector
    ( 
        2*(p[0] - center_[0]) / Foam::sqr(halfAxes_[0]),
        2*(p[1] - center_[1]) / Foam::sqr(halfAxes_[1]),
        2*(p[2] - center_[2]) / Foam::sqr(halfAxes_[2])
    );
}

}} // End namespace Foam::FrontTracking 

// ************************************************************************* //
