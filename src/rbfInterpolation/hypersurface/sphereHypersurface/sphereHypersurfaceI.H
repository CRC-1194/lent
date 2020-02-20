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
    Foam::sphereHypersurface

Description
    Sphere hyper-surface class.

Authors
    Tomislav Maric, maric@mma.tu-darmstadt.de
    Johannes Kromer, kromer@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis Group, TU Darmstadt

SourceFiles
    sphereHypersurfaceI.H
    sphereHypersurface.C

\*---------------------------------------------------------------------------*/

namespace Foam { namespace FrontTracking {

scalar sphereHypersurface::value(vector x) const
{
    return ((x - center_) & (x - center_)) - (radius_ * radius_); 
}

vector sphereHypersurface::gradient(vector x) const
{
    return vector
    (
        2 * (x[0] - center_[0]), 
        2 * (x[1] - center_[1]),
        2 * (x[2] - center_[2])
    ); 
}

}} // End namespace Foam::FrontTracking 

// ************************************************************************* //
