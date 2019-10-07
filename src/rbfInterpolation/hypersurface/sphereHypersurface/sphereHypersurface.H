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

#ifndef sphereHypersurface_H
#define sphereHypersurface_H

#include "hypersurfaceData.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam { namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class sphereHypersurface Declaration
\*---------------------------------------------------------------------------*/

class sphereHypersurface
{

    vector center_; 
    scalar radius_; 
    word name_;

public:

    // Constructors

        //- Construct from components
        sphereHypersurface(
            const vector& center, 
            const scalar& radius, 
            const word& name = "sphere"
        );

        sphereHypersurface(
            const scalar& radius, 
            const vector& center=vector(0.,0.,0.),
            const word& name = "sphere"
        );

    // Member Functions

        inline word name() const
        {
            return name_;
        }

        // Access
        
        inline scalar value(vector x) const;

        inline scalar isoValue() const; 

        inline vector gradient(vector x) const;
};

}} // End namespace Foam::FrontTracking 

#include "sphereHypersurfaceI.H"

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
