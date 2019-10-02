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

#ifndef ellipsoidHypersurface_H
#define ellipsoidHypersurface_H

#include "hypersurfaceData.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam { namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class ellipsoidHypersurface Declaration
\*---------------------------------------------------------------------------*/

class ellipsoidHypersurface
{
    vector center_; 
    vector halfAxes_;
    word name_;

public:

    // Constructors

        //- Construct from components
        ellipsoidHypersurface(
            const vector& center, 
            const vector& halfAxes, 
            const word name = "ellipsoid"
        );

    // Member Functions
    
        word name() const
        {
            return name_; 
        }
        
        inline scalar value(point p) const;

        inline vector gradient(point p) const;

        inline vector center() const
        {
            return center_; 
        }
        
        inline vector halfAxes() const
        {
            return halfAxes_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}} // End namespace Foam::FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ellipsoidHypersurfaceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
