/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::foamIsoSurfaceFrontReconstructor

Description
    Abstract base class for the heaviside function calculation from a signed 
    distance field.

SourceFiles
    foamIsoSurfaceFrontReconstructor.C

Authors
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

\*---------------------------------------------------------------------------*/

#ifndef foamIsoSurfaceFrontReconstructor_H
#define foamIsoSurfaceFrontReconstructor_H

#include "frontReconstructor.H"
#include "normalConsistency.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class foamIsoSurfaceFrontReconstructor Declaration
\*---------------------------------------------------------------------------*/

class foamIsoSurfaceFrontReconstructor
    :
        public frontReconstructor 
{

    scalar mergeTolerance_; 
    Switch regularize_; 

    tmp<normalConsistency> consistencyAlgTmp_;  

public:

    TypeName ("foamIsoSurface"); 

    // Constructors

        foamIsoSurfaceFrontReconstructor(const dictionary& configDict);

    // Destructor
        virtual ~foamIsoSurfaceFrontReconstructor();

    // Member Functions
        
        labelList reconstructFront(
            triSurfaceFront& front, 
            const volScalarField& signedDistance,
            const pointScalarField& pointSignedDistance
        ) const;  

        virtual void forceConsistentNormalOrientation(
            isoSurface& iso, 
            const volScalarField& signedDistance
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
