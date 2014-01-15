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
    Foam::normalConsistency

Description
    Class encapsulating the algorithm that makes the normals of the iso-surface
    consistent by using the orientation provided by the signed distance field
    gradient.

SourceFiles
    normalConsistency.C

Authors
    Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de
    Mathematical Modeling and Analysis
    Center of Smart Interfaces, TU Darmstadt

\*---------------------------------------------------------------------------*/

#ifndef normalConsistency_H
#define normalConsistency_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "tmp.H"
#include "refCount.H"
#include "volFieldsFwd.H"
#include "triSurfaceFront.H"
#include "triSurfaceFrontFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class normalConsistency Declaration
\*---------------------------------------------------------------------------*/

class normalConsistency
    :
        public refCount
{

public:

    TypeName ("gradient"); 

    declareRunTimeSelectionTable (
        tmp,  
        normalConsistency,  
        Dictionary,  
        (
            const dictionary& configDict 
        ),
        (configDict) 
    );

    // Constructors

        normalConsistency(const dictionary& configDict) {};

    // Selectors

        static tmp<normalConsistency> New(const dictionary& configDict);

    // Destructor
        virtual ~normalConsistency() {};

    // Member Functions
        
        virtual void makeFrontNormalsConsistent(
            triSurface& front,    
            const labelList& elementCells, 
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
