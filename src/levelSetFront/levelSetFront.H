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

Class
    Foam::levelSetFront

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany


Description

SourceFiles
    levelSetFront.C

\*---------------------------------------------------------------------------*/

#ifndef levelSetFront_H
#define levelSetFront_H

#include "triSurface.H"
#include "isoSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace frontTracking
{

/*---------------------------------------------------------------------------*\
                         Class levelSetFront Declaration
\*---------------------------------------------------------------------------*/

class levelSetFront
:
    public triSurface 
{

public:

    // Constructors

        //- Construct null
        levelSetFront();

        //- Construct as isoSurface from a distance field
        levelSetFront(const volScalarField& psi, 
                      const scalarField& psiPoint);

        //- Construct as copy
        levelSetFront(const levelSetFront&);


    //- Destructor
    ~levelSetFront();

    // Member Functions
    
        // Edit
        void reconstruct(const volScalarField& cellIsoVals, 
                         const scalarField& pointIsoVals, 
                         const bool regularise=true, 
                         const scalar mergeTol=1E-6);

        // Displace points.
        
    // Member operators 
    void operator=(const isoSurface& rhs);
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frontTracking

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
