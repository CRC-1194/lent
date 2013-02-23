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
    Foam::frontTrackingCalculator

Description

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    tomislav.maric@gmx.com
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany


SourceFiles
    frontTrackingCalculator.C

\*---------------------------------------------------------------------------*/

#ifndef frontTrackingCalculator_H
#define frontTrackingCalculator_H

#include "fvMesh.H"
#include "frontToMeshDistanceCalculator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace frontTracking 
{

/*---------------------------------------------------------------------------*\
                         Class frontTrackingCalculator Declaration
\*---------------------------------------------------------------------------*/

template<class MeshAndFrontConnection>
class frontTrackingCalculator
{
    // Private data
    
    // RTS selected calculators
    autoPtr<meshAndFrontDistanceCalculator<MeshAndFrontConnection> > distanceCalculator_;

public:

    // Static data members

    // Constructors

        //- Construct from components
        frontTrackingCalculator (); //const fvMesh& mesh, const levelSetFront& front);

        //- Construct as copy
        //frontTrackingCalculator(const frontTrackingCalculator&);


    //- Destructor
    //~frontTrackingCalculator();


    // Member Functions

        // Calculate 
        void calcCentresToElementsDistance (volScalarField& Psi, 
            const MeshAndFrontConnection& connectivity); 

        void calcPointsToElementsDistance (scalarField& psi, 
            const MeshAndFrontConnection& connectivity); 

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        // TODO: define
        //void operator=(const frontTrackingCalculator&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frontTracking 
} // End namespace Foam

#ifdef NoRepository
#   include "frontTrackingCalculator.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
