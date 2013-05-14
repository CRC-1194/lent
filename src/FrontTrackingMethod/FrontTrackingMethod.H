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
    Foam::FrontTrackingMethod

Description
    Front tracking method policy host. 
    
Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    tomislav.maric@gmx.com
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

SourceFiles
    FrontTrackingMethod.C

\*---------------------------------------------------------------------------*/

#ifndef FrontTrackingMethod_H
#define FrontTrackingMethod_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace frontTracking 
{

/*---------------------------------------------------------------------------*\
                         Class FrontTrackingMethod Declaration
\*---------------------------------------------------------------------------*/

template<class DistanceMethod>
class FrontTrackingMethod
{
    // Private data
    
    // Calculates the distance between the front and the mesh
    //DistanceMethod distanceMethod_; 
    
public:

    // Member Functions

        // Calculate 
        template <class Connection>
        void calcCentresToElementsDistance(volScalarField& Psi, 
            const Connection& connection) const; 

        template <class Connection>
        void calcPointsToElementsDistance(scalarField& psi, 
            const Connection& connection) const; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frontTracking 
} // End namespace Foam

//#ifdef NoRepository
//#   include "FrontTrackingMethod.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
