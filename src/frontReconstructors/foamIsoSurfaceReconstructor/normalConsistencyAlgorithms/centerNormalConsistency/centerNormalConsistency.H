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
    Foam::centerNormalConsistency

Description
    Class encapsulating the algorithm that makes the normals of the iso-surface
    consistent by using the signed distance values stored in cell centers. 

    Using a barycentric tetrahedral decomposition results with a monotone signed
    distance field within the cell.

    (Not proved, its an intuitive conclusion: isosurface cannot cross one edge
     twice and all edges span from the geometrical cell center. TM) 

SourceFiles
    centerNormalConsistency.C

Authors
    Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de
    Mathematical Modeling and Analysis
    Center of Smart Interfaces, TU Darmstadt

\*---------------------------------------------------------------------------*/

#ifndef centerNormalConsistency_H
#define centerNormalConsistency_H

#include "normalConsistency.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class centerNormalConsistency Declaration
\*---------------------------------------------------------------------------*/

class centerNormalConsistency
    :
        public normalConsistency 
{

public:

    TypeName ("centerNormal"); 

    // Constructors
        centerNormalConsistency(const dictionary& configDict); 

    // Selectors

        static tmp<centerNormalConsistency> New(const dictionary& configDict);

    // Destructor
        virtual ~centerNormalConsistency() {};

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
