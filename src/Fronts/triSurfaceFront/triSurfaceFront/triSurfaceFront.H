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
    Foam::triSurfaceFront

Author
    Tomislav Maric
    tomislav.maric@gmx.com

Description

SourceFiles
    triSurfaceFront.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceFront_H
#define triSurfaceFront_H

#include "isoSurface.H"
#include "triSurfaceMesh.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FrontTracking
{

/*---------------------------------------------------------------------------*\
                         Class triSurfaceFront Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceFront
    :
        public regIOobject, 
        public triSurface 
{
    // Write format. 
    word writeFormat_;

    // Prepend zeros. 
    label prependZeros_; 

    // Compute the iso-surface.
    void computeIsoSurface(
        const volScalarField& cellsToElementsDist, 
        const scalarField& pointsToElementsDist, 
        const bool regularise=true, 
        const scalar mergeTol=1e-06
    );

    // Write the indexed front file.
    fileName zeroPaddedFileName(word extension) const;

    template<class Displacement>
    void executeMovePoints(const Displacement& d); 

    void forceConsistentNormalOrientation 
    (
        const isoSurface& reconstructedFront, 
        const volScalarField& cellsToElementsDist
    ); 

public:
       
    TypeName ("triSurfaceFront");

    // Constructors

    //- Construct read 
    triSurfaceFront(const IOobject& io, word writeFormat="vtk", label prependZeros=8);

    //- Destructor
    virtual ~triSurfaceFront() {};

    // Member Functions
    //- Edit
    void reconstruct(
        const volScalarField& cellsToElementsDist, 
        const scalarField& pointsToElementsDist, 
        const bool regularise=true, 
        const scalar mergeTol=1e-06
    );

    void reconstruct(
        const volScalarField& cellsToElementsDist, 
        const bool regularise=true, 
        const scalar mergeTol=1e-06
    );

    // Move points with a spatially constant displacement vector.
    void move(vector deltaV);

    // Move points with a displacement vector field.
    void move(const vectorField& deltaV);  

    //- Write
    
    bool write() const;  

    bool writeObject
    (
        IOstream::streamFormat fmt,
        IOstream::versionNumber ver,
        IOstream::compressionType cmp
    ) const;

    bool writeData(Foam::Ostream& os) const; 
        
    // Member operators 
    void operator=(const triSurfaceFront& rhs);
    void operator=(const isoSurface& rhs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# include "triSurfaceFrontI.C"

#endif

// ************************************************************************* //
