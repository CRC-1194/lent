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
    tomislav.maric@gmx.com
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

#include "isoSurface.H"
#include "triSurfaceMesh.H"

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
    public triSurfaceMesh 
{

    // Front is moving.
    bool moving_; 
    // Front is topologically changing.
    bool changing_;

    // Both file names are stored because IOobject is an ambiguous base
    // class of the triSurfaceMesh, preventing calls to IOobject::instance()
    // and IOobject::name() member functions. TM 04/12/12
    
    // Name.
    fileName name_; 

    // Instance. 
    fileName instance_; 

    // Prepend zeros. 
    label prependZeros_; 

    // Compute the iso-surface.
    void computeIsoSurface (
        const volScalarField& cellsToElementsDist, 
        const scalarField& pointsToElementsDist, 
        const bool regularise=true, 
        const scalar mergeTol=1e-06
    );

    // Write the indexed front file.
    void write(label index);

    // Access to static write control keywords.  
    
    // This data is protected in the Foam::Time, so it is redefined
    // here to allow for output checking during writing. 
    // The writing is separated for the frontMesh from the ObjectRegistry
    // by upstream design and is kept so here. TM 04/12/2012 
    word wcTimeStep()
    {
        static word wcTimeStep_ = "timeStep";
        return wcTimeStep_;
    }
    word wcAdjustableRunTime()
    {
        static word wcAdjustableRunTime_ = "adjustableRunTime";
        return wcAdjustableRunTime_;
    }
    
    // Dissalow: 
    
    //- Construct as copy
    levelSetFront(const levelSetFront&);

public:
       
    TypeName ("levelSetFront");

    // Constructors

    //- Construct read 
    levelSetFront(const IOobject& io, label prependZeros=8);

    //- Construct as isoSurface from a distance field
    //levelSetFront(const volScalarField& psi, 
        //const scalarField& psiPoint);

    //- Destructor
    //~levelSetFront();

    // Member Functions
    
    // Set/get the moving/changing flags. 
    bool isMoving() const;
    void setMoving(bool b); 

    bool isChanging() const; 
    void setChanging(bool b);

    // Edit
    void reconstruct (
        const volScalarField& cellsToElementsDist, 
        const scalarField& pointsToElementsDist, 
        const bool regularise=true, 
        const scalar mergeTol=1e-06
    );

    void reconstruct (
        const volScalarField& cellsToElementsDist, 
        const bool regularise=true, 
        const scalar mergeTol=1e-06
    );

    // Displace points.
    void move(const vectorField& Dv);  

    // Write
    void write(const Time& runTime);
        
    // Member operators 
    void operator=(const isoSurface& rhs);
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frontTracking

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
