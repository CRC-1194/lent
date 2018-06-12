/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 AUTHOR,AFFILIATION
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
    Foam::lentSolutionControl

Description
    Class to control the solution procedure of the discretized
    Navier Stokes system.
    Allows for more fine grained control than the pimple control class
    it inherits from.

    Required parameters to be set in fvSolution in addition to PIMPLE
    parameters:
        - updateMomentumFlux (yes/no)
        - adaptiveInnerIteration (yes/no)
        - velocityUpdateType (none/adaptive/alwaysOn)
    
    if velocityUpdateType == "adaptive":
        - phiChangeTolerance scalar
        
SourceFiles
    lentSolutionControl.C

\*---------------------------------------------------------------------------*/

#ifndef lentSolutionControl_H
#define lentSolutionControl_H

#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class lentSolutionControl Declaration
\*---------------------------------------------------------------------------*/

class lentSolutionControl
:
    public pimpleControl
{

    // Private data
    
    // Field references
    const surfaceScalarField& volumetricFluxes_;

    // User defined flags for solution control
    bool resolveConvectiveNonlinearity_;
    bool adaptiveInnerIteration_;
    // Valid parameters: "none", "adaptive", "alwaysOn"
    word explicitVelocityUpdateType_;
    scalar relTolVolFlux_;

    // Internal parameters for solution control
    bool explicitVelocityUpdate_;

    // Private Member Functions
    bool volFluxesAreConverged();
    bool pressureIsConverged();
    void checkParameterCompatibility() const;
    void displayConfiguration() const;

    
protected:

    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual void read();


public:

    // Static data members
    
        //- Run-time type information
        TypeName("lentSolutionControl");


    // Constructors
    lentSolutionControl(fvMesh& mesh, const surfaceScalarField& phi, const word& dictName="lentSC");


    //- Destructor
    virtual ~lentSolutionControl() = default;


    // Member Functions

        // Access

        // Solution control
        bool adaptiveCorrect();
        bool explicitVelocityUpdate();
        bool updateMomentumFlux() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
