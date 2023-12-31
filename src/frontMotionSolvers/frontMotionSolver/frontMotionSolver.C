/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::frontMotionSolver

SourceFiles
    frontMotionSolver.C

Authors:
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    Interface for the Front Tracking advection. 

Affiliations:
    Mathematical Modeling and Analysis Institute, Mathematics Department, 
    TU Darmstadt, Germany

Funding:
    German Research Foundation (DFG) - Project-ID 265191195 - SFB 1194

    German Research Foundation (DFG) - Project-ID MA 8465/1-1, 
    Initiation of International Collaboration 
    "Hybrid Level Set / Front Tracking methods for simulating 
    multiphase flows in geometrically complex systems"

\*---------------------------------------------------------------------------*/


#include "frontMotionSolver.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontMotionSolver, 0);
    defineRunTimeSelectionTable(frontMotionSolver, Dictionary)

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontMotionSolver::frontMotionSolver(const dictionary& configDict)
    :
        cellDisplacementTmp_(nullptr),
        frontDisplacementTmp_(nullptr),
        interpolation_(configDict)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

tmp<frontMotionSolver>
frontMotionSolver::New(const dictionary& configDict)
{
    const word name = configDict.get<word>("type");

    // Find the constructor pointer for the model in the constructor table.
    auto* ctorPtr = DictionaryConstructorTable(name);

    // If the constructor pointer is not found in the table.
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            configDict,
            "frontMotionSolver",
            name,
            *DictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    // Construct the model and return the autoPtr to the object.
    return tmp<frontMotionSolver> (ctorPtr(configDict));
}


// * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

void frontMotionSolver::initDisplacements(
    const triSurfaceFront& front,
    const volVectorField& cellVelocity
)
{
    const auto& mesh = cellVelocity.mesh(); 
    const auto& runTime = mesh.time(); 

    dimensionedVector zeroDisplacement(
        "zero",
        dimLength, 
        vector(0,0,0)
    );

    if (!cellDisplacementTmp_)
    {
        cellDisplacementTmp_ = tmp<volVectorField>(
            new volVectorField( 
                IOobject(
                    "cellDisplacement",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                zeroDisplacement
            )
        );
    } 

    if (!frontDisplacementTmp_)
    {
        frontDisplacementTmp_ = tmp<triSurfaceFrontPointVectorField>(
            new triSurfaceFrontPointVectorField( 
                IOobject(
                    "frontDisplacement",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                front,
                zeroDisplacement
            )
        );
    }
    else if (frontDisplacementTmp_->size() != front.nPoints())
    {
        frontDisplacementTmp_->resize(front.nPoints()); 
    }
}

void frontMotionSolver::evolveFront(
    triSurfaceFront& front,
    const volVectorField& cellVelocity 
)
{
    initDisplacements(front,cellVelocity); 

    // Overload this member function for different functionality.
    calcCellDisplacement(cellVelocity); 

    // Interpolate the displacement field from the cell to the front.
    const auto& deltaC = cellDisplacements();
    auto& deltaF = frontDisplacements();
    interpolation_.interpolate(deltaC, deltaF); 

    // Displace front points with front displacements.  
    front.displace(deltaF);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
