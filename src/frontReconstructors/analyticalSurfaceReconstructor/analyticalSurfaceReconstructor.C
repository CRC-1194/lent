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
    Foam::analyticalFrontReconstructor

SourceFiles
    analyticalFrontReconstructor.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Reconstruct a front from an analyticalSurface object using the
    frontConstructor class.

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

#include "addToRunTimeSelectionTable.H"

#include "frontConstructor.H"
#include "lentCommunication.H"

#include "analyticalSurfaceReconstructor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalSurfaceReconstructor, 0);
    addToRunTimeSelectionTable(frontReconstructor, analyticalSurfaceReconstructor, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalSurfaceReconstructor::analyticalSurfaceReconstructor(const dictionary& configDict)
:
    frontReconstructor(configDict),
    analyticalSurfaceTmp_(
        analyticalSurface::New(configDict.subDict("analyticalSurface"))
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void analyticalSurfaceReconstructor::reconstructFront(
    triSurfaceFront& front,
    const volScalarField& signedDistance,
    const pointScalarField&
) const
{
    // Update the communication map after reconstruction.
    const auto& mesh = signedDistance.mesh();
    auto& communication = const_cast<lentCommunication&>(
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
        ) 
    );
    frontConstructor frontCon(analyticalSurfaceTmp_, mesh);

    front = frontCon.createTriSurface();

    // Set the triangleToCell using the map produced by the iso-surface
    // reconstruction algorithm.
    communication.setTriangleToCell(frontCon.triangleToCell());
    communication.updateVertexToCell();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
