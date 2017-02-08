/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2.2.x                               
    \\  /    A nd           | Copyright held by original author
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

Author
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    Reconstruct a front from an analyticalSurface object using the
    frontConstructor class.
    
    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

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
    const pointScalarField& pointSignedDistance
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