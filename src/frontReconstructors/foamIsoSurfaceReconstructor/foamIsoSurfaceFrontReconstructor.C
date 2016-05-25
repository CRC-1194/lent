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
    Foam::foamIsoSurfaceFrontReconstructor

SourceFiles
    diffuseInterfaceProperties.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Abstract base class for the heaviside function calculation from a signed
    distance field.

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


#include "foamIsoSurfaceFrontReconstructor.H"
#include "addToRunTimeSelectionTable.H"
#include "taylorIsoSurface.H" 
//#include "isoSurfaceCell.H"// Alternative iso-surface reconstruction.
#include "lentCommunication.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(foamIsoSurfaceFrontReconstructor, 0);
    addToRunTimeSelectionTable(foamIsoSurfaceFrontReconstructor, foamIsoSurfaceFrontReconstructor, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

foamIsoSurfaceFrontReconstructor::foamIsoSurfaceFrontReconstructor(
   const dictionary& configDict
)
:
    frontReconstructor(configDict),
    mergeTolerance_(readScalar(configDict.lookup("mergeTolerance"))),
    regularize_(configDict.lookup("regularization")),
    consistencyAlgTmp_(normalConsistency::New(configDict.subDict("normalConsistency")))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void foamIsoSurfaceFrontReconstructor::reconstructFront(
    triSurfaceFront& front,
    const volScalarField& signedDistance,
    const pointScalarField& pointSignedDistance
) const
{
    taylorIsoSurface iso(
        signedDistance,
        pointSignedDistance,
        0, // FIXME: Leave as a compile time constant? TM.
        bool(regularize_), 
        mergeTolerance_
    );

    // Clean up degenerate triangles. Report the cleanup process.  
    // FIXME: Investigate this.
    //iso.cleanup(true);

    // Update the communication map after reconstruction.
    const auto& mesh = signedDistance.mesh();
    auto& communication = const_cast<lentCommunication&>(
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
        ) 
    );
    // Assign the new iso surface mesh to the front.  
    front = static_cast<triSurface&>(iso);

    // Set the triangleToCell using the map produced by the iso-surface
    // reconstruction algorithm.
    communication.setTriangleToCell(iso.meshCells());
    // Update vertex to cell using triangle to cell information.
    // - NOTE: Nearest information gets completely lost: a completely new
    // set of triangles is generated!  
    communication.updateVertexToCell();

    // Make normals consistent. 
    // FIXME: Remove consistencyAlgTmp_, not required anymore, 
    // normals consistently oriented. TM.
    consistencyAlgTmp_->makeFrontNormalsConsistent(
        front,
        signedDistance, 
        pointSignedDistance
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
