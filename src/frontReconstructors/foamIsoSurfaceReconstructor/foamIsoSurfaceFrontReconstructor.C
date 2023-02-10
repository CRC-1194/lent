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
    Foam::foamIsoSurfaceFrontReconstructor

SourceFiles
    foamIsoSurfaceFrontReconstructor.C

Author
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    Reconstruct a front from as an iso-surface. 

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


#include "foamIsoSurfaceFrontReconstructor.H"
#include "addToRunTimeSelectionTable.H"
#include "lentCommunication.H"
#include "isoSurfaceTopo.H"

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
    mergeTolerance_(configDict.get<scalar>("mergeTolerance")),
    regularize_(configDict.get<Switch>("regularization")),
    consistencyAlgTmp_(normalConsistency::New(configDict.subDict("normalConsistency"))),
    EV0_(0),
    dict_(configDict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
scalar foamIsoSurfaceFrontReconstructor::surfaceEnclosedVolume(triSurfaceFront& front) const
{
    scalar Vsurf = 0;

    const auto& p = front.localPoints();
    const auto& faces = front.localFaces();
    forAll(faces, I)
    {
        const auto& Sf = front.Sf()[I];
        const auto& f = faces[I];
        Vsurf += dot(-Sf, // Surface normals are oriented into the phase.
            (p[f[0]] + p[f[1]] +
                p[f[2]]));
    }

    Vsurf = 1. / 9. * mag(Vsurf);

    return Vsurf;
}

void foamIsoSurfaceFrontReconstructor::reconstructFront(
    triSurfaceFront& front,
    const volScalarField& signedDistance,
    const pointScalarField& pointSignedDistance
) const
{
    using algorithmType = isoSurfaceParams::algorithmType; 
    using filterType = isoSurfaceParams::filterType;
    //scalar EV = 0;
    scalar isoValue = 0;

    // Original choice for the iso-surface reconstruction 
    isoSurfaceParams isoParams
    (
        dict_,
        algorithmType::ALGO_DEFAULT, 
        (regularize_) ? filterType::DIAGCELL : filterType::NONE
    );

    if (signedDistance.time().timeIndex()>0)
    {
        isoSurfacePoint iso0(
            signedDistance,
            pointSignedDistance,
            0, 
            isoParams
        );
        triSurfaceFront front0 
        (
            IOobject(
                "front",
                "front0",
                signedDistance.mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )//,
          //  front
        );
        // front0= front;
        front0 = iso0;
        //front0.write();
        auto EV = surfaceEnclosedVolume(front0);
        const auto& triArea = front.magSf();
        auto sumSf = gSum(triArea);
        isoValue = (-EV+EV0_)/(sumSf + SMALL);
        Info<<"### the isoValue for volume conservation:" << isoValue <<endl;
        Info<<"### the EV0_, EV, sumSf and SMALL: " << EV0_ << ", " << EV << ", " << sumSf << ", " << SMALL << ". "<<endl;
    }


    isoSurfacePoint iso(
        signedDistance,
        pointSignedDistance,
        //0,
        isoValue, 
        isoParams
    );
    
    //isoSurfaceParams isoParams
    //(
    //    algorithmType::ALGO_DEFAULT,
    //    (regularize_) ? filterType::CELL : filterType::NONE
    //);
    //isoSurfaceTopo iso(
    //    signedDistance.mesh(), 
    //    signedDistance, 
    //    pointSignedDistance, 
    //    0, 
    //    isoParams
    //);

    // Update the communication map after reconstruction.
    const auto& mesh = signedDistance.mesh();
    auto& communication = const_cast<lentCommunication&>(
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
        ) 
    );

    // Assign the new iso surface mesh to the front.  
    front = iso;

    if(signedDistance.time().timeIndex()==0)
    {
        EV0_=surfaceEnclosedVolume(front);// volume at the zero time step
    }

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
    // consistencyAlgTmp_->makeFrontNormalsConsistent(
    //     front,
    //     signedDistance, 
    //     pointSignedDistance
    // );

    communication.update();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
