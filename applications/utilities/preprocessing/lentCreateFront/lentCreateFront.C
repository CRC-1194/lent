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

Authors
    Tobias  Tolle  (tolle@mma.tu-darmstadt.de)

Description
    Create a front by intersecting an analytical surface with the mesh
    and subsequent per-cell triangulaion. It is written to "front/front.stl".

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"

#include "analyticalSurface.H"
#include "lentMethod.H"

#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace FrontTracking;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    // Ensure that the front file exists
    std::ofstream frontFile("constant/front.stl");
    frontFile.close();
    
    triSurfaceFront front(
        IOobject(
            "front",
            "front",
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );

    lentMethod lent(front, mesh);

    tmp<analyticalSurface> analyticalSurfaceTmp(
        analyticalSurface::New(lent.dict().subDict("frontSurface"))
    );

    auto& surface = analyticalSurfaceTmp.ref();

    Info << "Constructing front from analytical surface..." << endl;

    surface.setDistance(signedDistance);
    surface.setDistance(pointSignedDistance);

    lent.reconstructFront(front, signedDistance, pointSignedDistance);

    surface.moveFrontToSurface(front);
    surface.makeNormalOrientationConsistent(front);

    front.triSurface::write("constant/front.stl");

    Info << "Done!" << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
    //----------------
    /*
    IOdictionary lentControlDict(
        IOobject(
           "lentSolution",
           "system",
           mesh.time(),
           IOobject::MUST_READ_IF_MODIFIED,
           IOobject::AUTO_WRITE
        )
    );



    frontConstructor buildFront(analyticalSurfaceTmp, mesh);

    triSurface front = buildFront.createTriSurface();

    // Regularize surface mesh
    frontSmoother smoother{0.33, 3, "pointsAndEdges"};
    smoother.smoothFront(front, mesh);

    analyticalSurfaceTmp.ref().moveFrontToSurface(front);

    front.write("constant/front.stl");

    return 0;
    */
}


// ************************************************************************* //
