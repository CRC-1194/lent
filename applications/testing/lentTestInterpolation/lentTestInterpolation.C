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

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Test application for the interpolation methods used by the LENT method.

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "interfaceProperties.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "lentMethod.H"
#include "lentInterpolation.H"
#include "lentGtest.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;
using namespace FrontTracking;
using namespace Test;

TEST_F(lentTests, lentTestInterpolation)
{
    extern int mainArgc;
    extern char** mainArgv;

    int argc = mainArgc;
    char** argv = mainArgv;

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    triSurfaceFront front(
        IOobject(
            "front",
            "front",
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    lentMethod lent(front, mesh);

    lent.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);

    lent.reconstructFront(front, signedDistance, pointSignedDistance);

    front.write(); 

    triSurfaceFrontPointVectorField vertexVectorSource 
    (
        IOobject(
            "vertexVectorSource", 
            runTime.timeName(), 
            front,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ), 
        front,
        dimensionedVector(
            "zero", 
            dimless/dimLength, 
            vector(1,0,0)
        )
    );

    // Set the vertexVectorSource to (vertex.x(), 0,0) for testing. 
    const pointField& frontVertices = front.points(); 
    forAll(frontVertices, vertexI)
    {
        vertexVectorSource[vertexI] = vector(frontVertices[vertexI].x(),0,0); 
    }

    volVectorField cellVectorTarget(
        IOobject(
            "cellVectorTarget", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ), 
        mesh, 
        dimensionedVector(
            "zero", 
            dimless/dimLength, 
            vector(0,0,0)
        )
    );

    lentInterpolation interpolation; 
    interpolation.interpolate(vertexVectorSource, cellVectorTarget);
    cellVectorTarget.write(); 

    Info<< "End\n" << endl;
}

int mainArgc;
char** mainArgv;

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    mainArgc = argc;
    mainArgv = argv;

    return RUN_ALL_TESTS();

    return 0;
}

// ************************************************************************* //
