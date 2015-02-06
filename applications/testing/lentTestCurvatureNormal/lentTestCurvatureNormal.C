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
    Test application for the interface reconstruction algorithm of the LENT method.

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

#include "fvCFD.H"
#include "triSurfaceFields.H"

#include "interfaceProperties.H"
#include "incompressibleTwoPhaseMixture.H"
#include "lentMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//using namespace FrontTracking;
//using namespace Test;

void curvatureNormals(triSurfaceVectorField& cn, const triSurface& front)
{
    cn = dimensionedVector("zero",
                            dimless/dimLength,
                            vector(0,0,0)
                          );

    // Get label list of all mesh vertices
    // Taken from file "PrimitivePatch.H", ll.347
    const labelList& vertices = front.meshPoints();

    // Following line should yield a list of labelLists, which
    // contains the mapping between a vertex and a face (triangle)
    const labelListList& adjacentFaces = front.pointFaces();

    forAll(vertices, V)
    {
        scalar Amix = 0;

        // Get all triangles adjacent to V
        const labelList& oneRingNeighborhood = adjacentFaces[V];

        // Sum up mixed area and contributions to mean curvature normal
        forAll(oneRingNeighborhood, T)
        {
            // Get other two vertices of T
            // --> Kann so nicht funktionieren, da T vom Typ 'label' ist.
            // Müsste wenn dann so aussehen:
            // label bla = faces[T].nextLabel(V);
            /*
            label R = T.nextLabel(V);
            label Q = T.prevLabel(V);
            */

        }


    }

}

/*
TEST_F(lentTests, lentReconstruction)
{
    extern int mainArgc;
    extern char** mainArgv;

    int argc = mainArgc;
    char** argv = mainArgv;


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    lentMethod lent(front, mesh);

    lent.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);

    lent.reconstructFront(front, signedDistance, pointSignedDistance);

    front.write();

    while (runTime.run()) {

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
}
*/

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"
    
    // Read and intialize front 
    triSurface front(
        /*
        IOobject(
            "front.stl",
            "front",
            runTime,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
        */
        "front/front.stl"
    );

    // Initialize field for curvature normals
    triSurfaceVectorField curvatureNormals
    (
        IOobject
        (
            "curvatureNormals",
            "curvatureNormals",
            front,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        front,
        dimensionedScalar
        (
            "zero,"
            pow(dimLength, -1),
            0
        );
    );


    return 0;
}

// ************************************************************************* //
