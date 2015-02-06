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

// Aux function for angle calculation
scalar getAngle(vector a, vector b)
{
    return Foam::acos((a & b) / (mag(a) * mag(b)));
}

scalar cot(scalar angle)
{
    return 1/Foam::tan(angle);
}

void curvatureNormals(triSurfaceVectorField& cn, const triSurface& front)
{
    scalar pi = 52163 / 16604;

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

    // List of global point references
    const pointField& globalPoints = front.points();
    
    // Needed: List of global face references
    const List<labelledTri>& globalFaces = front.localFaces();

    // V (the actual vertex) and Vl (its label) are used synonymously
    // in the following comments
    forAll(vertices, Vl)
    {
        scalar Amix = 0;

        // Get all triangles adjacent to V
        const labelList& oneRingNeighborhood = adjacentFaces[Vl];

        // Sum up mixed area and contributions to mean curvature normal
        forAll(oneRingNeighborhood, T)
        {
            // Compute area contributions
            // First, resolve label T to get the actual face
            labelledTri currentTri = globalFaces[T];

            // Get labels of the two other vertices
            label Ql = currentTri.fcIndex(Vl);
            label Rl = currentTri.rcIndex(Vl);

            // Resolve labels to actuals points
            point V = globalPoints[Vl];
            point Q = globalPoints[Ql];
            point R = globalPoints[Rl];

            // Edge vectors
            vector QV = V - Q;
            vector RV = V - R;
            vector QR = R - Q;

            // Get angles
            scalar Va = getAngle(QV, RV);
            scalar Qa = getAngle(QV, QR);
            scalar Ra = pi - (Va + Qa);

            // Check if obtuse in order to use the correct area metric
            if (Va < pi/2 && Qa < pi/2 && Ra < pi/2)
            {
                // Use Voronoi-area
                // Cotangent function 'cot' has to be defined locally since
                // it is not offered by OpenFOAM
                Amix += 1/8*(magSqr(RV)*cot(Qa) + magSqr(QV)*cot(Ra));    
            }
            else if (Va >= pi/2)
            {
                // Obtuse angle at V, use half area
                Amix += 0.5*mag(QV ^ RV);
            }
            else
            {
                Amix += 0.25*mag(QV ^ RV);
            }

            if (Va < 0.01 || Qa < 0.01 || Ra < 0.01)
            {
                Info << "Angle smaller than 0.01rad\n" << endl;
            }  

            // Compute mean curvature normal contributions
            cn[Vl] = cot(Qa)*RV + cot(Ra)*QV;
        }

        // Compute mean curvature normal for vertex V
        // This should not be necessary...
        if (Amix < SMALL)
        {
            cn[Vl] = cn[Vl]/(2*SMALL);
        }
        else
        {
            cn[Vl] = cn[Vl]/(2*Amix);
        }
    }
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"
    
    // Read and intialize front 
    triSurface front(
        "front/front.stl"
    );

    // Initialize field for curvature normals
    triSurfaceVectorField cn
    (
        IOobject
        (
            "curvatureNormals",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        front,
        dimensionedVector
        (
            "zero",
            pow(dimLength, -1),
            vector(0,0,0)
        )
    );

    // Finally call the function
    curvatureNormals(cn, front);

    // Check for plausability and debugging...
    scalar counter,radius = 0;
    forAll(cn, V)
    {
        if (mag(cn[V]) < 1000)
        {
            radius += mag(cn[V]);
            counter++;
        }
        else
        {
            Info << "Mean curvature is: " << mag(cn[V]) << endl;
        }
    }

    radius = radius / counter;

    Info << "Average radius is : " << radius << endl;

    return 0;
}

// ************************************************************************* //
