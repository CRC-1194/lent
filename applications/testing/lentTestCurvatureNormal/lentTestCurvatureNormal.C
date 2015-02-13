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

// Aux functions

// Calculate angle between two vectors
scalar getAngle(vector& a, vector& b)
{
    return Foam::acos((a & b) / (mag(a) * mag(b)));
}

// Cotangent function to resemble algorihm notion from paper
scalar cot(scalar angle)
{
    return 1.0/Foam::tan(angle);
}

// Triangle area functions
// Use heron's formula
scalar areaHeron(vector& A, vector& B, vector& C)
{
    scalar a = mag(A);
    scalar b = mag(B);
    scalar c = mag(C);
    scalar s = 0.5 * (a + b + c);
                
    return Foam::sqrt(s*(s-a)*(s-b)*(s-c));
}

// Use cross product property
scalar areaCross(vector& a, vector& b)
{
    return 0.5*mag(a ^ b);
}

// Quality metric in order to get to the core of the problem
// of completely wrong curvatures
bool badQuality(point& A, point& B, point& C)
{
    bool isBad = false;

    scalar pi = 3.141592653589793238; 

    // Edges
    vector a = C - B;
    vector b = C - A;
    vector c = B - A;

    // Angles
    scalar alpha = getAngle(b, c);
    scalar gamma = getAngle(a, b);
    scalar beta = pi - (alpha + gamma);

    // Area
    scalar area = areaHeron(a, b, c);

    // Actual quality check, hardcoded for now
    scalar minAngle = 0.17453;  // 10 degree
    scalar maxAngle = 2.4435;   // 140 degree

    scalar maxRatio = 5.0; // maximum ratio of edge lengths

    scalar minArea = 1e-5; // Should be smaller than the value computed by checkSTL
    scalar maxArea = 1e-1; // Should be bigger than value computed by checkSTL

    // Angle check
    if (alpha < minAngle || beta < minAngle || gamma < minAngle)
    {
        isBad = true;
    }

    if (alpha > maxAngle || beta > maxAngle || gamma > maxAngle)
    {
        isBad = true;
    }

    // Edge ratio check
    scalar minLength = mag(a);
    scalar maxLength = mag(b);

    minLength = mag(c) < minLength ? mag(c) : minLength;
    minLength = mag(b) < minLength ? mag(b) : minLength;
    maxLength = mag(c) > maxLength ? mag(c) : maxLength;
    maxLength = mag(a) > maxLength ? mag(a) : maxLength;

    isBad = maxRatio > maxLength/minLength ? false : true;

    // Area check
    if (area < minArea || area > maxArea) {isBad = true;}

    return isBad;
}

void curvatureNormals(triSurfaceVectorField& cn, const triSurface& front)
{
    scalar pi = 3.141592653589793238; 

    cn = dimensionedVector("zero",
                            dimless/dimLength,
                            vector(0,0,0)
                          );

    // TODO: modify so that mutliple patches, e.g. in the case of
    // bubble breakup, are supported

    // All of the following functions are taken from
    // "PrimitivePatch.H"
    // Get label list of all mesh vertices
    const labelList& vertices = front.meshPoints();

    // Get assignment point --> faces
    const labelListList& adjacentFaces = front.pointFaces();

    // List of local point references for current patch
    const pointField& localPoints = front.localPoints();
    
    // Needed: List of local face references
    const List<labelledTri>& localFaces = front.localFaces();

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
            label triLabel = oneRingNeighborhood[T];
            labelledTri currentTri = localFaces[triLabel];

            // Read point labels from triangle and determine which
            // one matches Vl to resolve point coordinates correctly
            label tri0 = currentTri[0];
            label tri1 = currentTri[1];
            label tri2 = currentTri[2];

            point V;
            point Q;
            point R;

            if (tri0 == Vl)
            {
                V = localPoints[Vl];
                Q = localPoints[tri1];
                R = localPoints[tri2];
            }
            else if (tri1 == Vl)
            {
                V = localPoints[Vl];
                Q = localPoints[tri0];
                R = localPoints[tri2];
            }
            else
            {
                V = localPoints[Vl];
                Q = localPoints[tri0];
                R = localPoints[tri1];
            }

            // Edge vectors
            vector VQ = Q - V;
            vector VR = R - V;
            vector QR = R - Q;

            // Get angles
            scalar Va = getAngle(VQ, VR);
            scalar Ra = getAngle(VR, QR);
            scalar Qa = pi - (Va + Ra);

            /*
            if (badQuality(V, Q, R))
            {
                Info << "Bad triangle detected:\n"
                     << "Label: l = " << oneRingNeighborhood[T] << "\n"
                     << "Points:\nV" << V << "\nQ" << Q << "\nR" << R
                     << " \n"
                     << "Angles: v = " << Va*57.3 << ", q = " << Qa*57.3
                     << ", r = " << Ra*57.3 << "\n" << endl;
            }
            */

            // Check if non-obtuse in order to use the correct area metric
            if (Va < pi/2 && Qa < pi/2 && Ra < pi/2)
            {
                // Use Voronoi-area
                // Cotangent function 'cot' has to be defined locally since
                // it is not offered by OpenFOAM
                Amix += 0.125*(magSqr(VR)*cot(Qa) + magSqr(VQ)*cot(Ra));    
            }
            else if (Va >= pi/2)
            {
                // Obtuse angle at V, use half area
                // Use Heron's formula for now until problem
                // with vector approach is solved
                Amix += 0.5 * areaHeron(VQ, VR, QR);
            }
            else
            {
                // Obtuse angle not at V, use quarter area
                // Use Heron's formula for now until problem
                // with vector approach is solved
                Amix += 0.25 * areaHeron(VQ, VR, QR);
            }

            // Compute mean curvature normal contributions
            cn[Vl] += cot(Qa)*VR + cot(Ra)*VQ;
        }
        cn[Vl] = cn[Vl]/(2*Amix);
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
    scalar counter, curvature = 0;
    scalar minCur = 1000;
    scalar maxCur = 0;

    scalar tmp = 0;

    forAll(cn, V)
    {
        tmp = mag(cn[V]);        
        curvature += tmp;

        minCur = tmp < minCur ? tmp : minCur;
        maxCur = tmp > maxCur ? tmp : maxCur;

        counter++;
    }

    curvature = curvature / counter;

    Info << "Average curvature is: " << curvature<< endl;
    Info << "Maximum curvature is: " << maxCur<< endl;
    Info << "Minimum curvature is: " << minCur << endl; 

    return 0;
}

// ************************************************************************* //
