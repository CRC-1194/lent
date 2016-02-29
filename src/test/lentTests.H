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
    Foam::lentTests

SourceFiles

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Test functions for the LENT algorithm.

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


#ifndef lentTests_H
#define lentTests_H

#include "fvCFD.H"
#include "triSurface.H"
#include "gtest.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
    namespace FrontTracking {
        namespace Test {

bool triSurfaceNormalsAreConsistent(const triSurface& front)
{
    const labelListList& edgeFaces = front.edgeFaces();
    const vectorField& faceNormals = front.faceNormals();

    forAll (edgeFaces, I)
    {
        const vector& n0 = faceNormals[edgeFaces[I][0]];
        const vector& n1 = faceNormals[edgeFaces[I][1]];

        if ((n0 & n1) < 0)
        {
            return false;
        }
    }
    return true;
}

class triSurfaceSizeChange
{
    int lastSize_;

    const triSurface& t_;

    double allowedChange_;

    public:

        triSurfaceSizeChange(const triSurface& t, double allowedChange=0.2)
            :
                lastSize_(0),
                t_(t),
                allowedChange_(allowedChange)
        {}

        void operator()()
        {
            label triSurfaceSize = t_.localFaces().size();

            if (lastSize_ > 0)
            {
                if (((double(lastSize_) - triSurfaceSize) / lastSize_) > allowedChange_)
                {
                    Info << "TEST : triSurfaceSizeChange " << allowedChange_ * 100
                        << " % : FAIL " << endl;
                }
                else
                {
                    Info << "TEST : triSurfaceSizeChange : PASS " << endl;
                }
            }

            lastSize_ = triSurfaceSize;
        }
};

// Auxilliary functions for testing marker field models
bool isBounded(const volScalarField& markerField, scalar lowerBound,
               scalar upperBound)
{
    bool bounded = true;

    forAll(markerField, I)
    {
        if (markerField[I] < lowerBound || markerField[I] > upperBound)
        {
            bounded = false;
        }
    }

    return bounded;
}

scalar markerFieldVolume(const volScalarField& markerField)
{
    scalar volume = 0.0;

    const fvMesh& mesh = markerField.mesh();

    forAll(markerField, I)
    {
        volume += markerField[I] * mesh.V()[I];
    }

    return volume;
}

scalar meshVolume(const volScalarField& markerField)
{
    dimensionedScalar volume(
        "zero",
        pow(dimLength, 3),
        0.0
    );

    const fvMesh& mesh = markerField.mesh();

    volume = sum(mesh.V());

    return volume.value();
}

point innerPoint(const triSurfaceFront& front)
{
    // FIXME: works only for convex shapes
    point lowerRef(GREAT, GREAT, GREAT);
    point upperRef(0.0, 0.0, 0.0);

    const Field<point>& vertices = front.localPoints();

    forAll(vertices, I)
    {
        forAll(vertices[I], K)
        {
            if (vertices[I][K] < lowerRef[K])
            {
                lowerRef[K] = vertices[I][K];
            }

            if (vertices[I][K] > upperRef[K])
            {
                upperRef[K] = vertices[I][K];
            }
        }
    }
    
    return 0.5*(lowerRef + upperRef);
}

scalar tetrahedralVolume(const labelledTri& face, const Field<point>& vertices,
                         point innerPoint)
{
    // use property of scalar triple product
    return fabs(1.0/6.0 * ((vertices[face[1]] - vertices[face[0]]) ^ 
                    (vertices[face[2]] - vertices[face[0]]))
                    & (innerPoint - vertices[face[0]]));
}

scalar frontVolume(const triSurfaceFront& front)
{
    scalar volume = 0.0;

    point pseudoCentre = innerPoint(front);
    const List<labelledTri>& faces = front.localFaces();
    const Field<point>& vertices = front.localPoints();

    forAll(faces, I)
    {
        volume += tetrahedralVolume(faces[I], vertices, pseudoCentre);
    }

    return volume;
}

scalar relativeVolumeDifference(const triSurfaceFront& front,
                                const volScalarField& markerField)
{
    scalar meshVol = meshVolume(markerField);
    scalar markerFieldVol = markerFieldVolume(markerField);
    scalar frontVol = frontVolume(front);

    // Ensure correct comparison of volumes --> markerField can assume 1 
    // inside or outside of the front
    if (   fabs(meshVol - markerFieldVol - frontVol)
         < fabs(markerFieldVol - frontVol))
    {
        return fabs(meshVol - markerFieldVol - frontVol) / frontVol;
    }
    else
    {
        return fabs(markerFieldVol - frontVol) / frontVol;
    }
}

// Classes for testing the application code.

class lentTests : public ::testing::Test
{
    protected:
        template<typename Front>
        void TEST_NORMAL_CONSISTENCY(Front const & front)
        {
            ASSERT_TRUE(triSurfaceNormalsAreConsistent(front));
        }

        // Marker field model tests
        void TEST_MARKERFIELD_BOUNDS(const volScalarField& markerField)
        {
            ASSERT_TRUE(isBounded(markerField, 0.0, 1.0));
        }

        // This test is bound to fail...
        void TEST_MARKERFIELD_GLOBAL_VOLUME(const volScalarField& markerField,
                                            const triSurfaceFront& front)
        {
            EXPECT_EQ(relativeVolumeDifference(front, markerField), 0.0);
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

} // End namespace FrontTracking

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
