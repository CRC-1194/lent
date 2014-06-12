/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::lentTests

Description
    Test functions for the LENT algorithm.

SourceFiles


Authors
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#ifndef lentTests_H
#define lentTests_H

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

    bool normalsAreConsistent = true; 


    forAll (edgeFaces, I)
    {
        const vector& n0 = faceNormals[edgeFaces[I][0]]; 
        const vector& n1 = faceNormals[edgeFaces[I][1]]; 

        if ((n0 & n1) < 0)
        {
            normalsAreConsistent = false; 
        }
    }
    return normalsAreConsistent;
}

class triSurfaceSizeChange 
{
    int lastSize_ = 0; 

    const triSurface& t_; 

    double allowedChange_; 

    public: 

        triSurfaceSizeChange(const triSurface& t, double allowedChange=0.2) 
            : 
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

// Classes for testing the application code.  

class lentTests : public ::testing::Test 
{
    protected:
        template<typename Front>
        void TEST_NORMAL_CONSISTENCY(Front const & front) 
        {
            ASSERT_TRUE(triSurfaceNormalsAreConsistent(front));
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test 

} // End namespace FrontTracking 

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
