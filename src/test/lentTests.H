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

#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
    namespace FrontTracking {
        namespace Test {

bool normalsAreConsistent(const triSurface& front)
{
    const labelListList& edgeFaces = front.edgeFaces();
    const vectorField& faceNormals = front.faceNormals();

    forAll (edgeFaces, I)
    {
        if (edgeFaces[I].size() > 1)
        {
            const vector& n0 = faceNormals[edgeFaces[I][0]];
            const vector& n1 = faceNormals[edgeFaces[I][1]];

            if ((n0 & n1) < 0)
            {
                return false;
            }
        }
    }
    return true;
}

bool normalsAreInconsistent(const triSurface& front)
{
    return ! normalsAreConsistent(front); 
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

} // End namespace FrontTracking

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
