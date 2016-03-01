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
    Foam::lentMarkerfieldTest

SourceFiles
    lentMarkerfieldTest.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    This class provides all tests related to the marker field model

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

#ifndef lentMarkerfieldTest_H
#define lentMarkerfieldTest_H

#include "fvCFD.H"

#include "triSurfaceFront.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FrontTracking
{

// Forward declaration of friend functions and operators
class lentMarkerfieldTest;

/*---------------------------------------------------------------------------*\
                         Class lentMarkerfieldTest Declaration
\*---------------------------------------------------------------------------*/

class lentMarkerfieldTest
{
    // Private data
    const volScalarField& markerField_;
    const triSurfaceFront& front_;
    scalar meshVol_;
    scalar markerfieldVol_;
    scalar frontVol_;
    point innerPoint_;

    // Private Member Functions
        void markerfieldVolume();
        void meshVolume();
        void frontVolume();

        void innerPoint();
        scalar tetrahedralVolume(const labelledTri& face,
                                 const Field<point>& vertices);

        // Foam default template stuff
        //- Disallow default bitwise copy construct
        lentMarkerfieldTest(const lentMarkerfieldTest&);

        //- Disallow default bitwise assignment
        void operator=(const lentMarkerfieldTest&);


public:

    // Constructors

        //- Construct from components
        lentMarkerfieldTest(const volScalarField& markerField,
                            const triSurfaceFront& front);


    //- Destructor
    ~lentMarkerfieldTest();


    // Member Functions

        // Tests
        bool boundedness() const;
        scalar globalVolume() const;
        void localVolume() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
