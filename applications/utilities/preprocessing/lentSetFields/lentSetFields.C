/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Application
    lentSetFields

Authors
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

Description
    Pre-processing application that sets two signed distance fields and a 
    heaviside marker field using an input surface mesh in the STL format. 
    Used as a pre-processing application for the LENT algorithm. 


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "interfaceProperties.H"
#include "incompressibleTwoPhaseMixture.H"

#include "lentMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace FrontTracking;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    triSurfaceFront front(
        IOobject(
            "front.stl",
            "front",
            runTime, 
            IOobject::MUST_READ, 
            IOobject::AUTO_WRITE
        )
    );

    lentMethod lent(front, mesh); 

    Info << "Calculating the search distance fields..."; 
    lent.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);
    Info << "Done." << endl;

    Info << "Calculating the signed distance fields..."; 
    lent.calcSignedDistances(
        signedDistance, 
        pointSignedDistance, 
        searchDistanceSqr, 
        pointSearchDistanceSqr,
        front
    ); 
    Info << "Done." << endl;

    Info << "Calculating the heaviside field..."; 
    lent.calcHeaviside(heaviside, signedDistance, searchDistanceSqr); 
    Info << "Done." << endl;
    
    runTime.writeNow(); 

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
