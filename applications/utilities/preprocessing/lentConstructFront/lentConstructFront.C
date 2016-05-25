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

Authors
    Tobias  Tolle  tolle@csi.tu-darmstadt.de

Description
    Create a front by intersecting an analytical surface with the mesh
    and subsequent per-cell triangulaion. It is written to "front/front.stl".

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

#include "analyticalSurface.H"
#include "frontConstructor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace FrontTracking;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    IOdictionary lentControlDict(
        IOobject(
           "lentSolution",
           "system",
           mesh.time(),
           IOobject::MUST_READ_IF_MODIFIED,
           IOobject::AUTO_WRITE
        )
    );

    Info << "Constructing front..." << endl;

    tmp<analyticalSurface> analyticalSurfaceTmp(
        analyticalSurface::New(lentControlDict.subDict("analyticalSurface"))
    );

    frontConstructor buildFront(analyticalSurfaceTmp, mesh);

    triSurface front = buildFront.createTriSurface();

    front.write("front/front.stl");

    Info << "Done!" << endl;

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
