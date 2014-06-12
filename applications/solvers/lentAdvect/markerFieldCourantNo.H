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
    Foam::diffuseInterfaceProperties

SourceFiles
    diffuseInterfaceProperties.C

Author
    Holger Marschall <marschall@csi.tu-darmstadt.de>
    All rights reserved.

Description
    Calculates and outputs the mean and maximum Courant Numbers.

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


scalar maxMarkerFieldCo
(
    readScalar(runTime.controlDict().lookup("maxMarkerFieldCo"))
);

scalar markerFieldCoNum = 0.0;
scalar meanMarkerFieldCoNum = 0.0;

if (mesh.nInternalFaces())
{
    scalarField sumPhi
    (
        pos(markerField - 0.01)*pos(0.99 - markerField)
       *fvc::surfaceSum(mag(phi))().internalField()
    );

    markerFieldCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    meanMarkerFieldCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Interface Courant Number mean: " << meanMarkerFieldCoNum
    << " max: " << markerFieldCoNum << endl;

// ************************************************************************* //
