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

Global
    heavisideCourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar maxHeavisideCo(
    readScalar(runTime.controlDict().lookup("maxHeavisideCo"))
);

scalar heavisideCoNum = 0.0;
scalar meanHeavisideCoNum = 0.0;

if (mesh.nInternalFaces()) {
    scalarField sumPhi(
        pos(heaviside - 0.01)*pos(0.99 - heaviside)
       *fvc::surfaceSum(mag(phi))().internalField()
    );

    heavisideCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    meanHeavisideCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Interface Courant Number mean: " << meanHeavisideCoNum
    << " max: " << heavisideCoNum << endl;

// ************************************************************************* //
