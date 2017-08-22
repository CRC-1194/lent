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
    Foam::frontCompactDivGradCurvatureModel

SourceFiles
    frontCurvatureModel.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description
    Compact version of the divGrad curvature Model. Uses the front-mesh
    communication to transfer the curvature values of interface cells
    to all cells in the narrow band

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

#include "frontCompactDivGradCurvatureModel.H"

#include "addToRunTimeSelectionTable.H"
#include "triSurfaceFields.H"

#include "lentCommunication.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontCompactDivGradCurvatureModel, 0);
    addToRunTimeSelectionTable(frontCurvatureModel, frontCompactDivGradCurvatureModel, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontCompactDivGradCurvatureModel::frontCompactDivGradCurvatureModel(const dictionary& configDict)
    :
        frontCurvatureModel(configDict)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> frontCompactDivGradCurvatureModel::cellCurvature(
    const fvMesh& mesh, 
    const triSurfaceFront& frontMesh
) const
{
    const Time& runTime = mesh.time();  

    triSurfaceScalarField curvatureBuffer
    {
        IOobject(
            "curvatureBuffer",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        frontMesh,
        dimensionedScalar{
            "zero",
            dimless/dimLength,
            0
        }
    };

    const auto& communication = mesh.lookupObject<lentCommunication>(
                                    lentCommunication::registeredName(frontMesh, mesh)
                                ); 
    const auto& triangleToCell = communication.triangleToCell();
    auto cellCurvatureFieldTmp = frontCurvatureModel::cellCurvature(mesh, frontMesh);
    auto& cellCurvatureField = cellCurvatureFieldTmp.ref();


    const volScalarField& curvatureInputField = 
        mesh.lookupObject<volScalarField>(curvatureInputFieldName()); 

    forAll(triangleToCell, index)
    {
        auto temp = cellCurvatureField[triangleToCell[index]];
        curvatureBuffer[index] = sign(temp) * 2.0 / (2.0/mag(temp)
                                            - curvatureInputField[triangleToCell[index]]);
    }
    cellCurvatureField *= 0.0;

    const auto& cellToTriangle = communication.cellsTriangleNearest();

    forAll(cellToTriangle, I)
    {
        const auto& hitObject = cellToTriangle[I];
        
        if (hitObject.hit())
        {
            cellCurvatureField[I] = curvatureBuffer[hitObject.index()];
        }
    }

    return cellCurvatureFieldTmp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
