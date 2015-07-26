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
    Foam::frontExactCurvatureModel

SourceFiles
    frontExactCurvatureModel.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description

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

#include "frontExactCurvatureModel.H"
#include "fvcAverage.H"
#include "fvcDiv.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontExactCurvatureModel::frontExactCurvatureModel(const dictionary& configDict)
    :
        frontCurvatureModel(configDict), 
        write_(configDict.lookupOrDefault<Switch>("write", "off"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> frontExactCurvatureModel::cellCurvature(
    const fvMesh& mesh, 
    const triSurfaceMesh& frontMesh
) const
{
    const Time& runTime = mesh.time(); 

    tmp<volScalarField> cellCurvatureTmp( 
        new volScalarField(
            IOobject(
                "cellCurvatureExact", 
                runTime.timeName(), 
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            //FIXME: Error in the curvature dimensions. Fix it. TM.
            dimensionedScalar("zero", pow(dimLength, -1), 0) 
        )
    );

    volScalarField& cellCurvature = cellCurvatureTmp(); 

    const volVectorField& C = mesh.C(); 
    const surfaceVectorField& Cf = mesh.Cf(); 

    // Set internal cell centered curvature field.
    forAll(cellCurvature, I)
    {
        cellCurvature[I] = curvatureAtPoint(C[I]);
    }

    // Set the boundary cell centered curvature field
    forAll(cellCurvature.boundaryField(), I)
    {
        fvPatchScalarField& cellCurvatureBoundary = cellCurvature.boundaryField()[I];
        const fvsPatchVectorField& CfBoundary = Cf.boundaryField()[I];

        forAll(cellCurvatureBoundary, J)
        {
            cellCurvatureBoundary[J] = curvatureAtPoint(CfBoundary[J]);
        }
    }

    if (write_)
    {
        cellCurvature.write(); 
    }

    return cellCurvatureTmp;
}

tmp<surfaceScalarField> frontExactCurvatureModel::faceCurvature(
    const fvMesh& mesh, 
    const triSurfaceMesh& frontMesh
) const
{
    const Time& runTime = mesh.time(); 
    tmp<surfaceScalarField> surfaceCurvatureTmp( 
        new surfaceScalarField(
            IOobject(
                "faceCurvatureExact", 
                runTime.timeName(), 
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", pow(dimLength, -2), 0) 
        )
    );

    surfaceScalarField& surfaceCurvature = surfaceCurvatureTmp(); 

    const surfaceVectorField& Cf = mesh.Cf(); 

    // Set internal surface centered curvature field.
    forAll(surfaceCurvature, I)
    {
        surfaceCurvature[I] = curvatureAtPoint(Cf[I]);
    }

    // Set the boundary surface centered curvature field
    forAll(surfaceCurvature.boundaryField(), I)
    {
        fvsPatchScalarField& surfaceCurvatureBoundary = surfaceCurvature.boundaryField()[I];
        const fvsPatchVectorField& CfBoundary = Cf.boundaryField()[I];

        forAll(surfaceCurvatureBoundary, J)
        {
            surfaceCurvatureBoundary[J] = curvatureAtPoint(CfBoundary[J]);
        }
    }

    if (write_)
    {
        surfaceCurvature.write(); 
    }

    return surfaceCurvatureTmp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
