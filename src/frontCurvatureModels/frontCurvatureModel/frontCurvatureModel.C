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
    Foam::frontCurvatureModel

SourceFiles
    frontCurvatureModel.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Interface for the front curvature models. 

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


#include "frontCurvatureModel.H"
#include "dictionary.H"
#include "fvcGrad.H" 
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcAverage.H" 
#include "surfaceInterpolate.H" 
#include "addToRunTimeSelectionTable.H" 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontCurvatureModel, 0);
    defineRunTimeSelectionTable(frontCurvatureModel, Dictionary);
    addToRunTimeSelectionTable(frontCurvatureModel, frontCurvatureModel, Dictionary);

// * * * * * * * * * * * * * * *  Private member functions * * * * * * * * * //
void frontCurvatureModel::computeCurvature(const fvMesh& mesh, const triSurfaceFront& front) const
{
    const volScalarField& curvatureInputField = 
        mesh.lookupObject<volScalarField>(curvatureInputFieldName()); 

    const surfaceVectorField& Sf = mesh.Sf();

    //Cell gradient of alpha
    const volVectorField curvGrad(fvc::grad(curvatureInputField, "curvatureGradient"));

    // Interpolated face-gradient of alpha
    surfaceVectorField curvGradF(fvc::interpolate(curvGrad));

    // Hardcoded stabilization of the gradient to avoid floating point
    // exception.
    dimensionedScalar deltaN
    (
        "deltaN",
        curvatureInputField.dimensions() / dimLength, 
        SMALL 
    );

    // Face unit interface normal
    surfaceVectorField curvGradFhat(curvGradF /(mag(curvGradF) + deltaN));

    auto& cellCurvature = cellCurvatureTmp_.ref();
    cellCurvature = -fvc::div(curvGradFhat & Sf); 
}

void frontCurvatureModel::initializeCellCurvatureField(const fvMesh& mesh) const
{
    const Time& runTime = mesh.time();  

    cellCurvatureTmp_ = 
        tmp<volScalarField> 
        (
            new volScalarField
            (
                IOobject(
                    "cell_curvature", 
                    runTime.timeName(), 
                    mesh,
                    IOobject::NO_READ, 
                    IOobject::AUTO_WRITE
                ), 
                mesh, 
                dimensionedScalar(
                    "zero", 
                    dimless/dimLength, 
                    0.0
                )
            )
        );
}

bool frontCurvatureModel::curvatureNeedsUpdate(const fvMesh& mesh) const
{
    const auto& runTime = mesh.time();

    return (runTime.timeIndex() != lastTimeUpdated_);
}

void frontCurvatureModel::curvatureUpdated(const fvMesh& mesh) const
{
    lastTimeUpdated_ = mesh.time().timeIndex();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontCurvatureModel::frontCurvatureModel(const dictionary& configDict)
    :
        curvatureInputFieldName_{
            configDict.lookupOrDefault<word>("curvatureField", "curvatureField")
        },
        lastTimeUpdated_{-1},
        cellCurvatureTmp_{}
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

tmp<frontCurvatureModel>
frontCurvatureModel::New(const dictionary& configDict)
{
    const word name = configDict.lookup("type");

    DictionaryConstructorTable::iterator cstrIter =
        DictionaryConstructorTablePtr_->find(name);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "frontCurvatureModel::New(const word& name)"
        )   << "Unknown frontCurvatureModel type "
            << name << nl << nl
            << "Valid frontCurvatureModels are : " << endl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return tmp<frontCurvatureModel> (cstrIter()(configDict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> frontCurvatureModel::cellCurvature(
    const fvMesh& mesh, 
    const triSurfaceFront& frontMesh
) const
{
    if (cellCurvatureTmp_.empty())
    {
        initializeCellCurvatureField(mesh);
    }

    if (curvatureNeedsUpdate(mesh))
    {
        computeCurvature(mesh, frontMesh);
        curvatureUpdated(mesh); 
    }

    return cellCurvatureTmp_; 
}

tmp<surfaceScalarField> frontCurvatureModel::faceCurvature(
    const fvMesh& mesh, 
    const triSurfaceFront& frontMesh
) const
{
    const auto& cellCurvatureField = cellCurvature(mesh, frontMesh).ref();
    return fvc::interpolate(cellCurvatureField); 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
