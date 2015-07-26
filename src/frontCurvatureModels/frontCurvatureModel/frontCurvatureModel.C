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
#include "fvcAverage.H" 
#include "surfaceInterpolate.H" 
#include "addToRunTimeSelectionTable.H" 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontCurvatureModel, 0);
    defineRunTimeSelectionTable(frontCurvatureModel, Dictionary);
    addToRunTimeSelectionTable(frontCurvatureModel, frontCurvatureModel, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontCurvatureModel::frontCurvatureModel(const dictionary& configDict)
    :
        curvatureInputFieldName_(
            configDict.lookupOrDefault<word>("curvatureField", "curvatureField")
        )
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
    const triSurfaceMesh& frontMesh
) const
{
    const volScalarField& curvatureInputField = 
        mesh.lookupObject<volScalarField>(curvatureInputFieldName()); 

    const surfaceVectorField& Sf = mesh.Sf();

    //Cell gradient of alpha
    const volVectorField curvGrad(fvc::grad(curvatureInputField, "curvatureGradient"));

    // TODO: uncomment
    if (debug)
    {
        curvGrad.write(); 
    }

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

    if (debug)
    {
        volScalarField curvature = fvc::div(curvGradFhat & Sf); 
        curvature.rename("curvature"); 
        curvature.write(); 
    }

    // FIXME: Check the sign of the curvature in the pU coupling files. TM
    return fvc::div(curvGradFhat & Sf); 
    
    // TODO: Analyze the difference between the two options.
    //return -fvc::div(
        //fvc::grad(curvatureInputField) / (mag(fvc::grad(curvatureInputField)) + deltaN)
    //); 
}

tmp<surfaceScalarField> frontCurvatureModel::faceCurvature(
    const fvMesh& mesh, 
    const triSurfaceMesh& frontMesh
) const
{
    return fvc::interpolate(cellCurvature(mesh,frontMesh)); 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
