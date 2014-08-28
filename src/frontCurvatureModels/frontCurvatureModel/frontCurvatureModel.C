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
#include "surfaceInterpolate.H" 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontCurvatureModel, 0);
    defineRunTimeSelectionTable(frontCurvatureModel, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontCurvatureModel::frontCurvatureModel(const dictionary& configDict, const Time& runTime)
    :
        runTime_(runTime), 
        meshName_(configDict.lookupOrDefault("meshName", word("mesh"))),
        mesh_(runTime_.lookupObject<fvMesh>(meshName_)),
        inputFieldName_(configDict.lookup("inputFieldName")), 
        inputField_(runTime_.lookupObject<volScalarField>(inputFieldName_)),
        delta_(readScalar(configDict.lookup("delta"))) 
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

tmp<frontCurvatureModel>
frontCurvatureModel::New(const dictionary& configDict, const Time& runTime)
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

    return tmp<frontCurvatureModel> (cstrIter()(configDict, runTime));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

frontCurvatureModel::~frontCurvatureModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> frontCurvatureModel::cellCurvature() const
{
    tmp<volScalarField> curvatureTmp(
        new volScalarField(
            IOobject(
                "cellCurvature", 
                runTime_.timeName(), 
                mesh_, 
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimless, 
                0
            )
        )
    );

    volScalarField& curvature  = curvatureTmp();  

    tmp<volVectorField> cellGrad = fvc::grad(inputField_); 

    surfaceVectorField faceGrad = fvc::interpolate(cellGrad); 

    surfaceVectorField faceGradNorm = (faceGrad / (mag(faceGrad) + delta_));

    tmp<surfaceScalarField> scalarFaceGradNormTmp = faceGradNorm & mesh_.Sf(); 

    curvature = -fvc::div(scalarFaceGradNormTmp()); 

    return curvatureTmp; 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
