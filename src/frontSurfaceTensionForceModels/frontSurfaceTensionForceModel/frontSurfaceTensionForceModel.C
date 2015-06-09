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
    Foam::frontSurfaceTensionForceModel

SourceFiles
    frontSurfaceTensionForceModel.C

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


#include "frontSurfaceTensionForceModel.H"
#include "dictionary.H"
#include "fvcGrad.H" 
#include "fvcDiv.H"
#include "fvcAverage.H" 
#include "surfaceInterpolate.H" 
#include "addToRunTimeSelectionTable.H" 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontSurfaceTensionForceModel, 0);
    defineRunTimeSelectionTable(frontSurfaceTensionForceModel, Dictionary);
    addToRunTimeSelectionTable(frontSurfaceTensionForceModel, frontSurfaceTensionForceModel, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontSurfaceTensionForceModel::frontSurfaceTensionForceModel(
    const dictionary& configDict, 
    const Time& runTime
)
    :
        runTime_(runTime), 
        meshName_(configDict.lookupOrDefault("meshName", word("region0"))),
        mesh_(runTime_.lookupObject<fvMesh>(meshName_))
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

tmp<frontSurfaceTensionForceModel>
frontSurfaceTensionForceModel::New(const dictionary& configDict, const Time& runTime)
{
    const word name = configDict.lookup("type");

    DictionaryConstructorTable::iterator cstrIter =
        DictionaryConstructorTablePtr_->find(name);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "frontSurfaceTensionForceModel::New(const word& name)"
        )   << "Unknown frontSurfaceTensionForceModel type "
            << name << nl << nl
            << "Valid frontSurfaceTensionForceModels are : " << endl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return tmp<frontSurfaceTensionForceModel> (cstrIter()(configDict, runTime));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<surfaceScalarField> frontSurfaceTensionForceModel::faceSurfaceTensionForce() const
{
    tmp<surfaceScalarField> faceSurfaceTensionForceTmp(
        new surfaceScalarField ( 
            IOobject (
                runTime_.timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimForce / dimVolume,  
                scalar(0) 
            ) 
        )
    );

    return faceSurfaceTensionForceTmp;  
}

tmp<volVectorField> frontSurfaceTensionForceModel::cellSurfaceTensionForce() const
{
    tmp<volVectorField> cellSurcellTensionForceTmp(
        new volVectorField ( 
            IOobject (
                runTime_.timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                dimForce / dimVolume,  
                vector(0,0,0)
            ) 
        )
    );

    return cellSurcellTensionForceTmp;  
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
