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
    Foam::csfSurfaceTensionForceModel

SourceFiles
    csfSurfaceTensionForceModel.C

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
       
:qa
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#include "csfSurfaceTensionForceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "fvcReconstruct.H"
#include "surfaceInterpolate.H"
#include "surfaceFields.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(csfSurfaceTensionForceModel, 0);
    addToRunTimeSelectionTable(frontSurfaceTensionForceModel, csfSurfaceTensionForceModel, Dictionary);

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * //

csfSurfaceTensionForceModel::csfSurfaceTensionForceModel(const dictionary& configDict)
    :
        curvatureFieldName_(configDict.lookup("curvatureField")), 
        filterFieldName_(configDict.lookup("filterField"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> csfSurfaceTensionForceModel::cellCurvature(
    const volScalarField& curvatureField 
) const
{
    const fvMesh& mesh = curvatureField.mesh(); 
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(curvatureField, "nHat"));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Hardcoded stabilization of the gradient to avoid floating point
    // exception.
    dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/pow(average(curvatureField.mesh().V()), 1.0/3.0)
    );

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN));

    return -fvc::div(nHatfv & Sf); 
} 

tmp<surfaceScalarField> csfSurfaceTensionForceModel::faceSurfaceTensionForce(
    const fvMesh& mesh,  
    const triSurfaceMesh& frontMesh 
) const
{
    // Read off the surface tension coefficient from the transport properties 
    // dictionary. 
    const Time& runTime = mesh.time(); 
    const dictionary& transportProperties = 
        runTime.lookupObject<dictionary>("transportProperties");
    const dimensionedScalar sigma = transportProperties.lookup("sigma");  
    
    const volScalarField& curvatureField = 
        mesh.lookupObject<volScalarField>(curvatureFieldName_); 
    
    const volScalarField& filterField = 
        mesh.lookupObject<volScalarField>(filterFieldName_); 

    return fvc::interpolate(
               sigma * cellCurvature(curvatureField)
           ) * 
           fvc::snGrad(filterField);
}

tmp<volVectorField> csfSurfaceTensionForceModel::cellSurfaceTensionForce(
    const fvMesh& mesh,  
    const triSurfaceMesh& frontMesh 
) const
{
    return fvc::reconstruct(faceSurfaceTensionForce(mesh, frontMesh));  
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
