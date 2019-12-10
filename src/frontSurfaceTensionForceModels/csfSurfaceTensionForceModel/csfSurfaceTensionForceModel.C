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
#include "fvcAverage.H"
#include "fvmLaplacian.H"
#include "gaussLaplacianScheme.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(csfSurfaceTensionForceModel, 0);
    addToRunTimeSelectionTable(frontSurfaceTensionForceModel, csfSurfaceTensionForceModel, Dictionary);

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * //

csfSurfaceTensionForceModel::csfSurfaceTensionForceModel(const dictionary& configDict)
    :
        curvatureBasedSurfaceTensionForceModel(configDict) 
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<surfaceScalarField> csfSurfaceTensionForceModel::faceSurfaceTensionForce(
    const fvMesh& mesh,  
    const triSurfaceFront& frontMesh 
) const
{
    const Time& runTime = mesh.time(); 

    const dictionary& transportProperties = 
        runTime.lookupObject<dictionary>("transportProperties");

    const dimensionedScalar sigma = transportProperties.lookup("sigma");  
    
    const volScalarField& filterField = mesh.lookupObject<volScalarField>(filterFieldName());     

    auto faceCurvatureFieldPtr = faceCurvature(mesh, frontMesh);
    const auto& faceCurvatureField = *faceCurvatureFieldPtr;
    // FIXME: for a reason I do not yet understand, the orientation of the
    // surface tension field on the cell faces is 'unoriented' and will cause
    // an error when being added/subtracted to other face force fields.
    // For now, just manually set the resulting field to 'oriented' as the
    // original snGrad(alpha) field (TT)
    auto faceSurfaceTensionTmp = sigma * faceCurvatureField * fvc::snGrad(filterField);
    faceSurfaceTensionTmp.ref().setOriented();
    
    return faceSurfaceTensionTmp;
}

tmp<volVectorField> csfSurfaceTensionForceModel::cellSurfaceTensionForce(
    const fvMesh& mesh,  
    const triSurfaceFront& frontMesh 
) const
{
    return fvc::reconstruct(faceSurfaceTensionForce(mesh, frontMesh) * mesh.magSf());  
}

tmp<fvMatrix<vector>> csfSurfaceTensionForceModel::surfaceTensionImplicitPart(
    const volVectorField& velocity,
    const volScalarField& markerField,
    const triSurfaceFront& front
) const
{
    // Lookup material properties
    const auto& runTime = velocity.mesh().time();
    const auto& mesh = velocity.mesh();

    const dictionary& transportProperties = 
        runTime.lookupObject<dictionary>("transportProperties");

    const dimensionedScalar sigma = transportProperties.lookup("sigma");  

    auto sigmaDeltaT = sigma*runTime.deltaT();

    // Below is a first implementation of a fully implicit discretization
    // approach for the Laplace-Beltrami operator in the level set context (TT).
    /*
    const auto& phi = mesh.lookupObject<volScalarField>("signedDistance");
    const auto& Sf = mesh.Sf();
    const auto& nf = mesh.Sf()/mesh.magSf();
    const auto& delta = mesh.deltaCoeffs();

    surfaceVectorField grad_phi_f{fvc::interpolate(fvc::grad(phi))};

    surfaceScalarField gammaSf{sigmaDeltaT*(grad_phi_f & Sf)*(grad_phi_f & nf) / 
            (magSqr(grad_phi_f) + SMALL)};
    surfaceScalarField gammaFullSf{sigmaDeltaT*mesh.magSf()};

    // Filter flux field gammaSf here
    surfaceScalarField filter{fvc::interpolate(mag(fvc::grad(markerField)))};

    gammaFullSf.dimensions() *= pow(dimLength, -1);
    gammaSf.dimensions() *= pow(dimLength, -1);

    forAll(filter, I)
    {
        gammaSf[I] *= filter[I];
        gammaFullSf[I] *= filter[I];
    }
    
    fv::gaussLaplacianScheme<vector,scalar> gaussLaplace{mesh};

    return gaussLaplace.fvmLaplacianUncorrected(gammaFullSf, delta, velocity) - gaussLaplace.fvmLaplacianUncorrected(gammaSf, delta, velocity);
    */

    // Lookup interface properties
    auto curvaturePtr = cellCurvature(velocity.mesh(), front);
    auto interfaceNormalPtr = curvatureModelRef().cellInterfaceNormals(velocity.mesh(), front);
    const auto& curvature = *curvaturePtr;
    const auto& normals = *interfaceNormalPtr;
    
    // Below: normal calculation consistent with explicit surface tension part (TT)
    //dimensionedScalar dSmall{"SMALL", pow(dimLength, -1), SMALL};
    //auto normalsTmp = fvc::grad(markerField)/(mag(fvc::grad(markerField)) + dSmall);
    //const auto& normals = normalsTmp.ref();

    // Define Laplace-Beltrami of velocity as the full Laplace-operator
    // (implicit) and subtract the normal part (explicit)
    auto gradUTmp = fvc::grad(velocity);
    const auto& gradU = gradUTmp.ref();

    auto normalLaplacian = fvc::div((normals&gradU)*normals)
            - curvature*((gradU - ((normals&gradU)*normals))&normals);

    return (fvm::laplacian(sigmaDeltaT*mag(fvc::grad(markerField)), velocity)
                - sigmaDeltaT*mag(fvc::grad(markerField))*normalLaplacian);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
