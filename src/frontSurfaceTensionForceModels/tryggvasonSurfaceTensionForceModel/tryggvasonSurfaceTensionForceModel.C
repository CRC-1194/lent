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
    Foam::tryggvasonSurfaceTensionForceModel

SourceFiles
    tryggvasonSurfaceTensionForceModel.C

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


#include "tryggvasonSurfaceTensionForceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "lentInterpolation.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(tryggvasonSurfaceTensionForceModel, 0);
    addToRunTimeSelectionTable(frontSurfaceTensionForceModel, tryggvasonSurfaceTensionForceModel, Dictionary);

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * //

tryggvasonSurfaceTensionForceModel::tryggvasonSurfaceTensionForceModel(const dictionary& configDict)
    :
        frontSurfaceTensionForceModel(configDict) 
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<surfaceScalarField> tryggvasonSurfaceTensionForceModel::faceSurfaceTensionForce(
    const fvMesh& mesh,  
    const triSurfaceFront& frontMesh 
) const
{
    tmp<surfaceVectorField> SfHat = mesh.Sf() / mesh.magSf(); 
    return fvc::interpolate(cellSurfaceTensionForce(mesh,frontMesh)) & SfHat; 
}

tmp<volVectorField> tryggvasonSurfaceTensionForceModel::cellSurfaceTensionForce(
    const fvMesh& mesh,  
    const triSurfaceFront& frontMesh 
) const
{
    const Time& runTime = mesh.time(); 

    const dictionary& transportProperties = 
        runTime.lookupObject<dictionary>("transportProperties");

    tmp<volVectorField> fSigmaCellTmp(
        new volVectorField(
            IOobject(
                "fSigmaCell", 
                runTime.timeName(), 
                mesh, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh, 
            dimensionedVector(
                "zero", 
                dimForce / dimVolume, 
                vector(0,0,0)
            )
        )
    );

    // TODO: Store as data member and resize with topological change. TM.
    triSurfaceFrontPointVectorField fSigmaFront
    (
        IOobject(
            "fSigmaFront", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ), 
        frontMesh, 
        dimensionedVector(
            "zero", 
            dimForce / dimVolume, // Necessary for interpolation. TM. 
            vector(0.0,0.0,0.0)
        )
    );


    // Get necessary references
    const pointField& vertices = frontMesh.localPoints();
    const List<labelledTri>& triangles = frontMesh.localFaces();

    // Iterate over all triangles and compute the contribution to each of
    // its vertices
    forAll(triangles, Tl)
    {
        labelledTri T = triangles[Tl];

        // Get vertex labels
        label l0 = T[0];
        label l1 = T[1];
        label l2 = T[2];

        // Get vertices
        point v0 = vertices[l0];
        point v1 = vertices[l1];
        point v2 = vertices[l2];

        // Set edge vectors. They oriented in such a way that they follow
        // the rotational direction of the triangle normal vector. Since
        // the normal vectors are consistently defined (see lentFOAM paper) 
        // this ensures consistency in the direction of the contributions
        vector v01 = v1 - v0;
        vector v12 = v2 - v1;
        vector v20 = v0 - v2;

        // The order is inverted compared to the lentFOAM paper to factor in
        // the inverted direction of vector v20
        vector normal = v20 ^ v01;
        normal = normal / mag(normal);

        // Inverted order of the cross product is due to different labelling
        // compared to Tryggvason book. With outward facing normals, the cross
        // product wirtten here reults in an inward "pull" 
        fSigmaFront[l0] += 0.5 * v12 ^ normal;
        fSigmaFront[l1] += 0.5 * v20 ^ normal;
        fSigmaFront[l2] += 0.5 * v01 ^ normal;
    }

    const dimensionedScalar sigma = transportProperties.get<dimensionedScalar>("sigma");  
    fSigmaFront *= sigma; 

    lentInterpolation interpolation; // FIXME: Data member, RTS alternatives. TM.

    volVectorField& fSigmaCell = fSigmaCellTmp.ref(); 
    interpolation.interpolate(fSigmaFront,fSigmaCell); 

    // Avoid dimension check when dividing by mesh.V(); 
    const scalarField& V = mesh.V(); 
    forAll(fSigmaCell, cellI)
    {
        fSigmaCell[cellI] /= V[cellI]; 
    }

    return fSigmaCellTmp; 
}

tmp<fvMatrix<vector>> tryggvasonSurfaceTensionForceModel::surfaceTensionImplicitPart(
    const volVectorField&,
    const volScalarField&,
    const triSurfaceFront&
) const
{
    notImplemented("surfaceTensionImplicitPart()");

    return tmp<fvMatrix<vector>>{};
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
