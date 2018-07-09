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
    Foam::analyticalSurface

SourceFiles
    analyticalSurface.C

Author
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    Abstract base class for definition of (simple) analytical surfaces.
    Intended for directly constructing a front or setting an exact signed
    distance field.

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

#include "analyticalSurface.H"

#include "pointMesh.H"
#include "pointPatchField.H"
#include "pointFieldsFwd.H"
#include "valuePointPatchField.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalSurface, 0);
    defineRunTimeSelectionTable(analyticalSurface, Dictionary);
    

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
OFstream analyticalSurface::outputStream(const word& fileName) const
{
    return OFstream{fileName, IOstream::ASCII, IOstream::currentVersion, IOstream::UNCOMPRESSED, true};
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalSurface::analyticalSurface(const dictionary& configDict)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
tmp<analyticalSurface> analyticalSurface::New(const dictionary& configDict)
{
    const word name = configDict.lookup("type");

    DictionaryConstructorTable::iterator cstrIter =
        DictionaryConstructorTablePtr_->find(name);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "analyticalSurface::New(const word& name)"
        )   << "Unknown analyticalSurface type "
            << name << nl << nl
            << "Valid analyticalSurface are : " << endl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return tmp<analyticalSurface> (cstrIter()(configDict));
}

// * * * * * * * * * * * * * Public member functions * * * * * * * * * * * * //
void analyticalSurface::setDistance(volScalarField& signedDistance) const
{
    const auto& cellCentres = signedDistance.mesh().C();

    setDistance<volScalarField, volVectorField>(signedDistance, cellCentres);
}

void analyticalSurface::setDistance(surfaceScalarField& signedDistance) const
{
    const auto& faceCentres = signedDistance.mesh().Cf();

    setDistance<surfaceScalarField, surfaceVectorField>(signedDistance, faceCentres);

    // Set the boundary values
    auto& boundaryPatches = signedDistance.boundaryFieldRef();

    forAll(boundaryPatches, I)
    {
        auto& bPatch = boundaryPatches[I];

        setDistance<fvsPatchField<scalar>, surfaceVectorField>(bPatch, faceCentres);
    }
}

void analyticalSurface::setDistance(pointScalarField& signedDistance) const
{
    const auto& meshVertices = signedDistance.mesh()().points();

    setDistance<pointScalarField, pointField>(signedDistance, meshVertices);

    // TODO: Does the code above really set the signed distance for all
    // vertices of the fvMesh including boundary points? (TT)
}

void analyticalSurface::moveFrontToSurface(triSurfaceFront& front) const
{
    auto& points = const_cast<pointField&>(front.points());

    for (auto& vertex : points)
    {
        vertex = this->normalProjectionToSurface(vertex);
    }

    front.clearGeom();
}

void analyticalSurface::makeNormalOrientationConsistent(triSurfaceFront& front, const bool outwardOrientation) const
{
    // enforce recomputation of front normals to ensure they are up-to-date
    // with the front vertex positions
    front.clearGeom();

    scalar orientation = 1.0;

    if (!outwardOrientation)
    {
        orientation = -1.0;
    }

    List<labelledTri>& triangles = static_cast<List<labelledTri>& > (front);
    auto& triangleNormals = const_cast<pointField&>(front.faceNormals());
    const auto& faceCentres = front.Cf();

    forAll(triangles, I)
    {
        auto normal = orientation*normalToPoint(faceCentres[I]);
        if ((normal & triangleNormals[I]) < 0.0)
        {
            triangles[I].flip();
            triangleNormals[I] *= -1.0;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
