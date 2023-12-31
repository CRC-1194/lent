/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
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
    Foam::FrontTracking::analyticalSurface

SourceFiles
    analyticalSurfaceI.H
    analyticalSurface.C

Authors
    Tomislav Maric (maric@mma.tu-darmstadt.de)
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Abstract base class for definition of (simple) analytical surfaces.
    Intended for directly constructing a front or setting an exact signed
    distance field.

Affiliations:
    Mathematical Modeling and Analysis Institute, Mathematics Department, 
    TU Darmstadt, Germany

Funding:
    German Research Foundation (DFG) - Project-ID 265191195 - SFB 1194

    German Research Foundation (DFG) - Project-ID MA 8465/1-1, 
    Initiation of International Collaboration 
    "Hybrid Level Set / Front Tracking methods for simulating 
    multiphase flows in geometrically complex systems"

\*---------------------------------------------------------------------------*/

#include "analyticalSurface.H"

#include "pointMesh.H"
#include "pointPatchField.H"
#include "pointFieldsFwd.H"
#include "valuePointPatchField.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalSurface, 0);
    defineRunTimeSelectionTable(analyticalSurface, Dictionary)
    

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
OFstream analyticalSurface::outputStream(const word& fileName) const
{
    return OFstream{fileName, IOstream::ASCII, IOstream::currentVersion, IOstream::UNCOMPRESSED, true};
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalSurface::analyticalSurface(const dictionary&)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
tmp<analyticalSurface> analyticalSurface::New(const dictionary& configDict)
{
    const word name = configDict.get<word>("type");

    // Find the constructor pointer for the model in the constructor table.
    auto* ctorPtr = DictionaryConstructorTable(name);

    // If the constructor pointer is not found in the table.
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            configDict,
            "analyticalSurface",
            name,
            *DictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    // Construct the model and return the autoPtr to the object.
    return tmp<analyticalSurface> (ctorPtr(configDict));
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

void analyticalSurface::moveFrontToSurface(triSurface& front) const
{
    auto& points = const_cast<pointField&>(front.points());

    for (auto& vertex : points)
    {
        vertex = this->normalProjectionToSurface(vertex);
    }

    front.clearGeom();
}

void analyticalSurface::makeNormalOrientationConsistent(triSurface& front, const bool outwardOrientation) const
{
    // Enforce recomputation of front normals to ensure they are up-to-date
    // with the front vertex positions (TT)
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
        }
    }

    // Delete all deman driven data to ensure correctness and consistency (TT)
    front.clearOut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
