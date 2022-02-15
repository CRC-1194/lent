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
    Foam::FrontTracking::surfaceSet

SourceFiles
    surfaceSetI.H
    surfaceSet.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    This class stores pointers to objects of type 'analyticalSurface'
    and set operations. Furthermore, it can set the resulting level set
    field to vol-, surface- or pointScalarFields.

    Furthermore, set operations are provided here in terms of overloaded
    operators:
        A + B : union of A and B
        A / B : intersection of A and B
        A - B : complement of B in A
        A % B : symmetric difference of A and B

    See also: https://en.wikipedia.org/wiki/Set_(mathematics)#Basic_operations

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

#include "pointMesh.H"
#include "pointFieldsFwd.H"

#include "surfaceSet.H"

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void surfaceSet::resizeDistanceBuffer() const
{
    if (signedDistanceBuffer_.size() != orderedSurfaces_.size())
    {
        signedDistanceBuffer_.resize(orderedSurfaces_.size());
    }
}

void surfaceSet::computeDistanceBuffer(const point& P) const
{
    for (unsigned int I = 0; I != orderedSurfaces_.size(); ++I)
    {
        signedDistanceBuffer_[I] = orderedSurfaces_[I]->signedDistance(P);
    }
}

scalar surfaceSet::computeLevelSetSign() const
{
    for (unsigned int I = 0; I != orderedOperations_.size(); ++I)
    {
        if (isOutsidePoint(I))
        {
            return 1.0;
        }
    }

    return -1.0;
}

scalar surfaceSet::computelevelSetValue() const
{
    scalar lSetValue{GREAT};

    for (const auto& distance : signedDistanceBuffer_)
    {
        if (mag(distance) < lSetValue)
        {
            lSetValue = mag(distance);
        }
    }

    return lSetValue;
}

bool surfaceSet::isOutsidePoint(const label& operationID) const
{
    const auto& operation = orderedOperations_[operationID];
    const auto& distA = signedDistanceBuffer_[operationID];
    const auto& distB = signedDistanceBuffer_[operationID + 1];

    if (operation == setOperation::unite)
    {
        return !pointIsInUnion(distA, distB);
    }
    else if (operation == setOperation::intersect)
    {
        return !pointIsInIntersection(distA, distB);
    }
    else if (operation == setOperation::complement)
    {
        return !pointIsInComplement(distA, distB);
    }
    else if (operation == setOperation::symm_diff)
    {
        return !pointIsInSymmDiff(distA, distB);
    }

    // This should never be reached (TT)
    return true;
}

bool surfaceSet::pointIsInUnion(const scalar& distA, const scalar& distB) const
{
    return (distA < 0.0) || (distB < 0.0);
}

bool surfaceSet::pointIsInIntersection(const scalar& distA, const scalar& distB) const
{
    return (distA < 0.0) && (distB < 0.0);
}

bool surfaceSet::pointIsInComplement(const scalar& distA, const scalar& distB) const
{
    return (distA < 0.0) && (distB > 0.0);
}

bool surfaceSet::pointIsInSymmDiff(const scalar& distA, const scalar& distB) const
{
    return pointIsInComplement(distA, distB) || pointIsInComplement(distB, distA);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
surfaceSet::surfaceSet(const analyticalSurface& surface)
:
    orderedSurfaces_{},
    orderedOperations_{} 
{
    analyticalSurface* rawPtr = const_cast<analyticalSurface*>(&surface);
    orderedSurfaces_.push_back(rawPtr);
}


surfaceSet::surfaceSet(const surfaceSet& rhs)
:
    orderedSurfaces_{rhs.orderedSurfaces_},
    orderedOperations_{rhs.orderedOperations_}
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar surfaceSet::levelSetValue(const point& P) const
{
    computeDistanceBuffer(P);
    
    return computeLevelSetSign()*computelevelSetValue();
}

void surfaceSet::setLevelSet(volScalarField& levelSetField) const
{
    const auto& cellCentres = levelSetField.mesh().C();

    setLevelSet<volScalarField, volVectorField>(levelSetField, cellCentres);
}

void surfaceSet::setLevelSet(surfaceScalarField& levelSetField) const
{
    const auto& faceCentres = levelSetField.mesh().Cf();

    setLevelSet<surfaceScalarField, surfaceVectorField>(levelSetField, faceCentres);

    // Set the boundary values
    auto& boundaryPatches = levelSetField.boundaryFieldRef();

    forAll(boundaryPatches, I)
    {
        auto& bPatch = boundaryPatches[I];

        setLevelSet<fvsPatchField<scalar>, surfaceVectorField>(bPatch, faceCentres);
    }
}

void surfaceSet::setLevelSet(pointScalarField& levelSetField) const
{
    const auto& meshVertices = levelSetField.mesh()().points();

    setLevelSet<pointScalarField, pointField>(levelSetField, meshVertices);

    // TODO: Does the code above really set the signed distance for all
    // vertices of the fvMesh including boundary points? (TT)
}



// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //
/*
void surfaceSet::operator+=(const surfaceSet& rhs)
{
}

void surfaceSet::operator/=(const surfaceSet&)
{
}

void surfaceSet::operator-=(const surfaceSet&)
{
}

void surfaceSet::operator%=(const surfaceSet&)
{
}
*/

void surfaceSet::operator+=(const analyticalSurface& rhs)
{
    //operator+=(surfaceSet{rhs});
    analyticalSurface* rawPtr = const_cast<analyticalSurface*>(&rhs);
    orderedSurfaces_.push_back(rawPtr);
    orderedOperations_.push_back(setOperation::unite);
}

void surfaceSet::operator/=(const analyticalSurface& rhs)
{
    //operator/=(surfaceSet{rhs});
    analyticalSurface* rawPtr = const_cast<analyticalSurface*>(&rhs);
    orderedSurfaces_.push_back(rawPtr);
    orderedOperations_.push_back(setOperation::intersect);
}

void surfaceSet::operator-=(const analyticalSurface& rhs)
{
    //operator-=(surfaceSet{rhs});
    analyticalSurface* rawPtr = const_cast<analyticalSurface*>(&rhs);
    orderedSurfaces_.push_back(rawPtr);
    orderedOperations_.push_back(setOperation::complement);
}

void surfaceSet::operator%=(const analyticalSurface& rhs)
{
    //operator%=(surfaceSet{rhs});
    analyticalSurface* rawPtr = const_cast<analyticalSurface*>(&rhs);
    orderedSurfaces_.push_back(rawPtr);
    orderedOperations_.push_back(setOperation::symm_diff);
}

// * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * //
surfaceSet operator+(const surfaceSet& lhs, const analyticalSurface& rhs)
{
    surfaceSet tmp{lhs};
    tmp += rhs;
    return tmp;
}

surfaceSet operator/(const surfaceSet& lhs, const analyticalSurface& rhs)
{
    surfaceSet tmp{lhs};
    tmp /= rhs;
    return tmp;
}

surfaceSet operator-(const surfaceSet& lhs, const analyticalSurface& rhs)
{
    surfaceSet tmp{lhs};
    tmp -= rhs;
    return tmp;
}

surfaceSet operator%(const surfaceSet& lhs, const analyticalSurface& rhs)
{
    surfaceSet tmp{lhs};
    tmp %= rhs;
    return tmp;
}

surfaceSet operator+(const analyticalSurface& lhs, const analyticalSurface& rhs)
{
    return surfaceSet{lhs} + rhs;
}

surfaceSet operator/(const analyticalSurface& lhs, const analyticalSurface& rhs)
{
    return surfaceSet{lhs} / rhs;
}

surfaceSet operator-(const analyticalSurface& lhs, const analyticalSurface& rhs)
{
    return surfaceSet{lhs} - rhs;
}

surfaceSet operator%(const analyticalSurface& lhs, const analyticalSurface& rhs)
{
    return surfaceSet{lhs} % rhs;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
