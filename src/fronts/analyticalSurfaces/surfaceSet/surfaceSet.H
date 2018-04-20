/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 AUTHOR,AFFILIATION
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::surfaceSet

Description

SourceFiles
    surfaceSetI.H
    surfaceSet.C
    surfaceSetIO.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceSet_H
#define surfaceSet_H

#include "analyticalSurface.H"

#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class surfaceSet Declaration
\*---------------------------------------------------------------------------*/
enum class setOperation{unite, intersect, complement, symm_diff};

class surfaceSet
{
    // Private data
    std::vector<analyticalSurface> orderedSurfaces_;
    std::vector<setOperation> orderedOperations_;
    mutable std::vector<scalar> signedDistanceBuffer_{};


    // Private Member Functions
    void resizeDistanceBuffer() const;
    void computeDistanceBuffer(const point& P) const;
    scalar computeLevelSetSign() const;
    scalar computelevelSetValue() const;

    bool isOutsidePoint(const label& operationID) const;
    bool pointIsInUnion(const scalar& distA, const scalar& distB) const;
    bool pointIsInIntersection(const scalar& distA, const scalar& distB) const;
    bool pointIsInComplement(const scalar& distA, const scalar& distB) const;
    bool pointIsInSymmDiff(const scalar& distA, const scalar& distB) const;


public:

    // Constructors
    surfaceSet() = default;
    surfaceSet(const analyticalSurface& surface);
    surfaceSet(const surfaceSet& rhs);


    //- Destructor
    ~surfaceSet() = default;


    // Member Functions
    scalar levelSetValue(const point& P) const;
    
    template<class Container, class PointField>
    void setLevelSet(Container& levelSetField, const PointField& points) const;

    void setLevelSet(volScalarField& levelSetField) const;
    void setLevelSet(surfaceScalarField& levelSetField) const;
    void setLevelSet(pointScalarField& levelSetField) const;


    // Member Operators
    surfaceSet& operator=(const surfaceSet&) = default;
    surfaceSet& operator=(surfaceSet&&) = default;

    // TODO: the following operators need a graph structure
    // rather than a sequential one to preserve the correct order
    // of operation execution
    /*
    void operator+=(const surfaceSet& rhs);
    void operator/=(const surfaceSet& rhs);
    void operator-=(const surfaceSet& rhs);
    void operator%=(const surfaceSet& rhs);
    */

    void operator+=(const analyticalSurface& rhs);
    void operator/=(const analyticalSurface& rhs);
    void operator-=(const analyticalSurface& rhs);
    void operator%=(const analyticalSurface& rhs);


    // Friend Functions

    // Friend Operators
    // TODO: depending on the implementation of the compound assignment
    // operators, these operators do not need access to the internals of the
    // surface sets and thus can be defined as non-friend functions
    /*
    friend surfaceSet operator+(const surfaceSet& lhs, const surfaceSet& rhs);
    friend surfaceSet operator/(const surfaceSet& lhs, const surfaceSet& rhs);
    friend surfaceSet operator-(const surfaceSet& lhs, const surfaceSet& rhs);
    friend surfaceSet operator%(const surfaceSet& lhs, const surfaceSet& rhs);
    */
};

surfaceSet operator+(const surfaceSet& lhs, const analyticalSurface& rhs);
surfaceSet operator/(const surfaceSet& lhs, const analyticalSurface& rhs);
surfaceSet operator-(const surfaceSet& lhs, const analyticalSurface& rhs);
surfaceSet operator%(const surfaceSet& lhs, const analyticalSurface& rhs);

surfaceSet operator+(const analyticalSurface& lhs, const analyticalSurface& rhs);
surfaceSet operator/(const analyticalSurface& lhs, const analyticalSurface& rhs);
surfaceSet operator-(const analyticalSurface& lhs, const analyticalSurface& rhs);
surfaceSet operator%(const analyticalSurface& lhs, const analyticalSurface& rhs);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
