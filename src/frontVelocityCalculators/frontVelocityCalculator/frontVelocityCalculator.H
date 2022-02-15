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
    Foam::frontVelocityCalculator

SourceFiles
    frontVelocityCalculator.C

Authors
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    Interface for the front velocity calculation.

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


#ifndef frontVelocityCalculator_H
#define frontVelocityCalculator_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "tmp.H"
#include "refCount.H"
#include "volFieldsFwd.H"
#include "triSurfaceFront.H"
#include "triSurfaceFields.H"
#include "frontMeshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class frontVelocityCalculator Declaration
\*---------------------------------------------------------------------------*/

class frontVelocityCalculator
    :
        public refCount
{

    tmp<frontMeshSearch> searchAlgTmp_;

public:

    TypeName ("frontVelocityCalculator");

    declareRunTimeSelectionTable (
        tmp,
        frontVelocityCalculator,
        Dictionary,
        (
            const dictionary& configDict
        ),
        (configDict)
    )

    // Constructors
    frontVelocityCalculator(const dictionary& configDict);

    // Selectors
    static tmp<frontVelocityCalculator> New(const dictionary& configDict);

    // Destructor
    virtual ~frontVelocityCalculator() = default;

    // Member Functions
    bool pointIsInCell(
        const point p,
        const label cellLabel,
        const fvMesh& mesh,
        scalar tolerance=-SMALL
    ) const
    {
        return searchAlgTmp_->pointIsInCell(p, cellLabel, mesh, tolerance);
    }

    label cellContainingPoint(
        const point& p,
        const fvMesh& mesh,
        const label seedCell = -1
    ) const
    {
        return searchAlgTmp_->cellContainingPoint(p, mesh, seedCell);
    }

    virtual void calcFrontVelocity(
        triSurfacePointVectorField& frontVelocity,
        const volVectorField& U,
        labelList& elementCells
    ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
