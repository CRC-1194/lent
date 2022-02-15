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
    Foam::frontMeshSearch

SourceFiles
    frontMeshSearch.C

Authors:
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    Abstract base class for supplemental mesh search operations.

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


#ifndef frontMeshSearch_H
#define frontMeshSearch_H

#include "fvMesh.H"
#include "cellSet.H"
#include "triSurfaceFront.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class frontMeshSearch Declaration
\*---------------------------------------------------------------------------*/

class frontMeshSearch
    :
        public refCount
{
    mutable scalar lastDistance_;
    //mutable cellSet visualizationCellSet_;
    //mutable label iterationCount_;
    //void appendLabelAndWriteCellSet(label cellLabel) const;

public:

    TypeName ("frontMeshSearch");

    declareRunTimeSelectionTable(
        tmp,
        frontMeshSearch,
        Dictionary,
        (
            const dictionary& configDict
        ),
        (configDict)
    )

    // Constructors

        //frontMeshSearch(const Time& runTime);

        frontMeshSearch();

        explicit frontMeshSearch(const dictionary& configDict);

    // Selectors

        static tmp<frontMeshSearch> New(const dictionary& configDict);

    // Destructor
        virtual ~frontMeshSearch();

    // Member Functions

        bool pointIsInCell(
            const point p,
            const label cellLabel,
            const fvMesh& mesh,
            const scalar tolerance=1e-14
        ) const;

        virtual label cellContainingPoint(
            const point& p,
            const fvMesh& mesh,
            const label seedCell = -1,
            const scalar tolerance=1e-14
        ) const;

        labelList pointCellStencil(
            label cellLabel,
            const fvMesh& mesh
        ) const;

    // Member Operators

        void operator=(const frontMeshSearch& rhs);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
