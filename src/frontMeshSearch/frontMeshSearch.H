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
    Foam::frontMeshSearch

SourceFiles
    frontMeshSearch.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Abstract base class for supplemental mesh search operations.
    fvMesh supports:

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
    );

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
            scalar tolerance=0
        ) const;

        virtual label cellContainingPoint(
            const point& p,
            const fvMesh& mesh,
            const label seedCell = -1
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
