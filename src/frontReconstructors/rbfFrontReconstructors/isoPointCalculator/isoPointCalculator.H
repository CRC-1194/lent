/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Tomislav Maric, 
     \\/     M anipulation  |                    TU Darmstadt
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
    Foam::isoPointCalculator

Description
    Calculates points on an iso-surface \Gamma := {x, f(x) = s}, where 
    s is the isovalue using an RBF interpolation. 

Authors
    Tomislav Maric, maric@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis Group, TU Darmstadt

SourceFiles
    isoPointCalculatorI.H
    isoPointCalculator.C

\*---------------------------------------------------------------------------*/

#ifndef isoPointCalculator_H
#define isoPointCalculator_H

#include <memory>
#include "pointFields.H"
#include "volFields.H"
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam { namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class isoPointCalculator Declaration
\*---------------------------------------------------------------------------*/

class isoPointCalculator 
{

public:

    // Sticking with OpenFOAM datastructures for now. TM. 
    using Points = std::vector<point>;  
    using Labels = std::vector<label>; 

    static constexpr double MAX = std::numeric_limits<double>::max(); 

    // Member Functions

        const Points& edgePoints() const
        {
            return edgePoints_; 
        }

        const Labels& edgeLabels() const
        {
            return edgeLabels_; 
        }

        const Points& cellPoints() const
        {
            return cellPoints_; 
        }

        const Labels& cellLabels() const
        {
            return cellLabels_; 
        }

        scalar linearRoot(scalar fiso, scalar f0, scalar f1) const
        {
            return (fiso - f0) / (f1 - f0); 
        }

        inline bool isoValueInInterval
        (
            const scalar isoValue, 
            const scalar phi0, 
            const scalar phi1
        ) const;

        inline bool equalUnderTolerance(
            const scalar s1, 
            const scalar s2, 
            const scalar tol
        ) const;

        void calcEdgePoints 
        (
            const pointScalarField& pointPhi, 
            const scalar isoValue = 0.
        );

        void calcCellPoints 
        (
            const volScalarField& cellPhi, 
            const pointScalarField& pointPhi, 
            const scalar isoValue = 0.
        );

protected:

    Points edgePoints_; 
    Labels edgeLabels_; 

    Points cellPoints_; 
    Labels cellLabels_; 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}} // End amespace Foam::FrontTracking  

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
