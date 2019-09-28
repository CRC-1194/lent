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
    Foam::rbfIsoPointCalculator

Description
    Calculates points on an iso-surface \Gamma := {x, f(x) = s}, where 
    s is the isovalue using an RBF interpolation. 

Authors
    Tomislav Maric, maric@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis Group, TU Darmstadt

SourceFiles
    rbfIsoPointCalculatorI.H
    rbfIsoPointCalculator.C

\*---------------------------------------------------------------------------*/

#ifndef rbfIsoPointCalculator_H
#define rbfIsoPointCalculator_H

#include "isoPointCalculator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam { namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class rbfIsoPointCalculator Declaration
\*---------------------------------------------------------------------------*/

class rbfIsoPointCalculator 
    : 
        public isoPointCalculator
{

    using Scalars = std::vector<double>;
    using Vectors = std::vector<vector>; 

    //Scalars contourPointRbfVals_;  
    //Vectors contourPointRbfGrads_; 

public:

    using isoPointCalculator::MAX;

    // Member Functions
    
        template<typename MeshRbfs> 
        void calcContourPoints 
        (
            const volScalarField& cellPhi, 
            const pointScalarField& pointPhi, 
            MeshRbfs const& rbfs,
            const scalar isoValue = 0.
        );

        template<typename MeshRbfs> 
        void rbfCorrectContourPoint(
            point& contourPoint, 
            label cellI,
            MeshRbfs const& rbf
        );

        template<typename MeshRbfs> 
        void rbfCorrectContourPoints(MeshRbfs const& rbf);

        //Scalars const& contourPointRbfValues() const
        //{
            //return contourPointRbfVals_; 
        //}

        //Vectors const& contourPointRbfGradients() const
        //{
            //return contourPointRbfGrads_; 
        //}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}} // End amespace Foam::FrontTracking  

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rbfIsoPointCalculator.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
