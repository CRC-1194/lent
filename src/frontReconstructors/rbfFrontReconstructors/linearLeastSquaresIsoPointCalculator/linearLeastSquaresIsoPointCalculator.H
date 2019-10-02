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
    Foam::linearLeastSquaresIsoPointCalculator

Description
    Calculates points on an iso-surface \Gamma := {x, f(x) = s}, where 
    s is the isovalue using an RBF interpolation. 

Authors
    Tomislav Maric, maric@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis Group, TU Darmstadt

SourceFiles
    linearLeastSquaresIsoPointCalculatorI.H
    linearLeastSquaresIsoPointCalculator.C

\*---------------------------------------------------------------------------*/

#ifndef linearLeastSquaresIsoPointCalculator_H
#define linearLeastSquaresIsoPointCalculator_H

#include <memory>
#include "pointFields.H"
#include "volFields.H"
#include <vector>
#include <Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam { namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class linearLeastSquaresIsoPointCalculator Declaration
\*---------------------------------------------------------------------------*/

class linearLeastSquaresIsoPointCalculator 
{

public:

    // Sticking with OpenFOAM datastructures for now. TM. 
    using Points = std::vector<point>;  
    using Labels = std::vector<label>;
    using MatrixX4d = Eigen::Matrix<double, Eigen::Dynamic, 4>;
    using VectorXd = Eigen::VectorXd;
    using DiagonalMatrixXd = Eigen::DiagonalMatrix<double, Eigen::Dynamic>;

    static constexpr double MAX = std::numeric_limits<double>::max(); 

    // Constructor
        linearLeastSquaresIsoPointCalculator(const bool weighted);

    // Member Functions

        const Points& edgePoints() const
        {
            return edgePoints_; 
        }

        const Labels& edgeLabels() const
        {
            return edgeLabels_; 
        }

        const Points& contourPoints() const
        {
            return contourPoints_; 
        }

        const Labels& pointCellLabels() const
        {
            return pointCellLabels_; 
        }

        scalar linearRoot(scalar fiso, scalar f0, scalar f1) const
        {
            return (fiso - f0) / (f1 - f0); 
        }

        word calculationType() const
        {
            if (weighted_)
            {
                return "leastSquaresWeighted";
            }
            else
            {
                return "leastSquaresNoWeighting";
            }
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

        void calcContourPoints 
        (
            const volScalarField& cellPhi, 
            const pointScalarField& pointPhi, 
            const scalar isoValue = 0.
        );

protected:

    Points edgePoints_; 
    Labels edgeLabels_; 

    Points contourPoints_; 
    Labels pointCellLabels_;

    const bool weighted_;

private:
    point projectToLeastSquaresPlane
    (
        point aPoint,
        label cellI,
        const volScalarField& cellPhi,
        const pointScalarField& pointPhi
    ) const;

    DiagonalMatrixXd computeWeights
    (
        const point aPoint,
        const label cellI,
        const volScalarField& cellPhi,
        const pointScalarField& pointPhi
    ) const;

    inline scalar weight(const scalar d, const scalar support) const;
    inline scalar weightInput
    (
        const point p,
        const point& centre,
        const scalar R
    ) const; 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}} // End amespace Foam::FrontTracking  

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
