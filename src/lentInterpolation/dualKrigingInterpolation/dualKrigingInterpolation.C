/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "dualKrigingInterpolation.H"

#include <algorithm>

namespace Foam {
namespace FrontTracking {


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void dualKrigingInterpolation::removeDuplicates(std::vector<label>& listOfLabels) const
{
    std::sort(listOfLabels.begin(), listOfLabels.end());
    std::vector<label>::iterator newEnd;
    newEnd = std::unique(listOfLabels.begin(), listOfLabels.end());
    listOfLabels.resize(std::distance(listOfLabels.begin(), newEnd));
}

std::vector<label> dualKrigingInterpolation::cellNeighbourhood(const label& cellLabel, const fvMesh& mesh) const
{
    std::vector<label> neighbourCellLabels{};

    const auto& cellPoints = mesh.cellPoints()[cellLabel];

    for (const auto& pointLabel : cellPoints)
    {
        const auto& connectedCells = mesh.pointCells()[pointLabel];

        for (const auto& aCellLabel : connectedCells)
        {
            neighbourCellLabels.push_back(aCellLabel);
        }
    }

    removeDuplicates(neighbourCellLabels);

    return neighbourCellLabels;
}

void dualKrigingInterpolation::cacheCellCentres(const std::vector<label>& cellLabels, const fvMesh& mesh)
{
    const auto& cellCentres = mesh.C();

    for (unsigned int index = 0; index < cellLabels.size(); ++index)
    {
        cachedCellCentres_[index] = cellCentres[cellLabels[index]];
    }
}

void dualKrigingInterpolation::adaptSystemSize(const unsigned int nPoints)
{
    // The 4 arises from using a linear drift
    unsigned int nDriftCoeffs = 4;
    
    // The 10 arises from using a full quadratic drift
    //unsigned int nDriftCoeffs = 10;
    
    auto newDim = nPoints + nDriftCoeffs;

    cachedCellCentres_.resize(nPoints);

    dualKrigingSystem_.resize(newDim, newDim);
    rhs_.resize(newDim);

    driftCoefficients_.resize(nDriftCoeffs);
    fluctuationCoefficients_.resize(nPoints);
}

void dualKrigingInterpolation::setupKrigingSystem()
{
    // Reassignment is required for setting the system to zero as
    // multiplication fails for NaN values arising from resizing
    dualKrigingSystem_ = Eigen::MatrixXd::Zero(dualKrigingSystem_.rows(), dualKrigingSystem_.cols());

    // TODO: Optimization: resulting matrix symmetric. This fact is 
    // not exploited here
    auto nPoints = cachedCellCentres_.size();

    // Setup covariance block
    for (unsigned int I = 0; I < nPoints; ++I)
    {
        for (unsigned int K = 0; K < nPoints; ++K)
        {
            dualKrigingSystem_(I,K) = coVar(cachedCellCentres_[I], cachedCellCentres_[K]);
        }
    }

    // Setup one columns
    for (unsigned int I = 0; I < nPoints; ++I)
    {
        dualKrigingSystem_(nPoints, I) = 1.0;
        dualKrigingSystem_(I, nPoints) = 1.0;
    }

    // Setup drift blocks using the point coordinates
    for (unsigned int I = 0; I < 3; ++I)
    {
        for (unsigned int K = 0; K < nPoints; ++K)
        {
            // linear coefficients x, y, z
            dualKrigingSystem_(nPoints+1+I, K) = cachedCellCentres_[K][I];
            dualKrigingSystem_(K, nPoints+1+I) = cachedCellCentres_[K][I];

            /*
            // quadratic coefficients x^2, y^2, z^2
            dualKrigingSystem_(nPoints+7+I, K) = cachedCellCentres_[K][I]
                                                    * cachedCellCentres_[K][I];
            dualKrigingSystem_(K, nPoints+7+I) = cachedCellCentres_[K][I]
                                                    * cachedCellCentres_[K][I];
            
            // bilinear coefficients xy, yz, zx
            auto B = I+1;

            if (I == 2)
            {
                B = 0;
            }

            dualKrigingSystem_(nPoints+4+I, K) = cachedCellCentres_[K][I]
                                                    * cachedCellCentres_[K][B];
            dualKrigingSystem_(K, nPoints+4+I) = cachedCellCentres_[K][I]
                                                    * cachedCellCentres_[K][B];
                                                    * */
        }
    }
}

void dualKrigingInterpolation::setupRightHandSide(const std::vector<label>& cellLabels, const volVectorField& phi, const label index)
{
    rhs_ = Eigen::VectorXd::Zero(rhs_.rows());

    for (unsigned int I = 0; I < cellLabels.size(); ++I)
    {
        rhs_(I) = phi[cellLabels[I]][index];
    }
}

scalar dualKrigingInterpolation::coVar(const point& A, const point& B) const
{
    // TODO: I assume this will not work as intended...
    // linear covariance
    return mag(A - B);
    // Quadratic covariance
    //return magSqr(A - B);
}

vector dualKrigingInterpolation::drift(const point& vertex) const
{
    // First entry of coefficients is an absolute value
    vector driftValue{driftCoefficients_[0]};

    for (unsigned int I = 0; I < 3; ++I)
    {
        // linear terms
        driftValue += driftCoefficients_[I+1]*vertex[I];

        /*
        // quadratic terms
        driftValue += driftCoefficients_[I+7]*vertex[I]*vertex[I];

        // bilinear terms
        auto B = I+1;

        if (I == 2)
        {
            B = 0;
        }

        driftValue += driftCoefficients_[I+4]*vertex[I]*vertex[B];
        */
    }

    return driftValue;
}

vector dualKrigingInterpolation::fluctuation(const point& vertex) const
{
    vector fluctuationValue{0,0,0};

    for (unsigned int I = 0; I < cachedCellCentres_.size(); ++I)
    {
        fluctuationValue += fluctuationCoefficients_[I]*coVar(vertex, cachedCellCentres_[I]);
    }

    return fluctuationValue;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void dualKrigingInterpolation::computeDualKrigingParameters
(
    const label& cellLabel,
    const fvMesh& mesh,
    const volVectorField& phi
)
{
    auto cellLabels = cellNeighbourhood(cellLabel, mesh);
    auto nPoints = cellLabels.size();

    adaptSystemSize(nPoints);
    cacheCellCentres(cellLabels, mesh);

    setupKrigingSystem();

    // Re-use the decomposition. No need to compute it multiple times
    auto decomposedSytem = dualKrigingSystem_.colPivHouseholderQr();

    for (unsigned int I = 0; I < 3; ++I)
    {
        setupRightHandSide(cellLabels, phi, I);

        Eigen::VectorXd coeffs = decomposedSytem.solve(rhs_);

        for (unsigned int K = 0; K < nPoints; ++K)
        {
            fluctuationCoefficients_[K][I] = coeffs(K);
        }

        for (unsigned int K = nPoints; K < coeffs.rows(); ++K)
        {
            driftCoefficients_[K-nPoints][I] = coeffs(K);
        }
    }
}

vector dualKrigingInterpolation::interpolateTo(const point& vertex) const
{
    return drift(vertex) + fluctuation(vertex);
}

void dualKrigingInterpolation::selfTest(const volVectorField& phi, const label& cellLabel) const
{
    auto cellIDs = cellNeighbourhood(cellLabel, phi.mesh());
    const auto& C = phi.mesh().C();

    for (const auto& cellID : cellIDs)
    {
        auto est = interpolateTo(C[cellID]);

        if (mag(est - phi[cellID]) > SMALL)
        {
            Info << "Value = " << phi[cellID] << "; interpolated = "
                 << est << '\n';
        }
    } 
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
