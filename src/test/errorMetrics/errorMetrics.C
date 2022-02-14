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
    Foam::errorMetrics

SourceFiles
    errorMetrics.C

Author
    Tobias Tolle   tolle@mma.tu-darmstadt.de

Description
    Computes various scalar error metrics from a given error set, e.g. the
    difference between a computed and an exact velocity field.

\*---------------------------------------------------------------------------*/

#include "errorMetrics.H"

#include <algorithm>

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
scalar errorMetrics::powerMeanError(scalar x) const
{
    scalar result = 0.0;

    for (const auto& error : errorSet_)
    {
        result += std::pow(error, x);
    }

    result /= static_cast<scalar>(errorSet_.size());

    return std::pow(result, 1.0/x);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
errorMetrics::errorMetrics(const List<scalar>& errorSet)
:
    errorSet_{}
{
    errorSet_.resize(errorSet.size());

    forAll(errorSet, I)
    {
        errorSet_[I] = errorSet[I];
    }

    std::sort(errorSet_.begin(), errorSet_.end());
}

errorMetrics::errorMetrics(const std::vector<scalar>& errorSet)
:
    errorSet_{errorSet}
{
    std::sort(errorSet_.begin(), errorSet_.end());
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar errorMetrics::arithmeticMeanError() const
{
    return powerMeanError(1.0);
}

scalar errorMetrics::quadraticMeanError() const
{
    return powerMeanError(2.0);
}

scalar errorMetrics::maximumError() const
{
    scalar max = 0.0;

    for (const auto& error : errorSet_)
    {
        if (mag(error) > max)
        {
            max = error;
        }
    }
    return max;
}

scalar errorMetrics::medianError() const
{
    const auto& size = errorSet_.size();

    if (size%2 == 1)
    {
        return errorSet_[(size - 1)/2];
    }
    else
    {
        return (errorSet_[size/2] + errorSet_[size/2 - 1])/2.0;
    }
}

scalar errorMetrics::standardDeviation() const
{
    scalar stdDev = 0.0;

    auto mean = arithmeticMeanError();

    for (const auto& x : errorSet_)
    {
        stdDev += (x - mean)*(x - mean);
    }

    return sqrt(stdDev/static_cast<scalar>(errorSet_.size()));
}

std::map<scalar, label> errorMetrics::errorDistribution( const label& resolution) const
{
    std::map<scalar, label> distribution{};

    // Remember that data set is sorted, so no need to search for min and max
    auto minError = errorSet_[0];
    auto maxError = errorSet_[errorSet_.size() - 1];

    scalar increment = (maxError - minError)/resolution;
    scalar key = minError + 0.5*increment;
    scalar limit = minError + increment;

    distribution[key] = 0;

    for (const auto& error : errorSet_)
    {
        if (error < limit)
        {
            ++distribution[key];
        }
        else
        {
            // Ensure empty intervals are considered correctly
            do
            {
                key += increment;
                limit += increment;
                distribution[key] = 0;

            }
            while (error >= limit);

            ++distribution[key];
        }
    }

    return distribution;
}

std::map<scalar, scalar> errorMetrics::errorDistributionNormalized(const label& resolution) const
{
    std::map<scalar, scalar> normalizedDistribution{};

    auto distribution = errorDistribution(resolution);

    for (const auto& entry : distribution)
    {
        normalizedDistribution[entry.first] = scalar(entry.second)/scalar(errorSet_.size());
    }

    return normalizedDistribution;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
