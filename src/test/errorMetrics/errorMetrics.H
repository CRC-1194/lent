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

#ifndef errorMetrics_H
#define errorMetrics_H

#include <cmath>
#include <map>
#include <vector>

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class errorMetrics Declaration
\*---------------------------------------------------------------------------*/

// TODO: make this class a template (TT)
class errorMetrics
{
    // Private data
    std::vector<scalar> errorSet_;

    // Private Member Functions
    scalar powerMeanError(scalar x) const;


public:

    // Constructors
    errorMetrics(const List<scalar>& errorSet);
    errorMetrics(const std::vector<scalar>& errorSet);


    // Member Functions
    scalar arithmeticMeanError() const;
    scalar quadraticMeanError() const;
    scalar maximumError() const;
    scalar medianError() const;
    scalar standardDeviation() const;
    std::map<scalar, label> errorDistribution(const label& resolution) const;
    std::map<scalar, scalar> errorDistributionNormalized(const label& resolution) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
