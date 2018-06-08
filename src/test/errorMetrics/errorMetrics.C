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
    Foam::errorMetrics

SourceFiles
    errorMetrics.C

Author
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    Computes various scalar error metrics from a given error set, e.g. the
    difference between a computed and an exact velocity field.
    
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

#include "errorMetrics.H"

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
scalar errorMetrics::powerMeanError(scalar x) const
{
    scalar result = 0.0;

    forAll(errorSet_, I)
    {
        result += std::pow(errorSet_[I], x);
    }

    result /= errorSet_.size();

    return std::pow(result, 1.0/x);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
errorMetrics::errorMetrics(const List<scalar>& errorSet)
:
    errorSet_(errorSet)
{
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
    return max(errorSet_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
