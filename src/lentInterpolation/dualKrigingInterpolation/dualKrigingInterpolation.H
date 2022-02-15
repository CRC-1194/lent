/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | 
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

Class
    Foam::dualKrigingInterpolation

SourceFiles
    dualKrigingInterpolationI.H
    dualKrigingInterpolation.C
    dualKrigingInterpolationIO.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

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

#ifndef dualKrigingInterpolation_H
#define dualKrigingInterpolation_H

#include "fvMesh.H"
#include "fvCFD.H"

#include <vector>
#include <Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class dualKrigingInterpolation Declaration
\*---------------------------------------------------------------------------*/

class dualKrigingInterpolation
{
    // Private data
    std::vector<vector> cachedCellCentres_;

    Eigen::MatrixXd dualKrigingSystem_;

    Eigen::VectorXd rhs_;

    std::vector<vector> driftCoefficients_;
    std::vector<vector> fluctuationCoefficients_;


    // Private Member Functions
    void removeDuplicates(std::vector<label>& listOfLabels) const;
    std::vector<label> cellNeighbourhood(const label& cellLabel, const fvMesh& mesh) const;
    void cacheCellCentres(const std::vector<label>& cellLabels, const fvMesh& mesh);
    void adaptSystemSize(const label nPoints);

    void setupKrigingSystem();
    void setupRightHandSide(const std::vector<label>& cellLabels, const volVectorField& phi, const direction index);

    scalar coVar(const point& A, const point& B) const;
    vector drift(const point& vertex) const;
    vector fluctuation(const point& vertex) const;

public:


    // Constructors
    dualKrigingInterpolation() = default;


    //- Destructor
    ~dualKrigingInterpolation() = default;


    // Member Functions
    void computeDualKrigingParameters
         (
            const label& cellLabel,
            const fvMesh& mesh,
            const volVectorField& phi
         );
    vector interpolateTo(const point& vertex) const;

    void selfTest(const volVectorField& phi, const label& cellLabel) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
