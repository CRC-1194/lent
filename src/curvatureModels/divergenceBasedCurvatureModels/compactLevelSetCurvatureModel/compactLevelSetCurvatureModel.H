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
    Foam::compactLevelSetCurvatureModel

Description

SourceFiles
    compactLevelSetCurvatureModelI.H
    compactLevelSetCurvatureModel.C
    compactLevelSetCurvatureModelIO.C

Authors:
    Tobias Tolle   (tolle@mma.tu-darmstadt.de)
 
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

#ifndef compactLevelSetCurvatureModel_H
#define compactLevelSetCurvatureModel_H

#include "divergenceBasedCurvatureModel.H"

#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class compactLevelSetCurvatureModel Declaration
\*---------------------------------------------------------------------------*/

class compactLevelSetCurvatureModel
:
    public divergenceBasedCurvatureModel
{

private:

    // Private data
    word distanceCorrection_;


    // Private Member Functions
    virtual void computeCurvature(const fvMesh&, const triSurfaceFront&) const;

    void applySphereCorrection(
            const fvMesh& mesh,
            const triSurfaceFront& front
    ) const;

    void applyInterpolationCorrection(
            const fvMesh& mesh,
            const triSurfaceFront& front
    ) const;

    std::vector<label> findNeighbourCells(const fvMesh& mesh, const label& cellLabel) const;



public:

    TypeName ("compactDivGrad");

    // Constructors
    compactLevelSetCurvatureModel(const dictionary& configDict);


    //- Destructor
    virtual ~compactLevelSetCurvatureModel() = default;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
