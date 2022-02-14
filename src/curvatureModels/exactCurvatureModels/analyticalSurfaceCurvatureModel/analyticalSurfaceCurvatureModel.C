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
    Foam::analyticalSurfaceCurvatureModel

SourceFiles
    analyticalSurfaceCurvatureModel.C

Description

    Curvature model that uses an analytical surface description to
    provide the curvature.
    Expects a subdict named "frontSurface" with the parameters for an
    analyticalSurface.

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

#include "analyticalSurfaceCurvatureModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalSurfaceCurvatureModel, 0);
    addToRunTimeSelectionTable(curvatureModel, analyticalSurfaceCurvatureModel, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

analyticalSurfaceCurvatureModel::analyticalSurfaceCurvatureModel(const dictionary& configDict)
:
    exactCurvatureModel{configDict},
    surfaceTmp_{analyticalSurface::New(configDict.subDict("frontSurface"))}
{}

analyticalSurfaceCurvatureModel::analyticalSurfaceCurvatureModel(const dictionary& configDict, const analyticalSurface& surface)
:
    exactCurvatureModel{configDict},
    surfaceTmp_{surface}
{}

// * * * * * * * * * * * * * * * * Public member functions * * * * * * * * * //
scalar analyticalSurfaceCurvatureModel::curvatureAtPoint(const point& P) const
{
    return surfaceTmp_->curvatureAt(P);
}

vector analyticalSurfaceCurvatureModel::normalAtPoint(const point& P) const
{
    return surfaceTmp_->normalToPoint(P);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
