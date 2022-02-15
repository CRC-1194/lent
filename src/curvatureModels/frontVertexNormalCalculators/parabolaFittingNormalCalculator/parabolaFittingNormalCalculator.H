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
    Foam::parabolaFittingNormalCalculator

SourceFiles
    parabolaFittingNormalCalculator.C

Authors:
    Tobias Tolle   (tolle@mma.tu-darmstadt.de)

Description

    Compute the front vertex normals from a parabola fitted at each vertex.
    The parameters of the parabola f(x,y)=ax^2 + by^2 + cxy + dx + ey
    are found by solving a least-squares problem considering all vertices
    of the two-ring neighbourhood.

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

#ifndef parabolaFittingNormalCalculator_H
#define parabolaFittingNormalCalculator_H

#include "quaternion.H"

#include "frontVertexNormalCalculator.H"

#include <Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    using matrixXd = Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using vectorXd = Eigen::Matrix<scalar, Eigen::Dynamic, 1>;

/*---------------------------------------------------------------------------*\
                         Class parabolaFittingNormalCalculator Declaration
\*---------------------------------------------------------------------------*/

class parabolaFittingNormalCalculator
:
    public frontVertexNormalCalculator
{
private:
    quaternion rotationQuaternion(const triSurfaceFront& front, const label& pointLabel) const;
    vectorXd fitParabola(const std::vector<Foam::vector>& points) const;
    vector computeNormal(const vectorXd& coefficients, const quaternion& rotation) const;

public:

    TypeName ("parabolaFitting");

    // Constructors
    parabolaFittingNormalCalculator() = default;
    explicit parabolaFittingNormalCalculator(const dictionary& configDict);

    //- Destructor
    virtual ~parabolaFittingNormalCalculator() = default;

    // Member Functions
    virtual tmp<triSurfaceFrontPointVectorField> vertexNormals(const fvMesh&, const triSurfaceFront&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
