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
    Foam::parabolaFitting2DNormalCalculator

SourceFiles
    parabolaFitting2DNormalCalculator.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description

    Compute the front vertex normals from a parabola fitted at each vertex.
    The parameters of the parabola f(x,y)=ax^2 + by^2 + cxy + dx + ey
    are found by solving a least-squares problem considering all vertices
    of the two-ring neighbourhood.

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

#ifndef parabolaFitting2DNormalCalculator_H
#define parabolaFitting2DNormalCalculator_H

#include "quaternion.H"

#include "frontVertexNormalCalculator.H"

#include <Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    using matrixXd = Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using vectorXd = Eigen::Matrix<scalar, Eigen::Dynamic, 1>;

/*---------------------------------------------------------------------------*\
                         Class parabolaFitting2DNormalCalculator Declaration
\*---------------------------------------------------------------------------*/

class parabolaFitting2DNormalCalculator
:
    public frontVertexNormalCalculator
{
private:
    mutable tensor projector_;
    mutable quaternion globalToReference_;
    mutable tmp<triSurfaceFrontPointVectorField> vertexNormalTmp_;

    void initializeVertexNormal(const fvMesh&, const triSurfaceFront&) const;

    void computeProjectionData(const fvMesh& mesh) const;
    quaternion rotationQuaternion(const triSurfaceFront& front, const label& pointLabel) const;
    vectorXd fitParabola(const std::vector<Foam::vector>& points) const;
    vector computeNormal(const vectorXd& coefficients, const quaternion& rotation) const;

public:

    TypeName ("parabolaFitting2D");

    // Constructors
    parabolaFitting2DNormalCalculator() = default;
    explicit parabolaFitting2DNormalCalculator(const dictionary& configDict);

    //- Destructor
    virtual ~parabolaFitting2DNormalCalculator() = default;

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
