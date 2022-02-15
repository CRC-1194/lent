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
    Foam::frontTriangleCurvatureModel

SourceFiles
    frontTriangleCurvatureModel.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description

    Curvature model based on the surface tension model described in the
    2012 paper of Tukovic and Jasak

\*---------------------------------------------------------------------------*/

#ifndef frontTriangleCurvatureModel_H
#define frontTriangleCurvatureModel_H

#include "frontBasedCurvatureModel.H"
#include "frontVertexNormalCalculator.H"
#include "lentCommunication.H"
#include "triSurfaceFrontFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class frontTriangleCurvatureModel Declaration
\*---------------------------------------------------------------------------*/

class frontTriangleCurvatureModel
:
    public frontBasedCurvatureModel
{
    // Private data
    tmp<frontVertexNormalCalculator> normalCalculatorTmp_;

    // Private Member Functions
    virtual void computeCurvature(const fvMesh&, const triSurfaceFront&) const;

public:

    TypeName ("frontTriangle");
    
    // Constructors
    frontTriangleCurvatureModel(const dictionary& configDict);
        
    //- Destructor
    virtual ~frontTriangleCurvatureModel() = default;


    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
