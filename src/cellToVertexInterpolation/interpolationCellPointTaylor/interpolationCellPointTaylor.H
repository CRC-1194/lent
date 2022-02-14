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
    Foam::interpolationCellPointTaylor

Description
    Linear extrapolation from the cell center. 

Authors:
    Tomislav Maric (maric@mma.tu-darmstadt.de)
 
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

#ifndef interpolationCellPointTaylor_H
#define interpolationCellPointTaylor_H

#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class interpolationCellPointTaylor Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolationCellPointTaylor
:
    public interpolation<Type>
{
public: 

    typedef GeometricField
    <
        typename outerProduct<vector,Type>::type, 
        fvPatchField, 
        volMesh
    > GradFieldType; 

private: 

    tmp<GradFieldType> psiGrad_;  

public:

    //- Runtime type information
    TypeName("cellPointTaylor");


    // Constructors

        //- Construct from components
        interpolationCellPointTaylor
        (
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );


    // Member Functions

        //- Interpolate field to the given point in the given cell
        inline Type interpolate
        (
            const vector& position,
            const label cellI,
            const label faceI = -1
        ) const;

        //- Interpolate field to the given point in the tetrahedron
        //  defined by the given indices.
        inline Type interpolate
        (
            const vector& position,
            const tetIndices& tetIs,
            const label faceI = -1
        ) const;

        const GradFieldType& gradFieldRef() const
        {
            return psiGrad_(); 
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interpolationCellPointTaylorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolationCellPointTaylor.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
