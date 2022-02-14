/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "error.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::interpolationCellPointTaylor<Type>::interpolate
(
    const vector& position,
    const label cellI,
    const label
) const
{
    const auto& gradPsi_= psiGrad_();  
    const auto& mesh = gradPsi_.mesh(); 
    const auto& C = mesh.C(); 

    return (this->psi_[cellI] + ((position - C[cellI]) & gradPsi_[cellI])); 
}


template<class Type>
inline Type Foam::interpolationCellPointTaylor<Type>::interpolate
(
    const vector& position,
    const tetIndices& tetIs,
    const label faceI
) const
{
    notImplemented("interpolationCellPointTaylor::interpolate(const vector&, const tetIndices& const label)");
    return Type(); 
}


// ************************************************************************* //
