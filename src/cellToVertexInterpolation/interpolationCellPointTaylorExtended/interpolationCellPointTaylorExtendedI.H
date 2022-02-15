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
    Foam::interpolationCellPointTaylorExtended

Description
    Linear extrapolation from the cell center using second-order derivatives. 

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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::interpolationCellPointTaylorExtended<Type>::interpolate
(
    const vector& position,
    const label cellI,
    const label
) const
{
    auto result = pTraits<Type>::zero; 

    const auto& gradPsi_= this->gradFieldRef();  
    const auto& mesh = gradPsi_.mesh(); 
    const auto& C = mesh.C(); 
    const auto& cellCells = mesh.cellCells(); 
    const auto& neighborCells = cellCells[cellI];

    forAll (neighborCells, nCellI)
    {
        const auto nCell = neighborCells[nCellI]; 
        const auto nDelta = position - C[nCell]; 
        result += this->psi_[nCell] + (nDelta & gradPsi_[nCell]);  
    }

    // Arithmetic average.
    // TODO: IDW based average? Volume based average? TM.   
    result /= neighborCells.size(); 

    return result; 
}

// ************************************************************************* //
