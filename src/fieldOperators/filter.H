/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

Description
    Field filters. 

Author
    Tomislav Maric 

\*---------------------------------------------------------------------------*/

#ifndef filter_H
#define filter_H

#include "DynamicField.H"
#include "dimensionedType.H"

namespace Foam {

template
<
    typename Type, 
    template<typename> class PatchField, 
    typename Mesh, 
    typename FilterFunction
>
DynamicField<Type>
filter(
    GeometricField<Type, PatchField, Mesh> const& vf, 
    FilterFunction filterFunction = true_ref_lambda 
)
{
    DynamicField<Type> result;

    forAll(vf, cellI)
    {
        if (filterFunction(vf[cellI]))
        {
            result.append(vf[cellI]); 
        }
    }

    return result;
}; 

template
<
    typename Type, 
    template<typename> class PatchField, 
    typename Mesh, 
    typename FilterFunction
>
DynamicField<Type>
filter(
    tmp<GeometricField<Type, PatchField, Mesh> > vfTmp,
    FilterFunction filterFunction 
)
{
    return filter(vfTmp(), filterFunction); 
}

} // End namespace Foam 

// ************************************************************************* //

#endif 
