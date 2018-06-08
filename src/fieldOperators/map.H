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
    Geometric field map operations.

Author 
    Tomislav Maric

\*---------------------------------------------------------------------------*/

#ifndef map_H
#define map_H

#include "lambda.H" 

namespace Foam {


template
<
    typename Type, 
    template<typename> class PatchField, 
    typename Mesh, 
    typename FilterFunction
>
tmp<GeometricField<Type, PatchField,Mesh>>  
map(
    GeometricField<Type, PatchField, Mesh> const& vf, 
    Type value, 
    FilterFunction filterFunction = true_ref_lambda 
)
{
    tmp<GeometricField<Type, PatchField, Mesh> > resultTmp(
        new GeometricField<Type, PatchField, Mesh>(
            IOobject(
                "map(" + vf.name() + ")", 
                vf.time().timeName(), 
                vf.mesh(), 
                IOobject::NO_READ, 
                IOobject::AUTO_WRITE
            ),
            vf.mesh(), 
            dimensioned<Type>("zero", dimless, pTraits<Type>::zero)
        )
    ); 

    GeometricField<Type, PatchField, Mesh>& result = resultTmp.ref(); 

    forAll(vf, cellI)
    {
        if (filterFunction(vf[cellI]))
        {
            result[cellI] = value; 
        }
    }

    // TODO: Map boundary fields. TM

    return resultTmp;
}; 

template
<
    typename Type, 
    template<typename> class PatchField, 
    typename Mesh, 
    typename FilterFunction
>
tmp<GeometricField<Type, PatchField, Mesh>>  
map(
    tmp<GeometricField<Type, PatchField, Mesh> > vfTmp, 
    Type value, 
    FilterFunction filterFunction
)
{
    return map(vfTmp.ref(), value, filterFunction); 
}

} // End namespace Foam 

// ************************************************************************* //

#endif 
