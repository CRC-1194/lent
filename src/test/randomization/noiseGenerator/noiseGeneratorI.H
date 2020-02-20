/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "UList.H"

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template
<
    class tensorType,
    typename std::enable_if<!std::is_class<tensorType>::value>::type* //= nullptr
>
tensorType noiseGenerator::noise(const scalar& maxMagnitude) const
{
    tensorType result = pTraits<tensorType>::zero;

    result = scalarNoise(maxMagnitude);

    return result;
}


template
<
    class tensorType,
    typename std::enable_if<std::is_class<tensorType>::value>::type* //= nullptr
>
tensorType noiseGenerator::noise(const scalar& maxMagnitude) const
{
    tensorType result = pTraits<tensorType>::zero;

    for (direction I = 0; I != pTraits<tensorType>::nComponents; ++I)
    {
        result[I] = scalarNoise(maxMagnitude);
    }

    return result;
}


template
<
    class tensorType,
    typename std::enable_if<std::is_class<tensorType>::value>::type* //= nullptr
>
tensorType noiseGenerator::noise(const tensorType& magnitudes) const
{
    tensorType result = pTraits<tensorType>::zero;

    for (direction I = 0; I != pTraits<tensorType>::nComponents; ++I)
    {
        result[I] = scalarNoise(magnitudes[I]);
    }

    return result;
}


template
<
    typename tensorType,
    template<typename> class container
>
void noiseGenerator::addNoiseTo(container<tensorType>& field, const scalar magnitude) const
{
    for (tensorType& fieldValue : field)
    {
        fieldValue += noise<tensorType>(magnitude);
    }
}

template
<
    typename tensorType,
    template<typename> class container,
    typename std::enable_if<std::is_class<tensorType>::value>::type* //= nullptr
>
void noiseGenerator::addNoiseTo(container<tensorType>& field, const tensorType& magnitudes) const
{
    for (tensorType& fieldValue : field)
    {
        fieldValue += noise<tensorType>(magnitudes);
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
