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

Class
    Foam::noiseGenerator

Description

SourceFiles
    noiseGeneratorI.H
    noiseGenerator.C
    noiseGeneratorIO.C

\*---------------------------------------------------------------------------*/

#ifndef noiseGenerator_H
#define noiseGenerator_H

#include "scalar.H"

#include "randomNumberGenerator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class noiseGenerator Declaration
\*---------------------------------------------------------------------------*/

class noiseGenerator
{
private:
    // Private data
    mutable randomNumberGenerator rng_;


    // Private Member Functions
    scalar scalarNoise(const scalar& maxMagnitude) const;

public:

    // Constructors
    noiseGenerator() = default;


    //- Destructor
    ~noiseGenerator() = default;


    // Member Functions
    // Specialization for tensorType = scalar
    template
    <
        class tensorType,
        typename std::enable_if<!std::is_class<tensorType>::value>::type* = nullptr
    >
    tensorType noise(const scalar& maxMagnitude) const;

    // tensorType = vector, tensor
    template
    <
        class tensorType,
        typename std::enable_if<std::is_class<tensorType>::value>::type* = nullptr
    >
    tensorType noise(const scalar& maxMagnitude) const;

    // Use a noise magnitude prescribed per component
    template
    <
        class tensorType,
        typename std::enable_if<std::is_class<tensorType>::value>::type* = nullptr
    >
    tensorType noise(const tensorType& magnitudes) const;

    // Noise generation for fields
    
    // Uniform noise for each tensor component
    template
    <
        typename tensorType,
        template<typename> class container
    >
    void addNoiseTo(container<tensorType>& field, const scalar magnitude) const;
    
    // Per component noise
    template
    <
        typename tensorType,
        template<typename> class container,
        typename std::enable_if<std::is_class<tensorType>::value>::type* = nullptr
    >
    void addNoiseTo(container<tensorType>& field, const tensorType& magnitudes) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "noiseGeneratorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
