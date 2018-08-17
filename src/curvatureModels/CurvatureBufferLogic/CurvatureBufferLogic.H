/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 AUTHOR,AFFILIATION
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
    Foam::CurvatureBufferLogic

Description

SourceFiles
    CurvatureBufferLogicI.H
    CurvatureBufferLogic.C
    CurvatureBufferLogicIO.C

\*---------------------------------------------------------------------------*/

#ifndef CurvatureBufferLogic_H
#define CurvatureBufferLogic_H

#include "fvMesh.H"
#include "pTraits.H"

#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class CurvatureBufferLogic Declaration
\*---------------------------------------------------------------------------*/
template<class fieldType, class meshType, class tensorType>
class CurvatureBufferLogic
{

private:

    // Private data
    word bufferFieldName_{"curvature"};
    mutable bool bufferInitialized_{false};
    mutable label lastTimeUpdated_{-1};
    mutable std::shared_ptr<fieldType> curvatureBufferPtr_{nullptr};

    // Private Member Functions
    void initializeBufferField(const meshType&) const;


public:

    // Constructors
    CurvatureBufferLogic(){};
    CurvatureBufferLogic(const word& bufferFieldName);


    //- Destructor
    ~CurvatureBufferLogic() = default;


    // Member Functions
    bool curvatureRequiresUpdate(const fvMesh&) const;
    void curvatureUpdated(const fvMesh&) const;
    void resizeBufferField(const label, const meshType&) const;
    void setBufferToZero(const meshType&) const;
    std::shared_ptr<fieldType> curvatureBuffer(const meshType&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CurvatureBufferLogicI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
