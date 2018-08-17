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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
template<class fieldType, class meshType, class tensorType>
void CurvatureBufferLogic<fieldType, meshType, tensorType>::initializeBufferField(const meshType& mesh) const
{
    if (!bufferInitialized_)
    {
        const Time& runTime = mesh.time();  

        curvatureBufferPtr_ = std::shared_ptr<fieldType>{ 
                new fieldType
                (
                    IOobject(
                        bufferFieldName_, 
                        runTime.timeName(), 
                        mesh,
                        IOobject::NO_READ, 
                        IOobject::AUTO_WRITE
                    ), 
                    mesh, 
                    dimensioned<tensorType>(
                        "zero", 
                        dimless/dimLength, 
                        pTraits<tensorType>::zero
                    )
                )
        };

        bufferInitialized_ = true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class fieldType, class meshType, class tensorType>
CurvatureBufferLogic<fieldType, meshType, tensorType>::CurvatureBufferLogic(const word& bufferFieldName)
:
    bufferFieldName_{bufferFieldName}
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class fieldType, class meshType, class tensorType>
bool CurvatureBufferLogic<fieldType, meshType, tensorType>::curvatureRequiresUpdate(const fvMesh& mesh) const
{
    return (mesh.time().timeIndex() != lastTimeUpdated_);
}

template<class fieldType, class meshType, class tensorType>
void CurvatureBufferLogic<fieldType, meshType, tensorType>::curvatureUpdated(const fvMesh& mesh) const
{
    lastTimeUpdated_ = mesh.time().timeIndex();
}

template<class fieldType, class meshType, class tensorType>
std::shared_ptr<fieldType> CurvatureBufferLogic<fieldType, meshType, tensorType>::curvatureBuffer(const meshType& mesh) const
{
    initializeBufferField(mesh);

    return curvatureBufferPtr_;
}

template<class fieldType, class meshType, class tensorType>
void CurvatureBufferLogic<fieldType, meshType, tensorType>::resizeBufferField(const label size, const meshType& mesh) const
{
    initializeBufferField(mesh);

    if (curvatureBufferPtr_->size() != size) 
    {
        curvatureBufferPtr_->resize(size);
    }
}

template<class fieldType, class meshType, class tensorType>
void CurvatureBufferLogic<fieldType, meshType, tensorType>::setBufferToZero(const meshType& mesh) const
{
    initializeBufferField(mesh);

    *curvatureBufferPtr_ = dimensioned<tensorType>
                            {
                                "zero",
                                curvatureBufferPtr_->dimensions(),
                                pTraits<tensorType>::zero
                            };
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
