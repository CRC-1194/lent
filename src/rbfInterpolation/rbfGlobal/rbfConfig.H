/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Tomislav Maric, TU Darmstadt
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

    Global enumerations: 
        - Stencil types. 
        - Scaling radius types.

\*---------------------------------------------------------------------------*/

#ifndef rbfConfig_H
#define rbfConfig_H

namespace Foam
{
namespace RBF
{
    // Types of support radii for the RBF. 
    enum class supportType 
    {
        GLOBAL_SUPPORT, // N-N point connectivity. 
        MAX_CENTROID_RADIUS, 
        MEAN_CENTROID_RADIUS,
        RMS_RADIUS
    }; 

    // Types of support radii for the RBF. 
    enum class stencilType {BCC, BCC_FVM}; 

} // End namespace RBF
} // End namespace Foam 

// ************************************************************************* //

#endif
