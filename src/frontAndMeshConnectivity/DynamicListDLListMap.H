/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::DynamicListDLListMap

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

Description

SourceFiles
    DynamicListDLListMap.C

\*---------------------------------------------------------------------------*/

#ifndef DynamicListDLListMap_H
#define DynamicListDLListMap_H

#include "DynamicList.H"
#include "DLList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace frontTracking
{

/*---------------------------------------------------------------------------*\
                         Class DynamicListDLListMap Declaration
\*---------------------------------------------------------------------------*/

class DynamicListDLListMap
{
    // Private data

    DynamicList<DLList<label> > map_; 

public:

    typedef DynamicList<DLList<label> >::const_iterator const_iterator;
    typedef DynamicList<DLList<label> >::iterator iterator;

    // Constructors
    
    //- Construct from components
    DynamicListDLListMap(label size=0);

    //- Copy construction
    DynamicListDLListMap(const DynamicListDLListMap&);

    //- Destructor
    ~DynamicListDLListMap();

    // Member Functions 

    DynamicList<DLList<label > >::const_iterator
    begin()
    {
        return map_.begin(); 
    }

    DynamicList<DLList<label > >::const_iterator
    end()
    {
        return map_.end(); 
    }

    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frontTracking
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DynamicListDLListMap.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
