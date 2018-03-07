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
    Foam::lentDummyTest

Description

SourceFiles
    lentDummyTestI.H
    lentDummyTest.C
    lentDummyTestIO.C

\*---------------------------------------------------------------------------*/

#ifndef lentDummyTest_H
#define lentDummyTest_H

#include "lentSubalgorithmTest.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class lentDummyTest Declaration
\*---------------------------------------------------------------------------*/

class lentDummyTest
:
    public lentSubalgorithmTest
{
    // Private data
    tmp<volScalarField> referenceFieldTmp_;
    tmp<volScalarField> computedFieldTmp_;


    // Private Member Functions
    virtual void randomSetup();
    virtual void perturbInputFields(){};
    virtual void computeApproximatedFields();
    virtual void evaluateMetrics();


public:

    // Constructors
    lentDummyTest(const fvMesh& mesh, triSurfaceFront& front);


    //- Destructor
    virtual ~lentDummyTest() = default;


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
