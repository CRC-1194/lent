/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::distanceNormalConsistency

SourceFiles
    distanceNormalConsistency.C

Authors
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description

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


#include "normalConsistency.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "lentCommunication.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(normalConsistency, 0);
    defineRunTimeSelectionTable(normalConsistency, Dictionary)
    addToRunTimeSelectionTable(normalConsistency, normalConsistency, Dictionary);

// * * * * * * * * * * * * * Private  member functions * * * * * * * * * * * //
void normalConsistency::runNormalConsistencyAlgorithm(
    triSurfaceFront& front,
    const volScalarField& signedDistance,
    const pointScalarField&
) const
{
    // Gradient based normal consistency algorithm.
    volVectorField distGrad{fvc::grad(signedDistance)};

    List<labelledTri>& triangles = static_cast<List<labelledTri>& > (front);
    const vectorField& triangleNormals = front.faceNormals();

    const fvMesh& mesh = signedDistance.mesh(); 
    const lentCommunication& communication = 
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
    ); 

    const auto& triangleToCell = communication.triangleToCell();  

    // For all faces
    forAll (triangles, triangleI)
    {
        scalar gradMag = mag(distGrad[triangleToCell[triangleI]]);

        if (gradMag >= SMALL)
        {
            distGrad[triangleToCell[triangleI]] /= gradMag;

            scalar normalMag = mag(triangleNormals[triangleI]);

            if (normalMag > SMALL)
            {
                vector triangleNormal = triangleNormals[triangleI] / mag(triangleNormals[triangleI]);

                if ((triangleNormal & distGrad[triangleToCell[triangleI]]) < 0)
                {
                    triangles[triangleI].flip();
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

tmp<normalConsistency>
normalConsistency::New(const dictionary& configDict)
{
    const word name = configDict.get<word>("type");

    // Find the constructor pointer for the model in the constructor table.
    auto* ctorPtr = DictionaryConstructorTable(name);

    // If the constructor pointer is not found in the table.
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            configDict,
            "normalConsistency",
            name,
            *DictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    // Construct the model and return the autoPtr to the object.
    return tmp<normalConsistency> (ctorPtr(configDict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void normalConsistency::makeFrontNormalsConsistent(
    triSurfaceFront& front,
    const volScalarField& signedDistance,
    const pointScalarField& pointSignedDistance // Not used by this algorithm. TM. 
) const
{
    runNormalConsistencyAlgorithm(front, signedDistance, pointSignedDistance);

    // Ensure update of demand driven data by completely clearing it
    front.clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
