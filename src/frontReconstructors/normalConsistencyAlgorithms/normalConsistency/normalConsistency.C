/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2.2.x                               
    \\  /    A nd           | Copyright held by original author
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
    Foam::diffuseInterfaceProperties

SourceFiles
    diffuseInterfaceProperties.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de
    Mathematical Modeling and Analysis
    Center of Smart Interfaces, TU Darmstadt

Description
    
    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

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
    defineRunTimeSelectionTable(normalConsistency, Dictionary);
    addToRunTimeSelectionTable(normalConsistency, normalConsistency, Dictionary);

// * * * * * * * * * * * * * Private  member functions * * * * * * * * * * * //
void normalConsistency::runNormalConsistencyAlgorithm(
    triSurfaceFront& front,
    const volScalarField& signedDistance,
    const pointScalarField& pointSignedDistance // Not used.
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
    const word name = configDict.lookup("type");

    DictionaryConstructorTable::iterator cstrIter =
        DictionaryConstructorTablePtr_->find(name);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "normalConsistency::New(const word& name)"
        )   << "Unknown normalConsistency type "
            << name << nl << nl
            << "Valid normalConsistencys are : " << endl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return tmp<normalConsistency> (cstrIter()(configDict));
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
