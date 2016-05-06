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
    Foam::frontMotionSolver

SourceFiles
    frontMotionSolver.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Interface for the front motion solution.

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#include "frontMotionSolver.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontMotionSolver, 0);
    defineRunTimeSelectionTable(frontMotionSolver, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontMotionSolver::frontMotionSolver(const dictionary& configDict)
    :
        cellDisplacementTmp_(),
        frontDisplacementTmp_(),
        interpolation_(configDict)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

tmp<frontMotionSolver>
frontMotionSolver::New(const dictionary& configDict)
{
    const word name = configDict.lookup("type");

    DictionaryConstructorTable::iterator cstrIter =
        DictionaryConstructorTablePtr_->find(name);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "frontMotionSolver::New(const word& name)"
        )   << "Unknown frontMotionSolver type "
            << name << nl << nl
            << "Valid frontMotionSolvers are : " << endl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return tmp<frontMotionSolver> (cstrIter()(configDict));
}


// * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

void frontMotionSolver::initDisplacements(
    const triSurfaceFront& front,
    const volVectorField& cellVelocity
)
{
    const auto& mesh = cellVelocity.mesh(); 
    const auto& runTime = mesh.time(); 

    dimensionedVector zeroDisplacement(
        "zero",
        dimLength, 
        vector(0,0,0)
    );

    if (cellDisplacementTmp_.empty())
    {
        cellDisplacementTmp_ = tmp<volVectorField>(
            new volVectorField( 
                IOobject(
                    "cellDisplacement",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                zeroDisplacement
            )
        );
    } 

    if (frontDisplacementTmp_.empty())
    {
        frontDisplacementTmp_ = tmp<triSurfaceFrontPointVectorField>(
            new triSurfaceFrontPointVectorField( 
                IOobject(
                    "frontDisplacement",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                front,
                zeroDisplacement
            )
        );
    } else
    {
        frontDisplacementTmp_->resize(front.nPoints()); 
        //frontDisplacementTmp_.ref() = zeroDisplacement;  
    }
}

void frontMotionSolver::evolveFront(
    triSurfaceFront& front,
    const volVectorField& cellVelocity 
)
{
    initDisplacements(front,cellVelocity); 

    // Overload this member function for different functionality.
    calcCellDisplacement(cellVelocity); 

    // Interpolate the displacement field from the cell to the front.
    const auto& deltaC = cellDisplacements();
    auto& deltaF = frontDisplacements();
    interpolation_.interpolate(deltaC, deltaF); 

    // Displace front points with front displacements.  
    pointField& frontPoints = const_cast<pointField&>(front.points());
    frontPoints += deltaF; 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
