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
    Foam::FrontTracking::frontSmoother

SourceFiles
    frontSmoother.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description
    Implementation of two volume-conservative smoothing algorithms for
    triangulated surfaces as proposed by Kuprat et al 2001
    Both algorithms have been extended so they are applicable to fronts
    with boundaries, e.g. in 2D cases or contact line problems

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

#ifndef frontSmoother_H
#define frontSmoother_H

#include <vector>

#include "volFieldsFwd.H"
#include "triSurfaceFront.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class frontSmoother Declaration
\*---------------------------------------------------------------------------*/
using vectorTuple = FixedList<vector,2>;

class frontSmoother
{

    // Private data
    scalar underrelaxationFactor_;
    label nSweeps_;
    word smoothingType_;

    // Private Member Functions
    vector computeA(const label& pointLabel, const triSurfaceFront& front) const;
    vector computeV(const label& edgeLabel, const triSurfaceFront& front) const;
    vectorTuple smoothEdge(const edge& relaxEdge, const triSurfaceFront& front) const;
    vector starBarycentre(const label& pointLabel, const triSurfaceFront& front) const;

    // TODO: this functionality should be moved to the triSurfaceFront class (TT)
    std::vector<label> singlePointBoundaryEdges(const triSurfaceFront& front) const;
    std::vector<label> internalEdges(const label nInternalEdges, const std::vector<label>& boundaryEdges) const;
    void updateLocalPoints(triSurfaceFront& front) const;
    void updateGlobalPoints(triSurfaceFront& front) const;
    label determineBoundaryPoint(const edge& relaxEdge, const triSurfaceFront& front) const;
    Tensor<scalar> movementRestriction(const label& pointLabel, const triSurfaceFront& front, const fvMesh& mesh) const;
    label containingFace(const label& pointLabel, const triSurfaceFront& front, const fvMesh& mesh) const;
    label containingEdge(const label& pointLabel, const label& faceLabel, const triSurfaceFront& front, const fvMesh& mesh) const;
    bool boundaryFacesAreCoplanar(const label& edgeLabel, const triSurfaceFront& front, const fvMesh& mesh) const;


public:

    // Constructors
    frontSmoother(const dictionary& configDict);


    // Member Functions
    // Algorithm no. 4 of Kuprat paper
    void smoothEdges(triSurfaceFront& front, const fvMesh& mesh) const;
    // Algorithm no. 3 of Kuprat paper
    void smoothPoints(triSurfaceFront& front, const fvMesh& mesh) const;

    // Invokes algorithm 3, algorithm 4 or both of them depending on the
    // "smooth" keyword specified in the lentSolution dictionary
    void smoothFront(triSurfaceFront& front, const fvMesh& mesh) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
