/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | Copyright (C) 2018 Tomislav Maric, TU Darmstadt
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

Author
    Tomislav Maric maric@mma.tu-darmstadt.de

Description
    Helper function templates for the rbfCellsInterpolationEigen test. 

\*---------------------------------------------------------------------------*/

#ifndef lentTestrbfCellsInterpolationEigenTemplates_H
#define lentTestrbfCellsInterpolationEigenTemplates_H

std::pair<label, label> findPatchLabel
(
    label faceLabel, 
    const polyBoundaryMesh& boundaryMesh
) 
{
    std::pair<label, label> patchLabelFaceLabel = std::make_pair(-1, -1);  
    forAll(boundaryMesh, patchI)
    {
        const auto& patch = boundaryMesh[patchI]; 
        if  ((faceLabel >= patch.start()) && 
             (faceLabel < patch.start() + patch.size()))
        {
            patchLabelFaceLabel.first = patchI; 
            patchLabelFaceLabel.second = faceLabel - patch.start(); 
            return patchLabelFaceLabel; 
        }
    }
    return patchLabelFaceLabel;
}

template<typename Field, typename Points, typename Surface> 
void surfaceSetField(
    Field& field, 
    Points const& points, 
    Surface const& surface
)
{
    forAll(points, pI)
    {
        field[pI] = surface.value(points[pI]);
    }
}

template<typename GeoField, typename GeoMesh>
tmp<GeoField> 
makeField
(
    const char* name, 
    GeoMesh const& mesh,
    dimensioned<typename GeoField::value_type> value =  
    dimensioned<typename GeoField::value_type>
    (
        "zero", 
        dimless, 
        pTraits<typename GeoField::value_type>::zero
    )
 )
{
    return tmp<GeoField>
    (
        new GeoField 
        (
            IOobject
            (
                name,
                mesh.time().timeName(), 
                mesh.time(),  
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ), 
            mesh, 
            value,
            "zeroGradient"
        )
    );
}

template<typename Type> 
decltype(auto) l2norm(Type&& t) 
{
    return Foam::sqrt(Foam::sum(Foam::sqr(std::forward<Type>(t)))); 
};

template<typename PointListList> 
void genRandomPointsInCells
(
    PointListList& randPointsInCells, 
    const fvMesh& mesh, 
    const label nPointsPerCell
)
{
    std::random_device rdint;  
    std::mt19937 genint(rdint()); 

    std::random_device rdreal;  
    std::mt19937 genreal(rdreal()); 
    std::uniform_real_distribution<> disreal(0, 1.0);

    const auto& meshCells = mesh.cells(); 
    const auto& meshFaces = mesh.faces(); 
    const auto& meshPoints = mesh.points(); 
    const auto& cellCenters = mesh.C(); 
    const auto& faceCenters = mesh.Cf(); 
    forAll(meshCells, cellI)
    {
        const auto& meshCell = meshCells[cellI]; 
        const label cellN = meshCell.size() - 1;
        std::uniform_int_distribution<> distcell(0, cellN); 

        for(label pointI = 0; pointI < nPointsPerCell; ++pointI)
        {
            // Pick a random local face label in a cell. 
            const label faceR = distcell(genint); 
            // Get the global face label. 
            const label faceGl = meshCell[faceR];  


            vector faceCenter (-VGREAT, -VGREAT, -VGREAT); 
            if (mesh.isInternalFace(faceGl))
            {
                faceCenter = faceCenters[faceGl];
            }
            else
            {
                // Find the boundary patch ID and the local face ID from the mesh. 
                const auto& boundaryMesh = mesh.boundaryMesh(); 
                const std::pair<label, label> patchLabelFaceLabel = 
                    findPatchLabel(faceGl, boundaryMesh); 

                const label& patchP = patchLabelFaceLabel.first;
                const label& faceP = patchLabelFaceLabel.second; 

                // Get the face center 
                const auto& faceCentersBoundaryField = faceCenters.boundaryField();  
                const auto& faceCentersPatchField = faceCentersBoundaryField[patchP]; 
                faceCenter = faceCentersPatchField[faceP];
            }

            // Pick a random point on the cell face.
            const auto& meshFace = meshFaces[meshCell[faceR]];
            const label faceN = meshFace.size() - 1; 
            std::uniform_int_distribution<> distface(0, faceN); 
            label pointR = distface(genint); 
            // Next point on the face respective to the randomly chosen point. 
            label pointS = pointR + 1; 
            // If the next point label is greater than the size of the face.
            if (pointS > faceN) 
                pointS = 0; // First face label (0) is the next point label.

            // Global face point labels. 
            const label pointGl0 = meshFace[pointR]; 
            const label pointGl1 = meshFace[pointS]; 

            // Point is in the mesh cell with the cell label cellID. 
            auto testPointInCell = [&mesh, &meshCells, &meshFaces, &meshPoints]
            (
                const vector& point, 
                const label cellID, 
                const char* pointName
            )
            {
                bool pointIsInCell = 
                    mesh.pointInCell(point, cellID);

                if (!pointIsInCell)
                {
                    Info << pointName << " not in cell, " 
                        << pointName << " = " << point << endl;

                    const auto cellPoints = meshCells[cellID].points(meshFaces, meshPoints);  

                    Info << "CellID = " << cellID << endl
                        << "cell bounding box = " 
                        << boundBox(min(cellPoints), max(cellPoints)) << endl; 

                }

                return pointIsInCell;
            };
            
            const auto& edgePoint0 = meshPoints[pointGl0];  
            const auto& edgePoint1 = meshPoints[pointGl1];

            testPointInCell(faceCenter, cellI, "faceCenter"); 
            testPointInCell(edgePoint0, cellI, "edgePoint0"); 
            testPointInCell(edgePoint1, cellI, "edgePoint1"); 

            // Generate random barycentric coordinates for the randomly chosen tetrahedron. 
            scalarField alphas(4, 0.); 
            forAll(alphas, cmptI)
                alphas[cmptI] = disreal(genreal);
            // Normalize the barycentric coordinates to keep the point in the tetrahedron. 
            alphas /= sum(alphas); 

            assert((sum(alphas) - 1. )/ sum(alphas) < 
                    2*std::numeric_limits<double>::epsilon()); 

            // Set the random point in a random tetrahedron.
            randPointsInCells[cellI][pointI] = 
                alphas[0]*cellCenters[cellI] +
                alphas[1]*faceCenter + 
                alphas[2]*edgePoint0 +  
                alphas[3]*edgePoint1;
            
            const auto& randPoint = randPointsInCells[cellI][pointI]; 

            auto randPointName = 
                Foam::string("randomPoint") + std::to_string(pointI);
            assert(testPointInCell(randPoint, cellI, randPointName.c_str()));
        }
    }

}


template
<
    typename Surface, 
    typename PointListLists, 
    typename RbfCellInterpolation
>
void evaluateLinfErrors
(
    volScalarField& Linfs, 
    PointListLists const& randPointsInCells, 
    Surface const& surface,
    RbfCellInterpolation const& cellRbfs
)
{
    forAll(Linfs, cellI)
    {
        const auto& randPoints = randPointsInCells[cellI]; 
        for(const auto& randPoint : randPoints)
        {
            const double exactVal = 
                surface.value(randPoint);


            const double errBcc = std::abs
            (
                cellRbfs.value(randPoint, cellI) - exactVal
            );

            if (errBcc > Linfs[cellI])
                Linfs[cellI] = errBcc; 
        }
    }
}

#endif

// ************************************************************************* //
