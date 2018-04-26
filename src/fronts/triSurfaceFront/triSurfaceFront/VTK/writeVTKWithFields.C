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
    Foam::triSurfaceFront

SourceFiles
    diffuseInterfaceProperties.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description

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

#include "triSurfaceFront.H"
#include "triSurfaceFrontFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// This function has been taken over from the triSurface class. However, writing
// the region number to the output has been discarded for. If it turns out to
// be required again, look it up in "writeVTK.C" of the OpenFOAM source
void triSurfaceFront::writeVTKWithFields(Ostream& os) const
{
    // Write header
    os  << "# vtk DataFile Version 2.0" << nl
        << "triSurface" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA"
        << nl;

    const pointField& ps = points();

    os  << "POINTS " << ps.size() << " float" << nl;

    // Write vertex coords
    forAll(ps, pointi)
    {
        if (pointi > 0 && (pointi % 10) == 0)
        {
            os  << nl;
        }
        else
        {
            os  << ' ';
        }
        os  << ps[pointi].x() << ' '
            << ps[pointi].y() << ' '
            << ps[pointi].z();
    }
    os  << nl;

    os  << "POLYGONS " << triSurface::size() << ' ' << 4*triSurface::size() << nl;

    //labelList faceMap;
    //surfacePatchList patches(calcPatches(faceMap));

    for (label facei = 0; facei < triSurface::size(); ++facei)
    {
        if (facei > 0 && (facei % 10) == 0)
        {
            os  << nl;
        }
        else
        {
            os  << ' ';
        }
        os  << "3 "
            << triSurface::operator[](facei)[0] << ' '
            << triSurface::operator[](facei)[1] << ' '
            << triSurface::operator[](facei)[2];
    }
    os  << nl;

    // TODO: look up is complicated by different field types -->
    // cannot simply distinguish by point based and triangle based fields
    //
    // First, test with vector fields and if it works as intended, extend it to
    // support all field types
    auto frontCellVectorFieldMap = this->lookupClass<triSurfaceFrontVectorField>();
    auto fieldNames = frontCellVectorFieldMap.toc();

    label numFieldsToWrite = 0;

    forAll(fieldNames, I)
    {
        const auto& frontField = frontCellVectorFieldMap[fieldNames[I]];

        if (frontField->writeOpt() != IOobject::NO_WRITE)
        {
            numFieldsToWrite++;
        }
    }

    // Header for cell (triangle) based data
    os << "CELL_DATA " << this->localFaces().size() << nl
           << "FIELD frontTriangleFields " << numFieldsToWrite << nl;

    forAll(fieldNames, I)
    {
        const auto& frontFieldPtr = frontCellVectorFieldMap[fieldNames[I]];

        if (frontFieldPtr->writeOpt() != IOobject::NO_WRITE)
        {
            const auto& frontField = *frontFieldPtr;

            // Field header
            // FIXME: 3 only valid for vector quantities. Scalar quantities need 1
            os << frontField.name() << " 3 " << frontField.size() << " float" << nl;

            forAll(frontField, K)
            {
                if (K > 0 && (K % 10) == 0)
                {
                    os << nl;
                }
                else
                {
                    os << ' ';
                }
                os << frontField[K][0] << ' '
                       << frontField[K][1] << ' '
                       << frontField[K][2];
            }

            os << nl;
        }
    }

    // Same procedure for point based vector fields
    auto frontPointVectorFieldMap = this->lookupClass<triSurfaceFrontPointVectorField>();
    fieldNames = frontPointVectorFieldMap.toc();

    numFieldsToWrite = 0;

    forAll(fieldNames, I)
    {
        const auto& frontField = frontPointVectorFieldMap[fieldNames[I]];

        if (frontField->writeOpt() != IOobject::NO_WRITE)
        {
            numFieldsToWrite++;
        }
    }

    // Header for cell (triangle) based data
    os << "POINT_DATA " << this->localFaces().size() << nl
           << "FIELD frontTriangleFields " << numFieldsToWrite << nl;

    forAll(fieldNames, I)
    {
        const auto& frontFieldPtr = frontPointVectorFieldMap[fieldNames[I]];

        if (frontFieldPtr->writeOpt() != IOobject::NO_WRITE)
        {
            const auto& frontField = *frontFieldPtr;

            // Field header
            // FIXME: 3 only valid for vector quantities. Scalar quantities need 1
            os << frontField.name() << " 3 " << frontField.size() << " float" << nl;

            forAll(frontField, K)
            {
                if (K > 0 && (K % 10) == 0)
                {
                    os << nl;
                }
                else
                {
                    os << ' ';
                }
                os << frontField[K][0] << ' '
                       << frontField[K][1] << ' '
                       << frontField[K][2];
            }

            os << nl;
        }
    }
    /*
    if (writeSorted)
    {
        label faceIndex = 0;

        forAll(patches, patchi)
        {
            // Print all faces belonging to this patch

            for
            (
                label patchFacei = 0;
                patchFacei < patches[patchi].size();
                patchFacei++
            )
            {
                if (faceIndex > 0 && (faceIndex % 10) == 0)
                {
                    os  << nl;
                }
                else
                {
                    os  << ' ';
                }

                const label facei = faceMap[faceIndex++];

                os  << "3 "
                    << operator[](facei)[0] << ' '
                    << operator[](facei)[1] << ' '
                    << operator[](facei)[2];
            }
        }
        os  << nl;


        // Print region numbers

        os  << "CELL_DATA " << size() << nl;
        os  << "FIELD attributes 1" << nl;
        os  << "region 1 " << size() << " float" << nl;

        faceIndex = 0;

        forAll(patches, patchi)
        {
            for
            (
                label patchFacei = 0;
                patchFacei < patches[patchi].size();
                patchFacei++
            )
            {
                if (faceIndex > 0 && (faceIndex % 10) == 0)
                {
                    os  << nl;
                }
                else
                {
                    os  << ' ';
                }

                const label facei = faceMap[faceIndex++];

                os  << operator[](facei).region();
            }
        }
        os  << nl;
    }
    else
    {
        forAll(*this, facei)
        {
            if (facei > 0 && (facei % 10) == 0)
            {
                os  << nl;
            }
            else
            {
                os  << ' ';
            }
            os  << "3 "
                << operator[](facei)[0] << ' '
                << operator[](facei)[1] << ' '
                << operator[](facei)[2];
        }
        os  << nl;

        os  << "CELL_DATA " << size() << nl;
        os  << "FIELD attributes 1" << nl;
        os  << "region 1 " << size() << " float" << nl;

        forAll(*this, facei)
        {
            if (facei > 0 && (facei % 10) == 0)
            {
                os  << nl;
            }
            else
            {
                os  << ' ';
            }
            os  << operator[](facei).region();
        }
        os  << nl;
    }
    */
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
