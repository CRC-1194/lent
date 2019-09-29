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

#ifndef lentTestRbfReconstructionHelpers_H
#define lentTestRbfReconstructionHelpers_H

#include "foamVtkLegacyAsciiFormatter.H"
#include "foamVtkOutput.H"
#include <sstream>
#include <iomanip>
#include <numeric>

// Legacy VTK output functions.
template<typename Points> 
void write_vtk_points(
        std::ostream& vtks,
        Points const& points, 
        const std::string& fileName,
        int index = 0 
)
{
    // Initialize the legacy ascii VTK formatter: 15 decimal point precision.
    Foam::vtk::legacyAsciiFormatter legacyFormat(vtks, 15); 

    // Write the file header based on the chosen format (legacy VTK).
    Foam::vtk::legacy::fileHeader(legacyFormat, fileName, Foam::vtk::fileTag::POLY_DATA);

    // Write point header
    Foam::vtk::legacy::beginPoints(vtks, points.size());

    // Write point coordinates. 
    for(const auto& point : points) 
    {
        for(char cI = 0; cI < 3; ++cI)
            vtks << point[cI] << " "; 
        vtks << "\n"; 
    }
    vtks << "\nPOINT_DATA " << points.size() << "\n";
}

template<typename Stream, typename Points> 
std::enable_if_t<std::is_base_of_v<std::ostream, Stream>, Stream>
vtk_stream( 
    std::string fileName,
    Points const& points, 
    int fileIndex = 0
)
{
    std::stringstream fileNameStream; 
    fileNameStream << fileName << "-" 
        << std::setw(3) << std::setfill('0') <<  fileIndex;  
    std::ofstream vtks (fileNameStream.str() + ".vtk"); 
    write_vtk_points(vtks, points, fileName, fileIndex);
    return vtks;
}

template<typename Stream, typename Scalars> 
std::enable_if_t<std::is_arithmetic_v<typename Scalars::value_type>>
write_to_vtk_stream(
        Stream& vtks,
        Scalars const& values, 
        const std::string& fieldName,
        int index = 0 
)
{
    // Format for legacy VTK scalar data output:
    // SCALARS volume_scalars float nCmpts
    // LOOKUP_TABLE default
    vtks << "SCALARS " << fieldName << " float " << "1\n"; 
    vtks << "LOOKUP_TABLE default\n";

    // Write point values.
    for (const auto val : values)
        vtks << val << "\n"; 
}


template<typename Stream, typename Vectors> 
std::enable_if_t<std::is_arithmetic_v<typename Vectors::value_type::value_type>>
write_to_vtk_stream(
        Stream& vtks,
        Vectors const& vecs, 
        const string& fieldName,
        int index = 0 
)
{
    // Format for legacy VTK vector data output:
    // VECTORS pnorm float
    vtks << "VECTORS " << fieldName << " float " << "1\n"; 

    // Write point values.
    for (const auto& vec : vecs)
    {
        for (const auto cmpt : vec)
            vtks << cmpt << " "; 
        vtks << "\n";
    }
}

template<typename Field, typename Points, typename Surface> 
void surfaceSetField(
    Field& field, 
    Points const& points, 
    Surface const& surface
)
{
    forAll(points, pI)
        field[pI] = surface.value(points[pI]);
}

struct surfaceTestFields
{
    const fvMesh& mesh_;
    const word name_; 
    dimensionedScalar initLinf_; 
    pointMesh pMesh_;  

    pointScalarField pointValues_;
    volScalarField volValues_;  

    volScalarField volLinfBcc_; 
    volScalarField volLinfBccc_; 

    public: 

        surfaceTestFields(const fvMesh& mesh, const word& name)
            :
                mesh_(mesh),
                name_(name), 
                initLinf_(dimensionedScalar("Linf", dimLength, -VGREAT)),
                pMesh_(mesh), 
                pointValues_ 
                (
                    IOobject
                    (
                        "pointValues" + name_, 
                        mesh.time().timeName(), 
                        mesh, 
                        IOobject::NO_READ, 
                        IOobject::AUTO_WRITE
                    ),
                    pMesh_, 
                    dimensionedScalar(name_, dimless, 0)
                ),
                volValues_ 
                (
                    IOobject
                    (
                        "volValues" + name_,  
                        mesh.time().timeName(), 
                        mesh, 
                        IOobject::NO_READ, 
                        IOobject::AUTO_WRITE
                    ),
                    mesh, 
                    dimensionedScalar(name_, dimless, 0)
                ), 
                volLinfBcc_ 
                (
                    IOobject
                    (
                        "volLinfBcc" + name_,  
                        mesh.time().timeName(), 
                        mesh, 
                        IOobject::NO_READ, 
                        IOobject::AUTO_WRITE
                    ),
                    mesh, 
                    initLinf_
                ),
                volLinfBccc_
                (
                    IOobject
                    (
                        "volLinfBccc" + name_,  
                        mesh.time().timeName(), 
                        mesh, 
                        IOobject::NO_READ, 
                        IOobject::AUTO_WRITE
                    ),
                    mesh, 
                    initLinf_
                )
        {}

        template<typename Surface>
        surfaceTestFields(const fvMesh& mesh, Surface const& surface, const word& name)
        :
            surfaceTestFields(mesh, name)
        {
            setValues(surface);  
        }

        void resetErrors()
        {
            volLinfBcc_ == initLinf_; 
            volLinfBccc_ == initLinf_; 
        }

        template<typename Surface>
        void setVolValues(Surface const& surface)
        {
            surfaceSetField(volValues_, mesh_.C(), surface);  
        }

        template<typename Surface> 
        void setPointValues(Surface const& surface)
        {
            surfaceSetField(pointValues_, mesh_.points(), surface);
        }

        template<typename Surface> 
        void setValues(Surface const& surface)
        {
            setVolValues(surface); 
            setPointValues(surface); 
        }

        const fvMesh& mesh()
        {
            return mesh_;
        }

        void writeValueFields() const
        {
            pointValues_.write(); 
            volValues_.write();
        }
};

struct errorNorms 
{
    using errorList = std::vector<double>;

    errorList errors_; 
    errorList absErrors_;

    template<typename Points> 
    errorNorms(Points const& points)
    :
        errors_(),
        absErrors_()
    {
        errors_.reserve(points.size());
    }

    void append(double error)
    {
        errors_.push_back(error);
        absErrors_.push_back(std::abs(error));
    }

    double Linf()
    {
        double Linf = std::numeric_limits<double>::max(); 

        for (const auto& absErr : absErrors_)
            if (absErr < Linf) Linf = absErr;

        return Linf;
    }

    double L1()
    {
        double absSum = 0;  
        for(const auto& absErr : absErrors_)
            absSum += absErr; 

        return absSum / absErrors_.size(); 
            
    }

    double L2()
    {
        double sqrSum = 0; 
        for (const auto& absErr : absErrors_)
            sqrSum += absErr; 

        return std::sqrt(sqrSum / absErrors_.size()); 
    }

    const errorList& errors() const
    {
        return errors_; 
    }

    const errorList& absErrors() const
    {
        return errors_; 
    }
};

template <typename Surface, typename Points>
errorNorms calcPointErrors(
    Points const& points, 
    Surface const& surface
)
{
    errorNorms pointErrors(points); 

    for (const auto& point: points)
        pointErrors.append(surface.value(point)); 

    return pointErrors;  
}

template
<
    typename Surface, 
    typename SurfaceFields, 
    typename Outstream
> 
void testLinearIsoPoints(
    Surface const& surface, 
    SurfaceFields const& surfFields,
    Outstream& os
)
{
    using clock = std::chrono::high_resolution_clock;
    using sec = std::chrono::duration<double>;
    Info << "Linear iso-surface reconstruction for a " << surface.name() << endl; 
    isoPointCalculator linearCalc; 
    auto t0 = clock::now();
    linearCalc.calcContourPoints(surfFields.volValues_, surfFields.pointValues_); 
    auto t1 = clock::now();
    sec seconds = t1 - t0;
    std::cout << "done in " << seconds.count() << " seconds." << std::endl;

    auto linEdgeErrors = calcPointErrors(linearCalc.edgePoints(), surface);
    auto edgeVtkStream = vtk_stream<std::ofstream>(
        "linear-edgePoints-" + surface.name(), 
        linearCalc.edgePoints()
    );
    write_to_vtk_stream(edgeVtkStream, linEdgeErrors.absErrors(), "edgeLinf");
    
    auto linCellErrors = calcPointErrors(linearCalc.contourPoints(), surface);

    auto cellVtkStream = vtk_stream<std::ofstream>(
        "linear-cellPoints-" + surface.name(), 
        linearCalc.contourPoints()
    );

    write_to_vtk_stream(cellVtkStream, linCellErrors.errors(), "cellLinf");

    os << surface.name() << "," 
        << linEdgeErrors.Linf() << "," << linEdgeErrors.L1() << "," << linEdgeErrors.L2() 
        << "," 
        << linCellErrors.Linf() << "," << linCellErrors.Linf() << "," << linCellErrors.Linf() 
        << seconds.count() << "\n";
    surfFields.writeValueFields();
}

template
<
    typename Surface, 
    typename SurfaceFields, 
    typename MeshRbfs, 
    typename Outstream
>
void testRbfIsoPoints(
    Surface const& surface, 
    SurfaceFields const& surfFields,
    MeshRbfs const& rbfs, 
    Outstream& os,
    std::string stencilName
)
{
    using clock = std::chrono::high_resolution_clock;
    using sec = std::chrono::duration<double>;
    Info << "RBF iso-surface reconstruction for a " << surface.name() << endl; 
    // FIXME: Make the RBF calculator use the MeshRbfs as an attribute.   
    // Otherwise, calcPointErrors(rbfCalc, rbfs, surface) : error if different
    // rbfs are used. TM.
    rbfIsoPointCalculator rbfCalc; 
    auto t0 = clock::now();
    rbfCalc.calcContourPoints(surfFields.volValues_, surfFields.pointValues_, rbfs); 
    auto t1 = clock::now();
    sec seconds = t1 - t0;
    std::cout << "done in " << seconds.count() << " seconds." << std::endl;
    auto rbfCellErrors = calcPointErrors(rbfCalc.contourPoints(), surface);

    std::string baseFileName = "rbf-" + rbfs.rbfName() + "-cellPoints-" + 
        surface.name() + "-" + stencilName;

    auto vtkFile  = vtk_stream<std::ofstream>(baseFileName, rbfCalc.contourPoints()); 
    write_to_vtk_stream(vtkFile, rbfCellErrors.errors(), "L_INF"); 

    for (int iter = 1; iter < 20; ++iter)
    {
        rbfCalc.rbfCorrectContourPoints(rbfs);
        rbfCellErrors = calcPointErrors(rbfCalc.contourPoints(), surface);

        auto vtkFileI = vtk_stream<std::ofstream>(
            baseFileName, 
            rbfCalc.contourPoints(), 
            iter
        ); 
        write_to_vtk_stream(vtkFileI, rbfCellErrors.errors(), "L_INF"); 
    }
    
    os << rbfs.rbfName() << "," << stencilName << "," 
        << surface.name() << "," << rbfCellErrors.Linf() 
        << "," << seconds.count() << ",";  
}

template 
<
    typename RbfKernel, 
    typename Surface, 
    typename Outstream
> 
void testRbfReconstruction
(
    Surface const& surface,  
    surfaceTestFields& surfaceFields,
    Outstream& os 
)
{
    using clock = std::chrono::high_resolution_clock;
    using sec = std::chrono::duration<double>;

    Info << "RBF Kernel = " << RbfKernel::name() << endl;

    // BCC
    Info << "Factorizing BCC RBF systems... \n"; 
    const auto& mesh = surfaceFields.mesh();
    auto t0 = clock::now();
    rbfCellsInterpolationEigen<RbfKernel> cellRbfsBcc(mesh, stencilType::BCC); 
    auto t1 = clock::now();
    sec factorSeconds = t1 - t0; 
    std::cout << "BCC factorization CPU time in seconds = " 
        << factorSeconds.count() << " seconds." << std::endl;
    Info << "Solving for field values with BBC RBFs... \n"; 
    auto t2 = clock::now();
    cellRbfsBcc.solve(surfaceFields.volValues_, surfaceFields.pointValues_); 
    auto t3 = clock::now();
    sec solveSeconds = t3 - t2;
    std::cout << "BCC solution CPU time in seconds = " 
        << solveSeconds.count() << " seconds." << std::endl;

    testRbfIsoPoints(surface, surfaceFields, cellRbfsBcc, os, "BCC");
    os << factorSeconds.count() << "," << solveSeconds.count() << "\n";
    // BCC
    
    // BCCC
    Info << "Factorizing BCCC RBF systems... \n"; 
    t0 = clock::now();
    rbfCellsInterpolationEigen<RbfKernel> cellRbfsBccc(mesh, stencilType::BCCC); 
    t1 = clock::now();
    factorSeconds = t1 - t0;
    std::cout << "BCCC stencil factorization CPU time in seconds = " 
        << factorSeconds.count() << " seconds." << std::endl;
    Info << "Solving for field values with BCCC RBFs... \n"; 
    t2 = clock::now();
    cellRbfsBccc.solve(surfaceFields.volValues_, surfaceFields.pointValues_); 
    t3 = clock::now();
    solveSeconds = t3 - t2;
    std::cout << "BCCC solution CPU time in seconds = " 
        << solveSeconds.count() << " seconds." << std::endl;
    testRbfIsoPoints(surface, surfaceFields, cellRbfsBccc, os, "BCCC");
    os << factorSeconds.count() << "," << solveSeconds.count() << "\n";
    // BCCC

    return; 
}

template 
<
    typename RbfTuple, 
    int N = 0,
    typename Surface, 
    typename Outstream
> 
void rbfReconstructLoop(
    Surface const& surface,  
    surfaceTestFields& surfaceFields,
    Outstream& os 
)
{
    testRbfReconstruction<std::tuple_element_t<N, RbfTuple>>(surface,surfaceFields,os);

    if constexpr (N + 1 <  std::tuple_size_v<RbfTuple>)
        rbfReconstructLoop<RbfTuple, N + 1>(surface,surfaceFields,os);
}

#endif

// ************************************************************************* //
