/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) held by orignal authors.
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

Application
    lentTestCurvatureModels 

Description

    Test LENT curvature models against the exact curvature. 

Authors
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "incompressibleTwoPhaseMixture.H"
#include "turbulenceModel.H"
#include "lentMethod.H"
#include "volMesh.H"
#include <fstream>

using namespace FrontTracking;

namespace Foam {

auto true_ref_lambda = [](const auto& x) { return true; };

template
<
    typename Type, 
    template<typename> class PatchField, 
    typename Mesh, 
    typename FilterFunction
>
tmp<GeometricField<Type, PatchField,Mesh>>  
mapToOnes(
    GeometricField<Type, PatchField, Mesh> const& vf, 
    FilterFunction filter = true_ref_lambda 
)
{
    tmp<GeometricField<Type, PatchField, Mesh> > resultTmp(
        new GeometricField<Type, PatchField, Mesh>(
            IOobject(
                "ones", 
                vf.time().timeName(), 
                vf.mesh(), 
                IOobject::NO_READ, 
                IOobject::AUTO_WRITE
            ),
            vf.mesh(), 
            dimensioned<Type>("zero", dimless, pTraits<Type>::zero)
        )
    ); 

    GeometricField<Type, PatchField, Mesh>& result = resultTmp(); 

    forAll(vf, cellI)
    {
        if (filter(vf[cellI]))
        {
            result[cellI] = pTraits<Type>::one;
        }
    }

    return resultTmp;
}; 

template
<
    typename Type, 
    template<typename> class PatchField, 
    typename Mesh, 
    typename FilterFunction
>
tmp<GeometricField<Type, PatchField, Mesh>>  
mapToOnes(tmp<GeometricField<Type, PatchField, Mesh> > vfTmp, FilterFunction filter)
{
    return mapToOnes(vfTmp(), filter); 
}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    std::fstream errorFile; 

    errorFile.open("curvatureErrors.dat", std::ios_base::app);

    triSurfaceFront front(
        IOobject(
            "front.stl",
            "front",
            runTime,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    triSurfaceMesh frontMesh(front); 

    lentMethod lent(front, mesh);

    lent.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);

    runTime++;

    lent.reconstructFront(front, signedDistance, pointSignedDistance);

    const dictionary& lentDict = lent.dict(); 

    tmp<frontCurvatureModel> exactCurvatureModelTmp = frontCurvatureModel::New(lentDict.subDict("exactCurvatureModel")); 
    const frontCurvatureModel& exactCurvatureModel = exactCurvatureModelTmp(); 
    tmp<volScalarField> cellCurvatureExactTmp = exactCurvatureModel.cellCurvature(mesh,frontMesh);  
    volScalarField& exactCurvature = cellCurvatureExactTmp();  

    const frontCurvatureModel& numericalCurvatureModel = lent.curvatureModel();  
    tmp<volScalarField> numericalCurvatureTmp = numericalCurvatureModel.cellCurvature(mesh,frontMesh);  
    volScalarField& numericalCurvature = numericalCurvatureTmp();  
    numericalCurvature.rename("numericalCurvature"); 

    // FIXME: How to filter out the curvature for testing? PEqn uses snGrad for this. How is the surface
    //volScalarField onesFilter (mapToOnes(mag(fvc::grad(markerField)), [](scalar x) { return x > SMALL; })); 
    volScalarField onesFilter (mapToOnes(markerField, [](scalar x) { return (x > 0) && (x < 1); })); 
    // TODO: cleanup required
    onesFilter.rename("ones"); 
    onesFilter.write(); 
    numericalCurvature *= onesFilter; 
    exactCurvature *= onesFilter; 

    Info << max(numericalCurvature).value() << " " << min(numericalCurvature).value() << endl; 

    volScalarField LinfField ("LinfCurvatureErr", mag(exactCurvature - numericalCurvature)); 
    dimensionedScalar Linf = max(LinfField);

    dimensionedScalar h = max(mag(mesh.delta())); 
    errorFile << h.value() << " " << Linf.value() << std::endl;

    // FIXME: Clean up the write calls. TM. 
    numericalCurvature.write(); 
    LinfField.write(); 
    front.write();
    runTime.write();

    Info <<"Maximal curvature error = " << Linf.value() << endl;

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
