include(FetchContent)

cmake_minimum_required(VERSION 3.14)

project(lent VERSION 1.0
    DESCRIPTION "Simulation of two-phase flows on unstructured meshes using a hybrid level set / front tracking method"
        LANGUAGES CXX)
enable_testing()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Use lib instead of lib64 for libraries.
set(CMAKE_INSTALL_LIBDIR "lib")

# Set the build type suffix based on the build type.
set(PROJECT_BUILD_TYPE "")
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    message ("Using 'Release' CMAKE_BUILD_TYPE.")
elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    message ("Using 'Debug' CMAKE_BUILD_TYPE.")
    set(PROJECT_BUILD_TYPE "Debug")
elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
    message ("Using 'RelWithDebInfo' CMAKE_BUILD_TYPE.")
    set(PROJECT_BUILD_TYPE "RelWithDebInfo")
else()
    message (FATAL_ERROR "Use '-DCMAKE_BUILD_TYPE=' to define a build type. 'Release', 'Debug' and 'RelWithDebInfo' are supported.")
    return(1)
endif()

#------------------------------------------------------------------------------
#   Dependencies
#------------------------------------------------------------------------------
find_package(OpenMP REQUIRED)
set(CMAKE_PREFIX_PATH "$ENV{EIGEN3_ROOT}")
find_package(Eigen3 REQUIRED)

# Two-Phase Validation
FetchContent_Declare(two-phase-validation
    GIT_REPOSITORY          https://www.gitlab.com/leia-methods/two-phase-validation.git 
    GIT_TAG                 master
    GIT_SHALLOW             YES
    GIT_PROGRESS            YES
    USES_TERMINAL_DOWNLOAD  YES
    USES_TERMINAL_UPDATE    YES
)
FetchContent_MakeAvailable(two-phase-validation)


# OpenFOAM
# Read build relevant OpenFOAM environment variables and cache them
set(of_wm_project_dir "$ENV{WM_PROJECT_DIR}" CACHE PATH "Path to OpenFOAM project folder.")
set(of_wm_arch "$ENV{WM_ARCH}" CACHE STRING "Architecture. Usually linux64.")
set(of_wm_arch_option "$ENV{WM_ARCH_OPTION}" CACHE STRING "Information if 32 or 64 bit operating system.")
set(of_wm_precision_option "$ENV{WM_PRECISION_OPTION}" CACHE STRING "Flag if to use single precision (SP) or double precision (DP).")
set(of_wm_label_size "$ENV{WM_LABEL_SIZE}" CACHE STRING "Size in bit to use as label type. Can be either 32 or 64.")
set(of_wm_compile_option "$ENV{WM_COMPILE_OPTION}" CACHE STRING "OpenFOAM build type: Opt, Debug, Prof.")
set(of_wm_compiler "$ENV{WM_COMPILER}" CACHE STRING "Compiler used for OpenFOAM build.")
set(of_wm_label_option "$ENV{WM_LABEL_OPTION}" CACHE STRING "Concrete Type used for label. Either Int32 or Int64.")
# WM_ARCH + WM_COMPILER + WM_PRECISION_OPTION + WM_LABEL_OPTION + WM_COMPILE_OPTION
set(of_wm_options "${of_wm_arch}${of_wm_compiler}${of_wm_precision_option}${of_wm_label_option}" CACHE STRING "Name of subfolder which contains compiled OpenFOAM libraries" FORCE)

# Determine whether specific build exists
# TODO: change build path according to build option passed to cmake command (TT)
if(IS_DIRECTORY "${of_wm_project_dir}/platforms/${of_wm_options}${of_wm_compile_option}")
    set(of_wm_options "${of_wm_options}${of_wm_compile_option}")
else()
    message(FATAL_ERROR "Path ${of_wm_project_dir}/platforms/${of_wm_options}${of_wm_compile_option} does not exist. Make sure OpenFOAM environment is set.")
endif()

set(of_lib_path "${of_wm_project_dir}/platforms/${of_wm_options}/lib" CACHE PATH "Path to compiled OpenFOAM libraries.")
set(of_src_path "${of_wm_project_dir}/src" CACHE PATH "Path to OpenFOAM/src folder")

# Use lib instead of lib64 for libraries.
set(CMAKE_INSTALL_LIBDIR "$ENV{FOAM_USER_LIBBIN}")
set(CMAKE_INSTALL_BINDIR "$ENV{FOAM_USER_APPBIN}")

message(STATUS "OpenFOAM lib path: ${of_lib_path}")
message(STATUS "OpenFOAM src path: ${of_src_path}")

# Required OpenFOAM libraries
find_library(OF_FINITEVOLUME finiteVolume PATHS ${of_lib_path})
find_library(OF_MESHTOOLS meshTools PATHS ${of_lib_path})
find_library(OF_OPENFOAM OpenFOAM PATHS ${of_lib_path})
find_library(OF_SAMPLING sampling PATHS ${of_lib_path})
find_library(OF_SURFMESH surfMesh PATHS ${of_lib_path})

# OpenFoam dependencies of the lent library
set(OF_LENTDEPS ${OF_FINITEVOLUME} ${OF_MESHTOOLS} ${OF_SAMPLING} ${OF_SURFMESH} ${OF_OPENFOAM} dl m)

# Solver dependencies
find_library(OF_TWOPHASEMIXTURE twoPhaseMixture PATHS ${of_lib_path})
find_library(OF_TWOPHASEPROPERTIES twoPhaseProperties PATHS ${of_lib_path})
find_library(OF_IMMISCIBLEINCOMPRESSIBLETWOPHASEMIXTURE immiscibleIncompressibleTwoPhaseMixture PATHS ${of_lib_path})
find_library(OF_TURBULENCEMODELS turbulenceModels PATHS ${of_lib_path})
find_library(OF_INCOMPRESSIBLETURBULENCEMODELS incompressibleTurbulenceModels PATHS ${of_lib_path})
find_library(OF_FVOPTIONS fvOptions PATHS ${of_lib_path})
find_library(OF_DYNAMICMESH dynamicMesh PATHS ${of_lib_path})
find_library(OF_DYNAMICFVMESH dynamicFvMesh PATHS ${of_lib_path})

set(OF_SOLVER_COMMON ${OF_TWOPHASEMIXTURE} ${OF_TWOPHASEPROPERTIES} ${OF_IMMISCIBLEINCOMPRESSIBLETWOPHASEMIXTURE}
    ${OF_TURBULENCEMODELS} ${OF_INCOMPRESSIBLETURBULENCEMODELS} ${OF_FVOPTIONS}
)

# TODO:
# Check that those variables are set. Optional: check for reasonable values (TT)

# Compile definitions required for OpenFOAM
add_compile_definitions(
    WM_LABEL_SIZE=${of_wm_label_size}
    WM_${of_wm_precision_option}
    WM_ARCH_OPTION=${of_wm_arch_option}
    ${of_wm_arch}
    OPENFOAM="$ENV{WM_PROJECT_VERSION}" # Figures out OF version. TM. 
    NoRepository
)

# Required to make linking to OpenFOAM libraries work
set(CMAKE_EXE_LINKER_FLAGS "-Xlinker --add-needed -Xlinker --no-as-needed")

if(MSVC)
    add_compile_options(/W4 /WX)
else(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb3 -O0 -Wall -Wextra -Wfloat-conversion -DDEBUG -pedantic -D_USE_MATH_DEFINES -DBOOST_MATH_INSTRUMENT")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wall -Wextra -Wfloat-conversion -pedantic -Wno-deprecated -D_USE_MATH_DEFINES")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -ggdb3 -O3 -Wall -Wextra -Wfloat-conversion -pedantic -DDEBUG -D_USE_MATH_DEFINES")
endif(MSVC)



#------------------------------------------------------------------------------
#   Library
#------------------------------------------------------------------------------
add_library("lent" SHARED
    src/cellToVertexInterpolation/interpolationCellPointTaylorExtended/makeInterpolationCellPointTaylorExtended.C
    src/cellToVertexInterpolation/interpolationCellPointTaylor/makeInterpolationCellPointTaylor.C

    src/curvatureModels/curvatureModel/curvatureModel.C
    src/curvatureModels/divergenceBasedCurvatureModels/averagingCurvatureModel/averagingCurvatureModel.C
    src/curvatureModels/divergenceBasedCurvatureModels/compactLevelSetCurvatureModel/compactLevelSetCurvatureModel.C
    src/curvatureModels/divergenceBasedCurvatureModels/divergenceBasedCurvatureModel/divergenceBasedCurvatureModel.C
    src/curvatureModels/exactCurvatureModels/analyticalSurfaceCurvatureModel/analyticalSurfaceCurvatureModel.C
    src/curvatureModels/exactCurvatureModels/constantCurvatureModel/constantCurvatureModel.C
    src/curvatureModels/exactCurvatureModels/exactCurvatureModel/exactCurvatureModel.C
    src/curvatureModels/frontBasedCurvatureModels/frontBasedCurvatureModel/frontBasedCurvatureModel.C
    src/curvatureModels/frontBasedCurvatureModels/frontTriangleCurvatureModel/frontTriangleCurvatureModel.C
    src/curvatureModels/frontBasedCurvatureModels/meyerCurvatureModel/meyerCurvatureModel.C
    src/curvatureModels/frontToMeshTransferModels/frontToMeshTransferModel/frontToMeshTransferModel.C
    src/curvatureModels/frontToMeshTransferModels/nearestTriangleVicinityTransferModel/nearestTriangleVicinityTransferModel.C
    src/curvatureModels/frontToMeshTransferModels/triangleInCellTransferModel/triangleInCellTransferModel.C
    src/curvatureModels/frontVertexNormalCalculators/edgeInversedNormalCalculator/edgeInversedNormalCalculator.C
    src/curvatureModels/frontVertexNormalCalculators/exactNormalCalculator/exactNormalCalculator.C
    src/curvatureModels/frontVertexNormalCalculators/frontVertexNormalCalculator/frontVertexNormalCalculator.C
    src/curvatureModels/frontVertexNormalCalculators/inversedDistanceNormalCalculator/inversedDistanceNormalCalculator.C
    src/curvatureModels/frontVertexNormalCalculators/parabolaFitting2DNormalCalculator/parabolaFitting2DNormalCalculator.C
    src/curvatureModels/frontVertexNormalCalculators/parabolaFittingNormalCalculator/parabolaFittingNormalCalculator.C
    src/curvatureModels/frontVertexNormalCalculators/sinusEdgeLengthNormalCalculator/sinusEdgeLengthNormalCalculator.C

    src/distanceFieldCalculators/distanceFieldCalculator/distanceFieldCalculator.C
    src/distanceFieldCalculators/optimizedOctreeDistanceCalculator/optimizedOctreeDistanceCalculator.C
    src/distanceFieldCalculators/triSurfaceMeshDistanceCalculator/triSurfaceMeshDistanceCalculator.C

    src/frontMeshSearch/frontMeshSearch.C

    src/frontMotionSolvers/eulerFrontMotionSolver/eulerFrontMotionSolver.C
    src/frontMotionSolvers/frontMotionSolver/frontMotionSolver.C
    src/frontMotionSolvers/taylorFrontMotionSolver/taylorFrontMotionSolver.C

    src/frontReconstructionModels/averageNormalAngleFrontReconstructionModel/averageNormalAngleFrontReconstructionModel.C
    src/frontReconstructionModels/frontReconstructionModel/frontReconstructionModel.C
    src/frontReconstructionModels/maxAreaFrontReconstructionModel/maxAreaFrontReconstructionModel.C
    src/frontReconstructionModels/maxNormalAngleFrontReconstructionModel/maxNormalAngleFrontReconstructionModel.C
    src/frontReconstructionModels/reconstructionHistory/reconstructionHistory.C
    src/frontReconstructionModels/timeStepFrontReconstructionModel/timeStepFrontReconstructionModel.C
    src/frontReconstructionModels/twoCriteriaFrontReconstructionModel/twoCriteriaFrontReconstructionModel.C

    src/frontReconstructors/analyticalSurfaceReconstructor/analyticalSurfaceReconstructor.C
    src/frontReconstructors/foamIsoSurfaceReconstructor/foamIsoSurfaceFrontReconstructor.C
    src/frontReconstructors/frontReconstructor/frontReconstructor.C
    src/frontReconstructors/normalConsistencyAlgorithms/analyticalSurfaceNormalConsistency/analyticalSurfaceNormalConsistency.C
    src/frontReconstructors/normalConsistencyAlgorithms/centerNormalConsistency/centerNormalConsistency.C
    src/frontReconstructors/normalConsistencyAlgorithms/convexFrontNormalConsistency/convexFrontNormalConsistency.C
    src/frontReconstructors/normalConsistencyAlgorithms/distanceNormalConsistency/distanceNormalConsistency.C
    src/frontReconstructors/normalConsistencyAlgorithms/normalConsistency/normalConsistency.C

    src/fronts/analyticalSurfaces/analyticalCircle/analyticalCircle.C
    src/fronts/analyticalSurfaces/analyticalEllipse/analyticalEllipse.C
    src/fronts/analyticalSurfaces/analyticalEllipsoid/analyticalEllipsoid.C
    src/fronts/analyticalSurfaces/analyticalPlane/analyticalPlane.C
    src/fronts/analyticalSurfaces/analyticalRandomizedCircle/analyticalRandomizedCircle.C
    src/fronts/analyticalSurfaces/analyticalRandomizedEllipse/analyticalRandomizedEllipse.C
    src/fronts/analyticalSurfaces/analyticalRandomizedEllipsoid/analyticalRandomizedEllipsoid.C
    src/fronts/analyticalSurfaces/analyticalRandomizedPlane/analyticalRandomizedPlane.C
    src/fronts/analyticalSurfaces/analyticalRandomizedSphere/analyticalRandomizedSphere.C
    src/fronts/analyticalSurfaces/analyticalSphere/analyticalSphere.C
    src/fronts/analyticalSurfaces/analyticalSurface/analyticalSurface.C
    src/fronts/analyticalSurfaces/surfaceSet/surfaceSet.C

    src/fronts/faceTriangulators/simpleTriangulator/simpleTriangulator.C

    src/fronts/frontConstructor/frontConstructor/frontConstructor.C

    src/fronts/frontMeshIntersection/cutCellVolumeCalculator/cutCellVolumeCalculator.C

    src/frontSmoother/frontSmoother.C

    src/fronts/triSurfaceFront/triSurfaceFrontFields/triSurfaceFrontFields.C
    src/fronts/triSurfaceFront/triSurfaceFront/triSurfaceFront.C
    src/fronts/triSurfaceFront/triSurfaceFront/VTK/writeVTKWithFields.C

    src/frontSurfaceTensionForceModels/csfSurfaceTensionForceModel/csfSurfaceTensionForceModel.C
    src/frontSurfaceTensionForceModels/curvatureBasedSurfaceTensionForceModel/curvatureBasedSurfaceTensionForceModel.C
    src/frontSurfaceTensionForceModels/frontSurfaceTensionForceModel/frontSurfaceTensionForceModel.C
    src/frontSurfaceTensionForceModels/tryggvasonSurfaceTensionForceModel/tryggvasonSurfaceTensionForceModel.C

    src/frontVelocityCalculators/barycentricFrontVelocityInterpolator/barycentricFrontVelocityInterpolator.C
    src/frontVelocityCalculators/frontVelocityCalculator/frontVelocityCalculator.C

    src/geometrical/markerFieldModels/barycenterPlaneMarkerFieldModel/barycenterPlaneMarkerFieldModel.C

    src/lentCommunication/lentCommunication.C

    src/lentInterpolation/dualKrigingInterpolation/dualKrigingInterpolation.C
    src/lentInterpolation/lentInterpolation.C

    src/lentMethod/lentMethod.C

    src/lentSolutionControls/lentSolutionControl/lentSolutionControl.C

    src/markerFieldModels/harmonicMarkerFieldModel/harmonicMarkerFieldModel.C
    src/markerFieldModels/markerFieldModel/markerFieldModel.C
    src/markerFieldModels/sharpMarkerFieldModel/sharpMarkerFieldModel.C
    src/markerFieldModels/singleCellHarmonicMarkerFieldModel/singleCellHarmonicMarkerFieldModel.C
    src/markerFieldModels/tetFillingLevelAlternativeMarkerFieldModel/tetFillingLevelAlternativeMarkerFieldModel.C
    src/markerFieldModels/tetFillingLevelMarkerFieldModel/tetFillingLevelMarkerFieldModel.C
    src/markerFieldModels/twoCellHarmonicMarkerFieldModel/twoCellHarmonicMarkerFieldModel.C

    src/narrowBandPropagation/naiveNarrowBandPropagation/naiveNarrowBandPropagation.C
    src/narrowBandPropagation/narrowBandPropagation/narrowBandPropagation.C

    src/test/errorMetrics/errorMetrics.C
    src/test/lentCurvatureTest/lentCurvatureTest.C
    src/test/lentDummyTest/lentDummyTest.C
    src/test/lentFrontSmoothingTest/lentFrontSmoothingTest.C
    src/test/lentMarkerFieldTest/lentGlobalMarkerFieldTest/lentGlobalMarkerFieldTest.C
    src/test/lentSubalgorithmTest/lentSubalgorithmTest.C
    src/test/lentVertexNormalCalculatorTest/lentVertexNormalCalculatorTest.C
    src/test/randomization/noiseGenerator/noiseGenerator.C
    src/test/randomization/randomNumberGenerator/randomNumberGenerator.C
)
target_include_directories("lent" PUBLIC
    src/cellToVertexInterpolation/interpolationCellPointTaylor
    src/cellToVertexInterpolation/interpolationCellPointTaylorExtended

    src/curvatureModels/CurvatureBufferLogic
    src/curvatureModels/curvatureModel
    src/curvatureModels/divergenceBasedCurvatureModels/averagingCurvatureModel
    src/curvatureModels/divergenceBasedCurvatureModels/compactLevelSetCurvatureModel
    src/curvatureModels/divergenceBasedCurvatureModels/divergenceBasedCurvatureModel
    src/curvatureModels/exactCurvatureModels/analyticalSurfaceCurvatureModel
    src/curvatureModels/exactCurvatureModels/constantCurvatureModel
    src/curvatureModels/exactCurvatureModels/exactCurvatureModel
    src/curvatureModels/frontBasedCurvatureModels/frontBasedCurvatureModel
    src/curvatureModels/frontBasedCurvatureModels/frontTriangleCurvatureModel
    src/curvatureModels/frontBasedCurvatureModels/meyerCurvatureModel
    src/curvatureModels/frontToMeshTransferModels/frontToMeshTransferModel
    src/curvatureModels/frontToMeshTransferModels/nearestTriangleVicinityTransferModel
    src/curvatureModels/frontToMeshTransferModels/triangleInCellTransferModel
    src/curvatureModels/frontVertexNormalCalculators/edgeInversedNormalCalculator
    src/curvatureModels/frontVertexNormalCalculators/exactNormalCalculator
    src/curvatureModels/frontVertexNormalCalculators/frontVertexNormalCalculator
    src/curvatureModels/frontVertexNormalCalculators/inversedDistanceNormalCalculator
    src/curvatureModels/frontVertexNormalCalculators/parabolaFitting2DNormalCalculator
    src/curvatureModels/frontVertexNormalCalculators/parabolaFittingNormalCalculator
    src/curvatureModels/frontVertexNormalCalculators/sinusEdgeLengthNormalCalculator

    src/distanceFieldCalculators/distanceFieldCalculator
    src/distanceFieldCalculators/optimizedOctreeDistanceCalculator
    src/distanceFieldCalculators/triSurfaceMeshDistanceCalculator

    src/fieldOperators

    src/frontMeshSearch

    src/frontMotionSolvers/eulerFrontMotionSolver
    src/frontMotionSolvers/frontMotionSolver
    src/frontMotionSolvers/taylorFrontMotionSolver

    src/frontReconstructionModels/averageNormalAngleFrontReconstructionModel
    src/frontReconstructionModels/frontReconstructionModel
    src/frontReconstructionModels/maxAreaFrontReconstructionModel
    src/frontReconstructionModels/maxNormalAngleFrontReconstructionModel
    src/frontReconstructionModels/reconstructionHistory
    src/frontReconstructionModels/timeStepFrontReconstructionModel
    src/frontReconstructionModels/twoCriteriaFrontReconstructionModel

    src/frontReconstructors/analyticalSurfaceReconstructor
    src/frontReconstructors/foamIsoSurfaceReconstructor
    src/frontReconstructors/frontReconstructor
    src/frontReconstructors/normalConsistencyAlgorithms/analyticalSurfaceNormalConsistency
    src/frontReconstructors/normalConsistencyAlgorithms/centerNormalConsistency
    src/frontReconstructors/normalConsistencyAlgorithms/convexFrontNormalConsistency
    src/frontReconstructors/normalConsistencyAlgorithms/distanceNormalConsistency
    src/frontReconstructors/normalConsistencyAlgorithms/normalConsistency

    src/fronts/analyticalSurfaces/analyticalCircle
    src/fronts/analyticalSurfaces/analyticalEllipse
    src/fronts/analyticalSurfaces/analyticalEllipsoid
    src/fronts/analyticalSurfaces/analyticalPlane
    src/fronts/analyticalSurfaces/analyticalRandomizedCircle
    src/fronts/analyticalSurfaces/analyticalRandomizedEllipse
    src/fronts/analyticalSurfaces/analyticalRandomizedEllipsoid
    src/fronts/analyticalSurfaces/analyticalRandomizedPlane
    src/fronts/analyticalSurfaces/analyticalRandomizedSphere
    src/fronts/analyticalSurfaces/analyticalSphere
    src/fronts/analyticalSurfaces/analyticalSurface
    src/fronts/analyticalSurfaces/surfaceSet

    src/fronts/faceTriangulators/simpleTriangulator

    src/fronts/frontConstructor/frontConstructor

    src/fronts/frontMeshIntersection/cutCellVolumeCalculator

    src/frontSmoother

    src/fronts/triSurfaceFront/triSurfaceFront
    src/fronts/triSurfaceFront/triSurfaceFrontFields

    src/frontSurfaceTensionForceModels/csfSurfaceTensionForceModel
    src/frontSurfaceTensionForceModels/curvatureBasedSurfaceTensionForceModel
    src/frontSurfaceTensionForceModels/frontSurfaceTensionForceModel
    src/frontSurfaceTensionForceModels/tryggvasonSurfaceTensionForceModel

    src/frontVelocityCalculators/barycentricFrontVelocityInterpolator
    src/frontVelocityCalculators/frontVelocityCalculator

    src/geometrical/markerFieldModels/barycenterPlaneMarkerFieldModel

    src/lentCommunication

    src/lentInterpolation

    src/lentInterpolation/dualKrigingInterpolation

    src/lentMethod

    src/lentSolutionControls/lentSolutionControl

    src/markerFieldModels/harmonicMarkerFieldModel
    src/markerFieldModels/markerFieldModel
    src/markerFieldModels/sharpMarkerFieldModel
    src/markerFieldModels/singleCellHarmonicMarkerFieldModel
    src/markerFieldModels/tetFillingLevelAlternativeMarkerFieldModel
    src/markerFieldModels/tetFillingLevelMarkerFieldModel
    src/markerFieldModels/twoCellHarmonicMarkerFieldModel

    src/narrowBandPropagation/naiveNarrowBandPropagation
    src/narrowBandPropagation/narrowBandPropagation

    src/test
    src/test/errorMetrics
    src/test/lentCurvatureTest
    src/test/lentDummyTest
    src/test/lentFrontSmoothingTest
    src/test/lentMarkerFieldTest/lentGlobalMarkerFieldTest
    src/test/lentSubalgorithmTest
    src/test/lentVertexNormalCalculatorTest
    src/test/randomization/noiseGenerator
    src/test/randomization/randomNumberGenerator
)
# Using PUBLIC here removes the need to specify again the required include directories
# for the applications (TT).
target_include_directories("lent" SYSTEM PUBLIC
    ${EIGEN3_INCLUDE_DIR}
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/sampling/lnInclude
    ${of_src_path}/surfMesh/lnInclude
    ${of_src_path}/OSspecific/POSIX/lnInclude
    ${of_src_path}/OpenFOAM/lnInclude
)

target_link_libraries("lent"
    ${OF_LENTDEPS}
)

# Install the library
install(TARGETS "lent" LIBRARY)



#------------------------------------------------------------------------------
#   Executables
#------------------------------------------------------------------------------
set(OF_SOLVER_INCLUDES 
    applications/solvers
    ${of_src_path}/transportModels/twoPhaseMixture/lnInclude
    ${of_src_path}/transportModels
    ${of_src_path}/transportModels/incompressible/lnInclude
    ${of_src_path}/transportModels/interfaceProperties/lnInclude
    ${of_src_path}/TurbulenceModels/turbulenceModels/lnInclude
    ${of_src_path}/TurbulenceModels/incompressible/lnInclude
    ${of_src_path}/transportModels/immiscibleIncompressibleTwoPhaseMixture/lnInclude
    ${of_src_path}/fvOptions/lnInclude
)

# Solvers
add_executable("lentFoam"
    applications/solvers/lentFoam/lentFoam.C
)
target_include_directories("lentFoam" SYSTEM PRIVATE
    ${OF_SOLVER_INCLUDES}
)
target_link_libraries("lentFoam"
    "lent"
    ${OF_LENTDEPS}
    ${OF_SOLVER_COMMON}
)

add_executable("lentAdvect"
    applications/solvers/lentAdvect/lentAdvect.C
)
target_include_directories("lentAdvect" SYSTEM PRIVATE
    ${OF_SOLVER_INCLUDES}
)
target_link_libraries("lentAdvect"
    "lent"
    ${OF_LENTDEPS}
    ${OF_SOLVER_COMMON}
)

add_executable("lentReconstruct"
    applications/solvers/lentReconstruct/lentReconstruct.C
)
target_include_directories("lentReconstruct" SYSTEM PRIVATE
    ${OF_SOLVER_INCLUDES}
)
target_link_libraries("lentReconstruct"
    "lent"
    ${OF_LENTDEPS}
    ${OF_SOLVER_COMMON}
)

# Utilities
add_executable("lentSetFields"
    applications/utilities/preprocessing/lentSetFields/lentSetFields.C
)
target_include_directories("lentSetFields" SYSTEM PRIVATE
    applications/utilities/preprocessing
    ${of_src_path}/fvOptions/lnInclude
)
target_link_libraries("lentSetFields"
    "lent"
    ${OF_LENTDEPS}
    ${OF_FVOPTIONS}
)

add_executable("lentCreateFront"
    applications/utilities/preprocessing/lentCreateFront/lentCreateFront.C
)
target_include_directories("lentCreateFront" SYSTEM PRIVATE
    applications/utilities/preprocessing
)
target_link_libraries("lentCreateFront"
    "lent"
    ${OF_LENTDEPS}
)

# Tests
# TODO: list of test applications is still incomplete
add_executable("lentTestAdvection"
    applications/testing/lentTestAdvection/lentTestAdvection.C
)
target_include_directories("lentTestAdvection" SYSTEM PRIVATE
    applications/solvers
    ${OF_SOLVER_INCLUDES}
)
target_link_libraries("lentTestAdvection"
    "lent"
    ${OF_LENTDEPS}
    ${OF_SOLVER_COMMON}
    twoPhaseValidation
)

add_executable("lentTestCurvature"
    applications/testing/lentTestCurvature/lentTestCurvature.C
)
target_include_directories("lentTestCurvature" SYSTEM PRIVATE
    applications/testing
)
target_link_libraries("lentTestCurvature"
    "lent"
    ${OF_LENTDEPS}
)

# Install executables
install(TARGETS
    # Solvers
    "lentFoam"
    "lentAdvect"
    "lentReconstruct"

    # Utilities
    "lentSetFields"
    "lentCreateFront"

    # Tests
    "lentTestAdvection"
    "lentTestCurvature"

    RUNTIME
)
