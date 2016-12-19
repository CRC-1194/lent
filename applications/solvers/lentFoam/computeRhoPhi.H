// (TT) This is part of the new approach to compute the mass fluxes rhoPhi
// at the faces. It computes the area fraction for faces which belong to an
// interface cell (meaning they may be intersected by the interface) and then
// computes the face density based on this area fraction. The area fractions
// are computed with the point signed distances of each face (see Detrixhe and
// Aslam paper)
//
// TODO: Refactor this (including the functinality from "alphaFace.H").
// Maybe this should be encapsulated in a class?

// Compute rhoPhi using the old approach since it is valid for all bulk cells
rhoPhi = fvc::interpolate(rho) * phi;

// Correct faces which are part of an interface cell (0 < alpha < 1)
const auto& faceOwner = mesh.owner();
const auto& faceNeighbour = mesh.neighbour();

const auto& meshVertices = mesh.points();
const auto& meshFaces = mesh.faces();

scalar rhoFace = 0;
scalar alphaFace = 0;
scalar alphaOwn = 0;
scalar alphaNei = 0;

forAll(rhoPhi, faceI)
{
    alphaOwn = markerField[faceOwner[faceI]];
    alphaNei = markerField[faceNeighbour[faceI]];

    if
    (
        (alphaOwn > 0.0 && alphaOwn < 1.0)
        ||
        (alphaNei > 0.0 && alphaNei < 1.0)
    )
    {
        alphaFace = computeAlphaFace(meshFaces[faceI],
                                     pointSignedDistance, meshVertices);
        rhoFace = alphaFace*rho1.value() + (1 - alphaFace)*rho2.value();

        rhoPhi[faceI] = rhoFace * phi[faceI];
    }
}
