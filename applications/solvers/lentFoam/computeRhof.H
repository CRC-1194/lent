// FIXME: This shouldn't be boilerplate, refactor. TM. 

// Compute rhoPhi using the old approach since it is valid for all bulk cells
rhof.oldTime() == rhof;
rhof == fvc::interpolate(rho);

// Correct faces which are part of an interface cell (0 < alpha < 1)
const auto& faceOwner = mesh.owner();
const auto& faceNeighbour = mesh.neighbour();
const auto& meshVertices = mesh.points();
const auto& meshFaces = mesh.faces();

forAll(rhoPhi, faceI)
{
    const scalar& alphaOwn = markerField[faceOwner[faceI]];
    const scalar& alphaNei = markerField[faceNeighbour[faceI]];

    if
    (
        ((alphaOwn > 0.0) && (alphaOwn < 1.0))
        ||
        ((alphaNei > 0.0) && (alphaNei < 1.0))
    )
    {
        const scalar alphaFace = computeAlphaFace(meshFaces[faceI],
                                                  pointSignedDistance, meshVertices);
        rhof[faceI] = alphaFace*rho1.value() + (1 - alphaFace)*rho2.value();
    }
}
