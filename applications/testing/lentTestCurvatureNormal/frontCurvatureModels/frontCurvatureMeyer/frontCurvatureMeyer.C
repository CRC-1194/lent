#include "frontCurvatureMeyer.H"

// Destructors
frontCurvatureMeyer::~frontCurvatureMeyer() {};

// Private member functions
frontCurvatureMeyer::calcCurvatureNormal()
{
    curvatureNormals = dimensionedVector("zero",
                            dimless/dimLength,
                            vector(0,0,0)
                          );

    // TODO: modify so that mutliple patches, e.g. in the case of
    // bubble breakup, are supported

    // All of the following functions are taken from
    // "PrimitivePatch.H"
    // Get label list of all mesh vertices
    const labelList& vertices = front.meshPoints();

    // Get assignment point --> faces
    const labelListList& adjacentFaces = front.pointFaces();

    // List of local point references for current patch
    const pointField& localPoints = front.localPoints();
    
    // Needed: List of local face references
    const List<labelledTri>& localFaces = front.localFaces();

    // V (the actual vertex) and Vl (its label) are used synonymously
    // in the following comments
    forAll(vertices, Vl)
    {
        scalar Amix = 0;

        // Get all triangles adjacent to V
        const labelList& oneRingNeighborhood = adjacentFaces[Vl];

        // Sum up mixed area and contributions to mean curvature normal
        forAll(oneRingNeighborhood, T)
        {
            // Compute area contributions
            // First, resolve label T to get the actual face
            label triLabel = oneRingNeighborhood[T];
            labelledTri currentTri = localFaces[triLabel];

            // Read point labels from triangle and determine which
            // one matches Vl to resolve point coordinates correctly
            label tri0 = currentTri[0];
            label tri1 = currentTri[1];
            label tri2 = currentTri[2];

            point V;
            point Q;
            point R;

            if (tri0 == Vl)
            {
                V = localPoints[Vl];
                Q = localPoints[tri1];
                R = localPoints[tri2];
            }
            else if (tri1 == Vl)
            {
                V = localPoints[Vl];
                Q = localPoints[tri0];
                R = localPoints[tri2];
            }
            else
            {
                V = localPoints[Vl];
                Q = localPoints[tri0];
                R = localPoints[tri1];
            }

            // Edge vectors
            vector VQ = Q - V;
            vector VR = R - V;
            vector QR = R - Q;

            // Get angles
            scalar Va = getAngle(VQ, VR);
            scalar Ra = getAngle(VR, QR);
            scalar Qa = pi - (Va + Ra);

            // Check if non-obtuse in order to use the correct area metric
            if (Va < pi/2 && Qa < pi/2 && Ra < pi/2)
            {
                // Use Voronoi-area
                // Cotangent function 'cot' has to be defined locally since
                // it is not offered by OpenFOAM
                Amix += 0.125*(magSqr(VR)*cot(Qa) + magSqr(VQ)*cot(Ra));    
            }
            else if (Va >= pi/2)
            {
                // Obtuse angle at V, use half area
                // Use Heron's formula for now until problem
                // with vector approach is solved
                Amix += 0.5 * areaHeron(VQ, VR, QR);
            }
            else
            {
                // Obtuse angle not at V, use quarter area
                // Use Heron's formula for now until problem
                // with vector approach is solved
                Amix += 0.25 * areaHeron(VQ, VR, QR);
            }

            // Compute mean curvature normal contributions
            curvatureNormals[Vl] += cot(Qa)*VR + cot(Ra)*VQ;
        }
        curvatureNormals[Vl] = curvatureNormals[Vl] / (2.0 * Amix);
    }
}
