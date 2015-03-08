// Aux functions

#include "fvCFD.H"
#include "triSurfaceFields.H"

// Calculate angle between two vectors
scalar getAngle(vector& a, vector& b)
{
    return Foam::acos((a & b) / (mag(a) * mag(b)));
}

// Cotangent function to resemble algorihm notion from paper
scalar cot(scalar angle)
{
    return 1.0/Foam::tan(angle);
}

// Triangle area functions
// Use heron's formula
scalar areaHeron(vector& A, vector& B, vector& C)
{
    scalar a = mag(A);
    scalar b = mag(B);
    scalar c = mag(C);
    scalar s = 0.5 * (a + b + c);
                
    return Foam::sqrt(s*(s-a)*(s-b)*(s-c));
}

// Use cross product property
scalar areaCross(vector& a, vector& b)
{
    return 0.5*mag(a ^ b);
}

// Aux function to get a reference point for consistent normal vectors
// Does only work for convex, closed surfaces
point getRefPoint(const triSurface& front)
{
    point refPoint = vector(0,0,0);

    // Approach: get extremal points for each axis (x,y,z), thereby putting
    // the front in a box. The center of this cuboid should also be located
    // inside the front (assuming the surface is convex...)
    const labelList& vertices = front.meshPoints();
    const pointField& localPoints = front.localPoints();

    point initValue = localPoints[0];

    scalar xmin = initValue[0];
    scalar xmax = initValue[0];
    scalar ymin = initValue[1];
    scalar ymax = initValue[1];
    scalar zmin = initValue[2];
    scalar zmax = initValue[2];

    forAll(vertices, Vl)
    {
        point V = localPoints[Vl];

        // TODO: replace by min / max functions from standard library
        xmin = V[0] < xmin ? V[0] : xmin;
        xmax = V[0] > xmax ? V[0] : xmax;

        ymin = V[1] < ymin ? V[1] : ymin;
        ymax = V[1] > ymax ? V[1] : ymax;

        zmin = V[2] < zmin ? V[2] : zmin;
        zmax = V[2] > zmax ? V[2] : zmax;
    }

    // Assemble reference point
    refPoint[0] = 0.5 * (xmin+xmax);
    refPoint[1] = 0.5 * (ymin+ymax);
    refPoint[2] = 0.5 * (zmin+zmax);

    return refPoint;
}
