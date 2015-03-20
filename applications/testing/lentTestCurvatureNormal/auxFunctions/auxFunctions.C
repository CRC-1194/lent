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

labelList orderVertices(labelledTri& tri, label V)
{
    labelList vertices(3);

    if (tri[0] == V)
    {
        vertices[0] = tri[0];
        vertices[1] = tri[1];
        vertices[2] = tri[2];
    }
    else if (tri[1] == V)
    {
        vertices[0] = tri[1];
        vertices[1] = tri[0];
        vertices[2] = tri[2];
    }
    else
    {
        vertices[0] = tri[2];
        vertices[1] = tri[0];
        vertices[2] = tri[1];
    }

    return vertices;
}
