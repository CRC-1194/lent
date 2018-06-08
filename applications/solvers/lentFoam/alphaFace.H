// (TT) Provides functionality to compute area fractions on faces.
// This is required by "computeRhoPhi.H"
//
// TODO: Refactoring!

#include "SortableList.H"
#include <assert.h>

scalar triArea(const point& p1, const point& p2, const point& p3)
{
    return 0.5*mag(((p2 - p1)^(p3 - p1)));
}

scalar alphaTri(SortableList<scalar> d)
{
    scalar alpha = 0.0;

    d.sort();

    if (d[2] <= 0)
    {
        alpha = 0;
    }
    else if (d[1] <= 0 && d[2] > 0)
    {
        alpha = d[2]*d[2] / ((d[2] - d[0])*(d[2] - d[1]));
    }
    else if (d[0] <= 0 && d[1] > 0)
    {
        alpha = 1 - d[0]*d[0] / ((d[1] - d[0])*(d[2] - d[0]));
    }
    else
    {
        alpha = 1;
    }

    return alpha;
}

scalar computeAlphaFace(const face& iFace, const pointScalarField& distance,
                        const pointField& points)
{
    scalar alpha = 0.0;
    scalar faceArea = 0.0;
    scalar area = 0.0;
    SortableList<scalar> triDistance(3);

    for (label I = 1; I < iFace.size()-1; I++)
    {
        triDistance[0] = distance[iFace[0]];
        triDistance[1] = distance[iFace[I]];
        triDistance[2] = distance[iFace[I+1]];

        area = triArea(points[iFace[0]], points[iFace[I]], points[iFace[I+1]]);
        alpha += alphaTri(triDistance)*area;
        faceArea += area;
    }

    alpha /= faceArea;
    assert(alpha >= 0.0 && "Negative face area fraction");
    assert(alpha <= 1.0 && "Face area fraction greater than 1");
    return alpha;
}
