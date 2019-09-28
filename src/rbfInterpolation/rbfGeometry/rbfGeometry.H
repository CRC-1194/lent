/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Tomislav Maric, TU Darmstadt
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

\*---------------------------------------------------------------------------*/

#ifndef rbfGeometry_H
#define rbfGeometry_H

namespace Foam 
{
namespace RBF
{

namespace Geometry
{

    template<typename PointA, typename PointB> 
    double distance_sqr(PointA const& pA, PointB const & pB)
    {
        return pow(pB[0] - pA[0], 2) + 
               pow(pB[1] - pA[1], 2) + 
               pow(pB[2] - pA[2], 2);
    }

    template<typename PointA, typename PointB> 
    double distance(PointA const& pA, PointB const & pB)
    {
        return std::sqrt(distance_sqr(pA, pB));
    }

    template<typename PointA, typename PointB> 
    PointA operator+(PointA const& p, PointB const&q)
    {
        return PointA {p[0] + q[0], p[1] + q[1], p[2] + q[2]}; 
    }

    template<typename Point> 
    void operator /=(Point& p, double d) 
    {
        p[0] = p[0] / d; 
        p[1] = p[1] / d; 
        p[2] = p[2] / d; 
    }

    template<typename Points> 
    typename Points::value_type centroid(Points const& points)
    {
        typename Points::value_type centroid {0.,0.,0.}; 

        for (const auto& point : points)
            centroid = centroid + point; 

        centroid /= points.size(); 

        return centroid; 
    }

    // - Compute the fill distance of the point set. 
    // - Used as 'h' for O(h^p) convergence tests.
    template<typename Points> 
    double fillDistance(Points const & points)
    {
        double maxi = -std::numeric_limits<double>::infinity();  
        using sizeType = decltype(points.size());  

        for (sizeType i = 0; i < points.size(); ++i)
        {
            double minj = std::numeric_limits<double>::infinity();
            for (sizeType j = 0; j < points.size(); ++j)
            {
                if (i != j)
                {
                    const double dist = distance(points[i], points[j]);
                    if (dist < minj)
                        minj = dist; 
                }
            }
            if (minj > maxi)
                maxi = minj;
        }

        return maxi; 
    }

    // Maximal radius between a point and a centroid.
    template<typename Points>
    double maxCentroidRadius(Points const& points) 
    {
        double rmax = -std::numeric_limits<double>::infinity(); 

        auto pointsCentroid = centroid(points);

        for (const auto& point : points)
        {
            const double dist = distance(pointsCentroid, point); 
            if (dist > rmax)
                rmax = dist; 
        }

        return rmax; 
    }

    // Mean point-to-centroid radius.
    template<typename Points>
    double meanCentroidRadius(Points const& points) 
    {
        double rcmean = 0.;  

        auto pointsCentroid = centroid(points);
        for (const auto& point : points)
            rcmean += distance(pointsCentroid, point); 

        rcmean /= points.size(); 

        return rcmean; 
    }


    // Root mean square radius of points.
    template<typename Points>
    double rmsRadius(Points const& points) 
    {
        double rrms = 0.;  

        for (const auto& pointI : points)
            for (const auto& pointJ : points)
            rrms += distance_sqr(pointI, pointJ); 

        rrms = std::sqrt(rrms / points.size());  

        return rrms; 
    }

} // End namespace Geometry 
} // End namespace RBF 
} // End namespace Foam  


// ************************************************************************* //

#endif
