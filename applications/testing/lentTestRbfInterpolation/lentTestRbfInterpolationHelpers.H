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

Description
    Helper functions / classes for the lentTestRbfInterpolation test application.

Author
    Tomislav Maric maric@mma.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

template<typename Points, typename Distribution, typename Generator>
void seedPoints(Points& points, Distribution& dis, Generator& gen)
{
    for (auto& point : points)
    {
        point[0] = dis(gen);  
        point[1] = dis(gen);  
        point[2] = dis(gen);  
    }
};

template<typename Values, typename Points>
void setValuesXnYnZn(Values& values, Points const& points, double N)
{
    const long unsigned int Npoints = points.size();
    const long unsigned int Nvalues = values.size(); 
    EXPECT_EQ(Npoints, Nvalues);
    for (long unsigned int i = 0; i < Npoints; ++i)
        values[i] = std::pow(points[i][0], N) 
            + std::pow(points[i][1], N) + std::pow(points[i][2], N);
};

template<typename Point, typename DIS, typename GEN>
Point rand_point(DIS& dis, GEN& gen)
{
    Point pt; 
    pt[0] = dis(gen);
    pt[1] = dis(gen);
    pt[2] = dis(gen);
    return pt;
}

class eigenPlane 
{
    Eigen::Vector3d point_; 
    Eigen::Vector3d normal_;

    public: 

        eigenPlane()
            :
                point_(), 
                normal_()
        {}

        const string name() const { return "plane"; }

        template<typename DIS, typename DISSIGN, typename GEN, typename SIGNS>
        void randomize(DIS& dis, DISSIGN& dissign, GEN& gen, SIGNS const& signs)
        {
            point_ = rand_point<Eigen::Vector3d>(dis, gen);
            normal_ = rand_point<Eigen::Vector3d>(dis, gen);
            normal_[0] = signs[dissign(gen)] * normal_[0];
            normal_[1] = signs[dissign(gen)] * normal_[1];
            normal_[2] = signs[dissign(gen)] * normal_[2];
            normal_.normalize();  
        }

        template<typename Vector>
        decltype(auto) value(Vector const& evalPoint) const
        {
            return (evalPoint - point_).dot(normal_); 
        }

        decltype(auto) normal() const
        {
            return normal_;
        }

        template<typename Vector>
        decltype(auto) grad(Vector const&) const
        {
            return normal(); 
        }

        decltype(auto) point() const
        {
            return point_;
        }
};

std::ostream& operator<<(std::ostream& os, const eigenPlane& plane)
{
    os << "PLANE, point " << plane.point() << ", normal " << plane.normal();
    return os; 
}

class eigenSphere 
{
    Eigen::Vector3d center_; 
    double          radius_;

    public: 

        eigenSphere()
            :
                center_(), 
                radius_()
        {}

        const string name() const { return "sphere"; }

        eigenSphere(Eigen::Vector3d center, double radius)
            :
                center_(center), 
                radius_(radius)
        {}

        eigenSphere(double radius, Eigen::Vector3d center)
            : 
                eigenSphere(center, radius)
        {}

        template<typename Vector>
        decltype(auto) value(Vector const& evalPoint) const
        {
            double d0 = evalPoint[0] - center_[0]; 
            double d1 = evalPoint[1] - center_[1]; 
            double d2 = evalPoint[2] - center_[2]; 

            return d0*d0 + d1*d1 + d2*d2 - radius_*radius_;  
        }

        template<typename Vector>
        decltype(auto) grad(Vector const& evalPoint) const
        {
            Vector result; 
            result[0] = 2*(evalPoint[0] - center_[0]); 
            result[1] = 2*(evalPoint[1] - center_[1]); 
            result[2] = 2*(evalPoint[2] - center_[2]); 
            return result; 
        }

        decltype(auto) center() const
        {
            return center_;
        }

        decltype(auto) radius() const
        {
            return radius_;
        }
};

std::ostream& operator<<(std::ostream& os, const eigenSphere& sphere)
{
    os << "SPHERE, center" << sphere.center() << ", radius " << sphere.radius();
    return os; 
}

template<typename Values, typename Points, typename Surface> 
void set_surface_values(
    Values& vals, 
    Points const& points, 
    Surface const& surf 
)
{
    for (decltype(points.size()) pI = 0; pI < points.size(); ++pI)
        vals[pI] = surf.value(points[pI]);
}

// ************************************************************************* //
