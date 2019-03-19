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

Description
    A global map of RBF kernels. 

\*---------------------------------------------------------------------------*/

#ifndef rbFunctions_H
#define rbFunctions_H

#include <functional>
#include <map>
#include <string>
#include <cmath>
#include <tuple>
#include <array>
#include "rbfGeometry.H"

namespace Foam { namespace RBF {

    using rbfFunctionType = std::function<double(double, double, double)>;
    using rbfTableType = std::map<std::string, rbfFunctionType>;

    template<typename Derived>
    struct rbfKernelBase
    {
        template<typename Vector> 
        inline Vector grad(Vector const& point, Vector const& pointi, const double rs = 1.0) const
        {
            Vector result(
                point[0] - pointi[0],
                point[1] - pointi[1],
                point[2] - pointi[2]
            );

            const double r = Geometry::distance(point, pointi);

            if (r < std::numeric_limits<double>::epsilon())
                result *= static_cast<Derived const*>(this)->ddr(r,rs) / 
                    (r + std::numeric_limits<double>::epsilon());
            else
            {
                result *= static_cast<Derived const*>(this)->ddr(r,rs) / r; 
            }

            return result; 
        }
    };

    // Not using OpenFOAM TypeName to avoid virtual functions for efficiency.
    // These function objects should be inlined completely. TM.
    struct WENDLAND_C0
        :
            public rbfKernelBase<WENDLAND_C0>
    {
        static const std::string name_; 
        static const std::string& name() { return name_; };  

        inline double operator()(const double r, const double rs) const
        {
            const double x = 1. - (r / rs);
            return x * x; 
        }

        inline double ddr(const double r, const double rs = 1.0) const
        {
            return - 1.0 / rs; 
        }

    };
    const std::string WENDLAND_C0::name_ = "WENDLAND_C0";
    
    struct WENDLAND_C2
        :
            public rbfKernelBase<WENDLAND_C2>
    {
        static const std::string name_; 
        static const std::string& name() { return name_; };  

        inline double operator()(const double r, const double rs=1.0) const
        { 
            const double eta = r / rs;
            return std::pow(1. - eta, 4) * (4.*eta + 1.); 
        }

        inline double ddr(const double r, const double rs = 1.0) const
        {
            double rrs = r - rs;
            double rrs3 = rrs*rrs*rrs;

            return  (20. * r * rrs3) / (rs*rs*rs*rs*rs);
        }
    };

    const std::string WENDLAND_C2::name_ = "WENDLAND_C2";
    
    struct WENDLAND_C4
        :
            public rbfKernelBase<WENDLAND_C4>
    {
        static const std::string name_; 
        static const std::string& name() { return name_; };  

        inline double operator()(const double r, const double rs=1.0) const  
        {
            const double eta = r / rs; 
            return std::pow(1. - eta, 6) * 
                (35.*std::pow(eta, 2) + 18.*eta + 3.);
        }

        inline double ddr(const double r, const double rs=1.0) const  
        {
            double rrs5 = pow(r-rs,5);
            double rs7 = rs*rs*rs*rs*rs*rs*rs;
            return 280.*r*r*rrs5/(rs7 * rs) + 56.*r*rrs5/rs7;
        }
    };

    const std::string WENDLAND_C4::name_ = "WENDLAND_C4";

    struct WENDLAND_C6
        :
            public rbfKernelBase<WENDLAND_C6>
    {
        static const std::string name_; 
        static const std::string& name() { return name_; };  

        inline double operator()(const double r, const double rs=1.0) const 
        {
            const auto eta = r / rs;
            return pow((1. - eta), 8) * (32.*std::pow(eta, 3) + 
                    25.*std::pow(eta, 2) + 8.*(eta) +  1.);
        }

        inline double ddr(const double r, const double rs=1.0) const 
        {
            double rrs7 = pow(r - rs,7);
            double rs9 = rs*rs*rs*rs*rs*rs*rs*rs*rs;
            return 352.0*r*r*r*rrs7/(rs9 * rs * rs) + 
                154.0*r*r*rrs7/(rs9 * rs) + 
                22.0*r*rrs7/rs9;
        }
    };

    const std::string WENDLAND_C6::name_ = "WENDLAND_C6";

    struct LINEAR
        :
            public rbfKernelBase<LINEAR>
    {
        static const std::string name_;  
        static const std::string& name() { return name_; };  

        inline double operator()(const double r, const double rs=1.0) const 
        {
            return r / rs;
        }

        inline double ddr(const double r, const double rs=1.0) const 
        {
            return 1.0 / rs;
        }
    };
    const std::string LINEAR::name_ = "LINEAR";

    struct CUBIC_SPLINE
        :
            public rbfKernelBase<CUBIC_SPLINE>
    {
        static const std::string name_; 
        static const std::string& name() { return name_; };  

        inline double operator()(const double r, const double rs=1.0) const
        {
            return std::pow(r / rs, 3.0);
        }

        inline double ddr(const double r, const double rs=1.0) const
        {
            return 3.0 * (r*r) / (rs*rs);
        }
    };
    const std::string CUBIC_SPLINE::name_ = "CUBIC_SPLINE";

    struct LINEAR_CUBIC_SPLINE
        :
            public rbfKernelBase<LINEAR_CUBIC_SPLINE>
    {
        static const std::string name_; 
        static const std::string& name() { return name_; };  

        inline double operator()(const double r, const double rs=1.0) const
        {
            const auto eta = r / rs;
            return eta - std::pow(eta, 3);
        }

        inline double ddr(const double r, const double rs=1.0) const
        {
            double rssqr = rs*rs; 
            return (-3*r*r + rssqr)/(rssqr * rs);
        }
    };
    const std::string LINEAR_CUBIC_SPLINE::name_ = "LINEAR_CUBIC_SPLINE";

    struct TRUNCATED_POWER
    {
        static const std::string name_; 
        static const std::string& name() { return name_; };  

        inline double operator()
        (
            const double r, 
            const double rs=1.0, 
            const double n = 1.
        ) const
        {
            return std::pow(1. - (r / rs), n);
        }
    };
    const std::string TRUNCATED_POWER::name_ = "TRUNCATED_POWER";

    struct MULTIQUADRIC
        :
            public rbfKernelBase<MULTIQUADRIC>
    {
        static const std::string name_; 
        static const std::string& name() { return name_; };  

        inline double operator()(const double r, const double rs=1.0) const
        {
            return std::sqrt(1. + std::pow(r / rs, 2));
        }

        inline double ddr(const double r, const double rs=1.0) const
        {
            double rssqr = rs*rs;
            return r / (rssqr*sqrt(r*r/rssqr + 1.0));
        }
    };
    const std::string MULTIQUADRIC::name_ = "MULTIQUADRIC";

    struct INVERSE_MULTIQUADRIC
        :
            public rbfKernelBase<INVERSE_MULTIQUADRIC>
    {
        static const std::string name_; 
        static const std::string& name() { return name_; };  
        
        inline double operator() (const double r, const double rs=1.0) const 
        {
            return 1.0 / (1 + (r / rs));
        }

        inline double ddr(const double r, const double rs=1.0) const 
        {
            double s = r+rs;
            return  -rs / (s*s);
        }
    };
    const std::string INVERSE_MULTIQUADRIC::name_ = "INVERSE_MULTIQUADRIC";

    struct INVERSE_QUADRATIC
        :
            public rbfKernelBase<INVERSE_QUADRATIC>
    {
        static const std::string name_;
        static const std::string& name() { return name_; };  

        inline double operator()(const double r, const double rs=1.0) const
        {
            double eta = r / rs;
            return 1.0  / (1.0 + eta*eta);
        }

        inline double ddr(const double r, const double rs=1.0) const
        {
            double rssqr = rs*rs;
            return -2*r*rssqr/((r*r + rssqr)*(r*r + rssqr));
        }
    }; 
    const std::string INVERSE_QUADRATIC::name_ = "INVERSE_QUADRATIC"; 

    struct THIN_PLATE_SPLINE
        :
            public rbfKernelBase<THIN_PLATE_SPLINE>
    {
        static const std::string name_; 
        static const std::string& name() { return name_; };  

        inline double operator()(const double r, const double rs=1.0) const 
        {
            const auto eta = r / rs; 
            return eta * eta * 
                std::log(eta + std::numeric_limits<double>::epsilon());
        }

        inline double ddr(const double r, const double rs=1.0) const 
        {
            const double eps = std::numeric_limits<double>::epsilon();
            return r*(r + 2*(eps*rs + r)*log((eps*rs + r)/rs))/(rs*rs*(eps*rs + r));
        }
    };
    const std::string THIN_PLATE_SPLINE::name_ = "THIN_PLATE_SPLINE"; 

    struct GAUSSIAN
        :
            public rbfKernelBase<GAUSSIAN>
    {
        static const std::string name_; 
        static const std::string& name() { return name_; };  

        inline double operator()(const double r, const double rs=1.0) const 
        {
            const auto eta = r / rs; 
            return std::exp(-eta * eta);
        }

        inline double ddr(const double r, const double rs=1.0) const 
        {
            double rssqr = rs*rs;  
            return -2.*r*exp(-(r*r)/rssqr)/rssqr;
        }
    };
    const std::string GAUSSIAN::name_ = "GAUSSIAN"; 

    using wendlandRbfTuple = std::tuple
    <
        WENDLAND_C0,
        WENDLAND_C2,
        WENDLAND_C4,
        WENDLAND_C6
    >;

    using rbfTupleAll = std::tuple
    <
        WENDLAND_C0,
        WENDLAND_C2,
        WENDLAND_C4,
        WENDLAND_C6,
        LINEAR,
        CUBIC_SPLINE,
        LINEAR_CUBIC_SPLINE,
        MULTIQUADRIC,
        INVERSE_MULTIQUADRIC,
        INVERSE_QUADRATIC,
        THIN_PLATE_SPLINE,
        GAUSSIAN
    >;

    using rbfTuple = std::tuple
    <
        WENDLAND_C2,
        WENDLAND_C4,
        LINEAR,
        CUBIC_SPLINE,
        LINEAR_CUBIC_SPLINE,
        MULTIQUADRIC,
        INVERSE_QUADRATIC,
        THIN_PLATE_SPLINE,
        GAUSSIAN
    >;

}} // End namespace Foam::RBF 

// ************************************************************************* //

#endif
