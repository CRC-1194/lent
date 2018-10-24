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

#ifndef rbFunctions_H
#define rbFunctions_H

#include <functional>
#include <map>
#include <string>
#include <cmath>

namespace RBF
{
    using rbfFunctionType = std::function<double(double, double, double)>;
    using rbfTableType = std::map<std::string, rbfFunctionType>;

    static const rbfTableType rbfKernels = 
    {
        {
            "WENDLAND_C0",
            [](double r, double rs=1.0, double n = 1.) 
            { 
                return std::pow(1. - (r / rs),2); 
            } 
        },
        {
            "WENDLAND_C2",
            [](double r, double rs=1.0, double n = 1.) 
            { 
                const auto eta = r / rs;
                return std::pow(1. - eta, 4) * (4.*eta + 1.); 
            }
        },
        {
            "WENDLAND_C4", 
            [](double r, double rs=1.0, double n = 1.) 
            {
                const auto eta = r / rs; 
                return std::pow(1. - eta, 6) * (35.*std::pow(eta, 2)  + 18.*eta + 3.);
            }
        },
        {
            "WENDLAND_C6",
            [](double r, double rs=1.0, double n = 1.) 
            {
                const auto eta = r / rs;
                return pow((1. - eta), 8) * (32.*std::pow(eta, 3) + 25.*std::pow(eta, 2) + 8.*(eta) +  1.);
            }

        },
        {
            "LINEAR",
            [](double r, double rs=1.0, double n = 1.) 
            {
                return r / rs;
            }
        },
        {
            "CUBIC_SPLINE",
            [](double r, double rs=1.0, double n = 1.) 
            {
                return std::pow(r / rs, 3);
            }
        },
        {
            "LINEAR_CUBIC_SPLINE",
            [](double r, double rs=1.0, double n = 1.) 
            {
                const auto eta = r / rs;
                return eta - std::pow(eta, 3);
            }
        },
        {
            "TRUNCATED_POWER",
            [](double r, double rs=1.0, double n = 1.) 
            {
                return std::pow(1. - (r / rs), n);
            }
        },
        {
            "MULTIQUADRIC",
            [](double r, double rs=1.0, double n = 1.) 
            {
                return std::sqrt(1. + std::pow(r / rs, 2));
            }
        },
        {
            "INVERSE_MULTIQUADRIC", 
            [](double r, double rs=1.0, double n = 1.) 
            {
                return std::pow(1 + (r / rs), -n);
            }
        },
        {
            "INVERSE_QUADRATIC",
            [](double r, double rs=1.0, double n = 1.) 
            {
                return 1.0  / (1.0 + std::pow(r / rs, 2));
            }
        }, 
        {
            "THIN_PLATE_SPLINE", 
            [](double r, double rs=1.0, double n = 1.) 
            {
                const auto eta = r / rs; 
                return eta * eta * std::log(eta + std::numeric_limits<double>::epsilon());
            }
        },
        {
            "GAUSSIAN", 
            [](double r, double rs=1.0, double n = 1.) 
            {
                const auto eta = r / rs; 
                return std::exp(-n * n * eta * eta);
            }
        }
    };
}

// ************************************************************************* //

#endif
