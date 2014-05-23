/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    lentAdvectRefined 

Description
    Interface advection with the LENT method coupled with local AMR in OpenFOAM. 
    
Authors
    Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de
    Mathematical Modeling and Analysis
    Center of Smart Interfaces, TU Darmstadt

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "interfaceProperties.H"
#include "incompressibleTwoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "lentMethod.H"

// Time Measurements
#include <chrono>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>
#include <unordered_map>

using namespace FrontTracking;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class timer 
{
    public: 

        typedef std::list<double> TimesContainer; 
        typedef std::unordered_map<std::string, TimesContainer> TimesMap;
        typedef TimesMap::const_iterator const_iterator;
        typedef std::list<std::string> NamesContainer;
        typedef std::chrono::high_resolution_clock Clock;
        typedef std::chrono::milliseconds milliseconds; 
        
        static const std::string totalTimeName() 
        {
            static std::string name ("Total Time"); 
            return name;  
        }

        timer() 
        {
            addName(totalTimeName()); 
        } 

        timer(const std::initializer_list<std::string>& i)
        {
            for (const auto& name : i)
            {
                times_[name].push_back(0); 
                names_.push_back(name); 
            }

            times_[totalTimeName()].push_back(0); 
            names_.push_back(totalTimeName()); 
        }

        template<typename Ostream>
        Ostream& writeHeader(Ostream& os) const
        {
            os << "#| ";
            for (auto const & name : names_)
            {
                os << name << " | ";
            }
            os << "\n"; 

            return os; 
        }

        template<typename Ostream>
        Ostream& writeTimes(Ostream& os) const 
        {
            for (auto const & name : names_)
            {
                auto timesIt = times_.find(name); 

                if (timesIt != times_.end())
                {
                    const auto& times = timesIt->second; 

                    os << times.back() << " "; 
                }
            }
            os << "\n"; 

            return os; 
        }

        template<typename Ostream>
        Ostream& writeAveragedTimes(Ostream& os) const
        {
            TimesContainer averageTimes; 

            for (const auto& name : names_)
            {
                auto timesIt = times_.find(name); 

                if (timesIt != times_.end())
                {
                    const auto& times = timesIt->second;

                    auto result = std::accumulate(
                        times.begin(), 
                        times.end(), 
                        TimesContainer::value_type(0)
                    ); 

                    result /= times.size(); 

                    averageTimes.push_back(result); 
                }
            }

            auto totalAveragedTime = std::accumulate(
                std::next(averageTimes.begin()), 
                averageTimes.end(), 
                TimesContainer::value_type(0)
            );

            auto it = averageTimes.begin(); 

            *it = totalAveragedTime; 

            for (const auto& time : averageTimes)
            {
                os << time << " "; 
            }

            return os;; 
        }

        void addName (const std::string& name)
        {
            auto it = find(names_.begin(), names_.end(), name); 

            if (it == names_.end())
            {
                names_.push_back(name); 
            }
        }

        void start(const std::string& name)
        {
            addName(name); 

            firstTime_ = Clock::now(); 
        }

        void stop(const std::string& name)
        {
            secondTime_ = Clock::now(); 

            auto diff = std::chrono::duration_cast<milliseconds>(secondTime_ - firstTime_); 
            auto diffSeconds = diff.count() / 1e03;

            auto& totalTimes = times_[totalTimeName()]; 
            auto lastTotalTime = totalTimes.back(); 
            totalTimes.push_back(lastTotalTime + diffSeconds); 

            auto& times = times_[name]; 
            times.push_back(diffSeconds); 
        } 

        const_iterator begin() const
        {
            return times_.begin(); 
        }

        const_iterator end() const
        {
            return times_.end(); 
        }


    private: 
    
        TimesMap times_;  
        NamesContainer names_;

        decltype(Clock::now()) firstTime_; 
        decltype(Clock::now()) secondTime_; 
};

template<typename OStream>
OStream& operator << (OStream& os, const timer& t)
{
    t.writeTimes(os); 
    return os; 
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    timer timing; 

    std::ofstream timingFile;
    timingFile.open("timing.dat"); 

    std::ofstream timingAveragedFile; 
    timingAveragedFile.open("timingAveraged.dat"); 

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    triSurfaceFront front(
        IOobject(
            "front.stl",
            "front",
            runTime, 
            IOobject::MUST_READ, 
            IOobject::AUTO_WRITE
        )
    );

    triSurfaceFrontGeoMesh frontMesh(front); 

    triSurfaceFrontVectorField frontVelocity(
        IOobject(
            "frontVelocity", 
            runTime.timeName(), 
            runTime, 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        front, 
        dimensionedVector(
            "zero",
            dimLength / dimTime, 
            vector(0,0,0)
        )
    ); 

    lentMethod lent(front, mesh); 

    lent.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);

    lent.reconstructFront(front, signedDistance, pointSignedDistance); 

    front.write(); 

    while (runTime.run())
    {
        #include "readTimeControls.H"

        runTime++;

        #include "CourantNo.H"
        #include "markerFieldCourantNo.H"
        #include "setDeltaT.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        timing.start("Mesh Update"); 
        mesh.update();
        timing.stop("Mesh Update");

        timing.start("Search Distance"); 
        lent.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);
        timing.stop("Search Distance"); 

        timing.start("Signed Distance"); 
        lent.calcSignedDistances(
            signedDistance, 
            pointSignedDistance,
            searchDistanceSqr,
            pointSearchDistanceSqr, 
            front
        ); 
        timing.stop("Signed Distance"); 

        timing.start("MarkerField"); 
        lent.calcMarkerField(markerField, signedDistance, searchDistanceSqr); 
        timing.stop("MarkerField"); 

        // FIXME: heisenbug in Debug mode: field checking probably fails TM, Mar 05 14
        //twoPhaseProperties.correct();

        timing.start("Reconstruction"); 
        lent.reconstructFront(front, signedDistance, pointSignedDistance); 
        timing.stop("Reconstruction"); 

        timing.start("Velocity Calculation"); 
        lent.calcFrontVelocity(frontVelocity, U); 
        timing.stop("Velocity Calculation"); 

        timing.start("Front Evolution"); 
        lent.evolveFront(front, frontVelocity); 
        timing.stop("Front Evolution"); 
        
        timing.start("Writing"); 
        runTime.write();
        timing.stop("Writing"); 

        if (runTime.timeIndex() == 1)
        {
            timing.writeHeader(timingFile); 
        }

        timing.writeTimes(timingFile); 

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    timing.writeHeader(timingAveragedFile); 
    timing.writeAveragedTimes(timingAveragedFile); 

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
