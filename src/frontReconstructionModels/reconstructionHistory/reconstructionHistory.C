/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
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

Class
    Foam::FrontTracking::reconstructionHistory

SourceFiles
    reconstructionHistory.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Class that records in which time steps the front is reconstructed
    and writes this to a log file

Affiliations:
    Mathematical Modeling and Analysis Institute, Mathematics Department, 
    TU Darmstadt, Germany

Funding:
    German Research Foundation (DFG) - Project-ID 265191195 - SFB 1194

    German Research Foundation (DFG) - Project-ID MA 8465/1-1, 
    Initiation of International Collaboration 
    "Hybrid Level Set / Front Tracking methods for simulating 
    multiphase flows in geometrically complex systems"

\*---------------------------------------------------------------------------*/

#include "OFstream.H"

#include "reconstructionHistory.H"

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void reconstructionHistory::writeHistory() const
{

    fileName dataFileName = time_.rootPath() + "/" + time_.globalCaseName() + "/"
        + "reconstructionHistory.dat";

    OFstream historyFile(dataFileName);

    historyFile << "# time step number | physical time | operation" << endl;

    for (unsigned int index = 0; index < timeStepNumber_.size(); ++index)
    {
        historyFile << timeStepNumber_[index] << ' ' << physicalTime_[index]
                    << ' ' << operation_[index]
                    << endl;
    }
}

void reconstructionHistory::addOperation(const word& operation)
{
    timeStepNumber_.push_back(time_.timeIndex());
    physicalTime_.push_back(time_.timeName());
    operation_.push_back(operation);

    writeHistory();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reconstructionHistory::reconstructionHistory(const Time& time)
:
    time_{time},
    timeStepNumber_{},
    physicalTime_{}
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void reconstructionHistory::frontReconstructed()
{
    addOperation("reconstruction");
}

void reconstructionHistory::frontSmoothed()
{
    void frontReconstructed();
    // Do not add entry for smoothing if front has been reconstructed
    // in the same time step. Reconstruction always implies smoothing
    if (timeStepNumber_.size() > 0 && timeStepNumber_.back() == time_.timeIndex())
    {
        // Do nothing
    }
    else
    {
        addOperation("smoothing");
    }    
}


// ************************************************************************* //

} // End namespace FrontTracking

// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //
