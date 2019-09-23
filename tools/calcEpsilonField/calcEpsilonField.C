/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    calcEpsilonField

Description
    Calculate and create the SGS epsilon field as well mean resolved epsilon field from a SOWFA LES simulation.

Source files:
    readFields.H

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info << "Time = " << runTime.timeName() << endl;

        #include "readFields.H"

        if (!IOobject("epsilonSGS", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nRetrieving SGS dissipation rate field, epsilonSGS..." << endl;
            volScalarField epsilonSGS
            (
                IOobject
                (
                    "epsilonSGS",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                turbulence->epsilon()
            );

            Info << "\nWriting SGS dissipation rate field, epsilonSGS..." << endl;
            epsilonSGS.write();
        }
        else
        {
            Info << "\nSGS dissipation rate field epsilonSGS already exists!" << endl;
        }


        // If UAvg exists, calculate mean resolved dissipation rate epsilonResolved
        if (IOobject("UAvg", runTime.timeName(), mesh).headerOk())
        {
            // If epsilonResolved is not calculated already, calculate it
            if (!IOobject("epsilonResolved", runTime.timeName(), mesh).headerOk())
            {
                Info<< "Reading mean velocity field, UAvg\n" << endl;
                volVectorField UAvg
                (
                    IOobject
                    (
                        "UAvg",
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                );

                Info << "\nRetrieving mean resolved dissipation rate field, epsilonResolved..." << endl;
                volScalarField epsilonResolved
                (
                    IOobject
                    (
                        "epsilonResolved",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    2.*turbulence->nu()*(symm(fvc::grad(UAvg)) && symm(fvc::grad(UAvg)))
                );

                Info << "\nWriting mean resolved dissipation rate field, epsilonResolved..." << endl;
                epsilonResolved.write();
            }
            else
            {
                Info << "\nMean resolved dissipation rate field epsilonResolved already exists!" << endl;
            }

        }
        else
        {
            Info << "\nMean velocity field UAvg not found! Not calculating mean resolved dissipation rate epsilonResolved!" << endl;
        }

    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
