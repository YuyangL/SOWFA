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
    calcTurbulenceProductionField

Description
    Calculate and create the turbulence production field. Can detect SOWFA LES and RANS case automatically.

Source files:
    readFields.H

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// #include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
// #include "LESModel.H"
// #include "RASModel.H"
#include "IOdictionary.H"
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

        // If UAvg exists, it's from SOWFA LES and calculate mean turbulence production
        if (IOobject("UAvg", runTime.timeName(), mesh).headerOk())
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

            Info<< "Reading mean SGS viscosity field, nuSGSmean\n" << endl;
            volScalarField nuSGSmean
            (
                IOobject
                (
                    "nuSGSmean",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            Info<< "Reading mean resolved Reynolds stress field, uuPrime2\n" << endl;
            volSymmTensorField uuPrime2
            (
                IOobject
                (
                    "uuPrime2",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            // If mean SGS turbulence production field, GSGSmean has not been calculated already, calculate it
            if (!IOobject("GSGSmean", runTime.timeName(), mesh).headerOk())
            {
                Info << "\nRetrieving mean SGS turbulence production field, GSGSmean..." << endl;
                volScalarField GSGSmean
                (
                    IOobject
                    (
                        "GSGSmean",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    2.*nuSGSmean*magSqr(symm(fvc::grad(UAvg)))
                );

                Info << "\nWriting mean SGS turbulence production field, GSGSmean..." << endl;
                GSGSmean.write();

                scalar gsgs_min = min(GSGSmean).value();
                scalar gsgs_max = max(GSGSmean).value();
                scalar gsgs_mean = GSGSmean.weightedAverage(mesh.V()).value();
                Info << "Min GSGSmean is " << gsgs_min << "; max is " << gsgs_max << "; weight mean is " << gsgs_mean << endl;
            }
            else
            {
                Info << "\nMean SGS turbulence production field, GSGSmean already exists!" << endl;
            }

            // If mean resolved turbulence production field, GResolved has not been calculated already, calculate it
            if (!IOobject("GResolved", runTime.timeName(), mesh).headerOk())
            {
                Info << "\nRetrieving mean resolved turbulence production field, GResolved..." << endl;
                volScalarField GResolved
                (
                    IOobject
                    (
                        "GResolved",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    -uuPrime2 && fvc::grad(UAvg)  // tau_ij u_i, tau_ij = -u_i'u_j'
                );

                Info << "\nWriting mean resolved turbulence production field, GResolved..." << endl;
                GResolved.write();

                scalar gresolved_min = min(GResolved).value();
                scalar gresolved_max = max(GResolved).value();
                scalar gresolved_mean = GResolved.weightedAverage(mesh.V()).value();
                Info << "Min GResolved is " << gresolved_min << "; max is " << gresolved_max << "; weight mean is " << gresolved_mean << endl;
            }
            else
            {
                Info << "\nMean resolved turbulence production field, GResolved already exists!" << endl;
            }

            // If mean total turbulence production field, GAvg has not been calculated already, calculate it
            if (!IOobject("GAvg", runTime.timeName(), mesh).headerOk())
            {
                Info << "\nRetrieving mean total turbulence production field, GAvg..." << endl;
                volScalarField GAvg
                (
                    IOobject
                    (
                        "GAvg",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    (-uuPrime2 && fvc::grad(UAvg))  // GResolved
                    + 2.*nuSGSmean*magSqr(symm(fvc::grad(UAvg)))  // GSGSmean
                );

                Info << "\nWriting mean total turbulence production field, GAvg..." << endl;
                GAvg.write();

                scalar gavg_min = min(GAvg).value();
                scalar gavg_max = max(GAvg).value();
                scalar gavg_mean = GAvg.weightedAverage(mesh.V()).value();
                Info << "Min GAvg is " << gavg_min << "; max is " << gavg_max << "; weight mean is " << gavg_mean << endl;
            }
            else
            {
                Info << "\nMean total turbulence production field, GAvg already exists!" << endl;
            }

            // If bij_pred exists, it's from SOWFA LES ML and calculate mean predicted turbulence production
            if (IOobject("bij_pred", runTime.timeName(), mesh).headerOk())
            {
                // If mean predicted turbulence production field G_pred has not been calculated already, calculate it
                if (!IOobject("G_pred", runTime.timeName(), mesh).headerOk())
                {
                    Info<< "Reading predicted anisotropy tensor field, bij_pred\n" << endl;
                    volSymmTensorField bij_pred
                    (
                        IOobject
                        (
                            "bij_pred",
                            runTime.timeName(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh
                    );

                    Info<< "Reading mean resolved TKE field kResolved\n" << endl;
                    volScalarField kResolved
                    (
                        IOobject
                        (
                            "kResolved",
                            runTime.timeName(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh
                    );

                    Info<< "Reading mean SGS TKE field kSGSmean\n" << endl;
                    volScalarField kSGSmean
                    (
                        IOobject
                        (
                            "kSGSmean",
                            runTime.timeName(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh
                    );

                    Info<< "Computing mean predicted ui'uj' tensor field uuPrime2_pred\n" << endl;
                    volSymmTensorField uuPrime2_pred(
                        IOobject(
                            "uuPrime2_pred",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE),
                        2./3.*(kResolved + kSGSmean)*I + 2.*(kResolved + kSGSmean)*bij_pred // ui'uj' = 1/3*2k*I + 2k*bij
                    );

                    Info << "\nComputing mean predicted turbulence production field G_pred..." << endl;
                    volScalarField G_pred
                    (
                        IOobject
                        (
                            "G_pred",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        -uuPrime2_pred && fvc::grad(UAvg)
                    );

                    Info << "\nWriting mean predicted turbulence production field G_pred..." << endl;
                    G_pred.write();
                }
                else
                {
                    Info << "\nMean predicted turbulence production field G_pred already exists!" << endl;
                }
            }
        }
        // Otherwise it's SOWFA RANS
        else
        {
            if (!IOobject("G", runTime.timeName(), mesh).headerOk())
            {
                // Read turbulent viscosity field
                Info<< "Reading Reynolds-averaged turbulent viscosity field, nut..." << endl;
                volScalarField nut
                (
                    IOobject
                    (
                        "nut",
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                );

                // Read the LES anisotropy tensor field from time 0
                Info << "Reading the LES anisotropy tensor field, bij from time 0..." << endl;
                const volSymmTensorField bij
                (
                    IOobject
                    (
                        "bij",
                        runTime.timeName(0),  // Read from time 0 only
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    // If not found, it should read as 0 as below
                    dimensionedSymmTensor
                    (
                        "bij",
                        dimensionSet(0, 0, 0, 0, 0, 0, 0),  // bij is dimless
                        symmTensor (0, 0, 0, 0, 0, 0)
                    )
                );

                IOdictionary transportDict
                (
                    IOobject
                    (
                        "transportProperties",
                        U.time().constant(),
                        U.db(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );

                dimensionedScalar mix_startTime
                (
                    transportDict.lookupOrDefault<dimensionedScalar>("mix_startTime", 2000)
                );
                dimensionedScalar mix_duration
                (
                    transportDict.lookupOrDefault<dimensionedScalar>("mix_duration", 2000)
                );
                dimensionedScalar mix_ratio_cap
                (
                    transportDict.lookupOrDefault<dimensionedScalar>("mix_ratio_cap", 0.5)
                );

                Info << "\nRetrieving (LES/ML blended) turbulence production field, G..." << endl;
                scalar mix_ratio = 0.;
                scalar t = runTime.value();

                // Mixing ratio of LES-RANS bij, capped on request
                if (t > mix_startTime.value())
                {
                    mix_ratio = min((t - mix_startTime.value())*mix_ratio_cap.value()/mix_duration.value(),
                    mix_ratio_cap.value());
                }

                Info << "LES/ML-RANS bij mixing ratio is " << mix_ratio << endl;
                volScalarField G
                (
                    IOobject
                    (
                        "G",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    ((1. - mix_ratio)*nut*twoSymm(fvc::grad(U))
                    + mix_ratio*(-2.*k*bij))
                    && fvc::grad(U)
                );

                Info << "\nWriting (LES/ML blended) turbulence production field, G..." << endl;
                G.write();

                scalar g_min = min(G).value();
                scalar g_max = max(G).value();
                scalar g_mean = G.weightedAverage(mesh.V()).value();
                Info << "Min G is " << g_min << "; max is " << g_max << "; weight mean is " << g_mean << endl;
            }
            else
            {
                Info << "\nReynolds-averaged turbulence production field, G already exists!" << endl;
            }

            // If bij_pred exists, it's from SOWFA RANS ML and calculate mean predicted turbulence production
            if (IOobject("bij_pred", runTime.timeName(), mesh).headerOk())
            {
                // If mean predicted turbulence production field G_pred has not been calculated already, calculate it
                if (!IOobject("G_pred", runTime.timeName(), mesh).headerOk())
                {
                    Info<< "Reading predicted anisotropy tensor field, bij_pred\n" << endl;
                    volSymmTensorField bij_pred
                    (
                        IOobject
                        (
                            "bij_pred",
                            runTime.timeName(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh
                    );

                    Info << "\nComputing mean predicted turbulence production field G_pred..." << endl;
                    volScalarField G_pred
                    (
                        IOobject
                        (
                            "G_pred",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        -(2./3.*k*I + 2.*k*bij_pred) && fvc::grad(U)
                    );

                    Info << "\nWriting mean predicted turbulence production field G_pred..." << endl;
                    G_pred.write();
                }
                else
                {
                    Info << "\nMean predicted turbulence production field G_pred already exists!" << endl;
                }
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
