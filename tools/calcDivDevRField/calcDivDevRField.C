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
    calcDivRField

Description
    Calculate and create the turbulent shear stress gradient
    -div(ui'uj') field. Can detect SOWFA LES and RANS case automatically.

Source files:
    readFields.H

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
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

        // If uuPrime2 exists, it's from SOWFA LES and calculate mean deviatoric total Ryenolds stress divergence vector field divDevR
        if (IOobject("uuPrime2", runTime.timeName(), mesh).headerOk())
        {
            // If mean deviatoric total Ryenolds stress divergence vector field divDevR has not been calculated already, calculate it
            if (!IOobject("divDevR", runTime.timeName(), mesh).headerOk())
            {
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

                Info << "\nComputing mean deviatoric SGS Reynolds stress field, devRSGS..." << endl;
                volSymmTensorField devRSGS
                (
                    IOobject
                    (
                        "devRSGS",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    -nuSGSmean*twoSymm(fvc::grad(UAvg))
                );

                Info << "\nComputing mean deviatoric total Reynolds stress divergence vector field divDevR..." << endl;
                volVectorField divDevR
                (
                    IOobject
                    (
                        "divDevR",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    -fvc::div(dev(uuPrime2) + devRSGS)  // devRResolved is RResolved - 2/3*kResolved*I
                );

                Info << "\nWriting mean deviatoric total Reynolds stress divergence vector field, divDevR..." << endl;
                divDevR.write();

                scalar divDevR_min = min(cmptMin(divDevR));
                scalar divDevR_max = max(cmptMax(divDevR));
                Info << "Min divDevR is " << divDevR_min << "; max is " << divDevR_max << endl;
            }
            else
            {
                Info << "\nMean deviatoric total Reynolds stress vector field divDevR already exists!" << endl;
            }

            // If predicted anisotropy tensor bij_pred exists, it's from SOWFA LES ML
            // and calculate mean predicted deviatoric Ryenolds stress divergence vector field divDevR_pred
            if (IOobject("bij_pred", runTime.timeName(), mesh).headerOk())
            {
                // If mean predicted deviatoric Ryenolds stress divergence vector field divDevR_pred has not been calculated already, calculate it
                if (!IOobject("divDevR_pred", runTime.timeName(), mesh).headerOk())
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

                    Info<< "Computing mean predicted deviatoric ui'uj' tensor field devR_pred\n" << endl;
                    volSymmTensorField devR_pred
                    (
                        IOobject
                        (
                            "devR_pred",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        2.*(kResolved + kSGSmean)*bij_pred // ui'uj'^D = 1/3*2k*I + 2k*bij - 1/3*2k*I
                    );

                    Info << "\nComputing mean predicted deviatoric Reynolds stress divergence vector field divDevR_pred..." << endl;
                    volVectorField divDevR_pred
                    (
                        IOobject
                        (
                            "divDevR_pred",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        -fvc::div(devR_pred)
                    );

                    Info << "\nWriting mean predicted deviatoric Reynolds stress divergence vector field, divDevR_pred..." << endl;
                    divDevR_pred.write();
                }
                else
                {
                    Info << "\nMean predicted deviatoric Reynolds stress vector field divDevR_pred already exists!" << endl;
                }
            }
        }
        // Otherwise it's SOWFA RANS
        else
        {
            if (!IOobject("divDevR_blend", runTime.timeName(), mesh).headerOk())
            {
                // Read turbulent kinetic energy field
                Info<< "Reading (LES/ML blended) Reynolds stress field Rij..." << endl;
                volSymmTensorField Rij
                (
                    IOobject
                    (
                        "Rij",
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                );

                // Info << "\nComputing Reynolds-averaged deviatoric Reynolds stress divergence vector field divDevR..." << endl;
                // volVectorField divDevR
                // (
                //     IOobject
                //     (
                //         "divDevR",
                //         runTime.timeName(),
                //         mesh,
                //         IOobject::NO_READ,
                //         IOobject::AUTO_WRITE
                //     ),
                //     -fvc::div(-nut*twoSymm(fvc::grad(U)))
                // );

                // Info << "\nWriting Reynolds-averaged deviatoric Reynolds stress divergence vector field divDevR..." << endl;
                // divDevR.write();

                // scalar divDevR_min = min(cmptMin(divDevR));
                // scalar divDevR_max = max(cmptMax(divDevR));
                // Info << "Min divDevR is " << divDevR_min << "; max is " << divDevR_max << endl;

                // From here blended divDevR is calculated
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

                Info << "\nComputing blended Reynolds-avgeraged deviatoric Reynolds stress divergence vector field divDevR_blend..." << endl;
                scalar mix_ratio = 0.;
                scalar t = runTime.value();

                // Mixing ratio of LES-RANS bij, capped on request
                if (t > mix_startTime.value())
                {
                    mix_ratio = min((t - mix_startTime.value())*mix_ratio_cap.value()/mix_duration.value(),
                    mix_ratio_cap.value());
                }

                Info << "LES/ML-RANS bij mixing ratio is " << mix_ratio << endl;
                volVectorField divDevR_blend
                (
                    IOobject
                    (
                        "divDevR_blend",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    -fvc::div(dev(Rij))
                );

                Info << "\nWriting blended Reynolds-avgeraged deviatoric Reynolds stress divergence vector field divDevR_blend..." << endl;
                divDevR_blend.write();

                scalar divDevR_blend_min = min(cmptMin(divDevR_blend));
                scalar divDevR_blend_max = max(cmptMax(divDevR_blend));
                Info << "Min divDevR_blend is " << divDevR_blend_min << "; max is " << divDevR_blend_max << endl;
            }
            else
            {
                Info << "\nReynolds-avgeraged deviatoric Reynolds stress divergence vector field divDevR and divDevR_blend already exist!" << endl;
            }

            // If predicted anisotropy tensor bij_pred exists, it's from SOWFA LES ML
            // and calculate mean predicted deviatoric Ryenolds stress divergence vector field divDevR_pred
            if (IOobject("bij_pred", runTime.timeName(), mesh).headerOk())
            {
                // If mean predicted deviatoric Ryenolds stress divergence vector field divDevR_pred has not been calculated already, calculate it
                if (!IOobject("divDevR_pred", runTime.timeName(), mesh).headerOk())
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

                    Info<< "Reading mean TKE field k\n" << endl;
                    volScalarField k
                    (
                        IOobject
                        (
                            "k",
                            runTime.timeName(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh
                    );

                    Info<< "Computing mean predicted deviatoric ui'uj' tensor field devR_pred\n" << endl;
                    volSymmTensorField devR_pred
                    (
                        IOobject
                        (
                            "devR_pred",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        2.*k*bij_pred // ui'uj'^D = 1/3*2k*I + 2k*bij - 1/3*2k*I
                    );

                    Info << "\nComputing mean predicted deviatoric Reynolds stress divergence vector field divDevR_pred..." << endl;
                    volVectorField divDevR_pred
                    (
                        IOobject
                        (
                            "divDevR_pred",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        -fvc::div(devR_pred)
                    );

                    Info << "\nWriting mean predicted deviatoric Reynolds stress divergence vector field, divDevR_pred..." << endl;
                    divDevR_pred.write();
                }
                else
                {
                    Info << "\nMean predicted deviatoric Reynolds stress vector field divDevR_pred already exists!" << endl;
                }
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
