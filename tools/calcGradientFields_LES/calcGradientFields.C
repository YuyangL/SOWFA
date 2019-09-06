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
    createGradientFields

Description
    Creates a gradient fields of UAvg, kResolved, kSGSmean, p_rghAvg, TAvg, which are relevant to ML.

Source files:
    readFields.H

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// #include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
// #include "RASModel.H"

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

        //Info<< "\nRetrieving field k from turbulence model" << endl;
        //const volScalarField k(turbulence->k());

        //Info<< "\nRetrieving field epsilon from turbulence model" << endl;
        //const volScalarField epsilon(turbulence->epsilon());

        // Info<< "\nRetrieving field R from turbulence model" << endl;
        // const volSymmTensorField R(turbulence->R());

        Info << "\nRetrieving mean resolved velocity gradient field grad_UAvg" << endl;
        volTensorField grad_UAvg
        (
            IOobject
            (
                "grad_UAvg",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(UAvg)
        );

        Info << "\nRetrieving mean resolved TKE gradient field grad_kResolved" << endl;
        volVectorField grad_kResolved
        (
            IOobject
            (
                "grad_kResolved",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(kResolved)
        );

        Info << "\nRetrieving mean SGS TKE gradient field grad_kSGSmean" << endl;
        volVectorField grad_kSGSmean
        (
            IOobject
            (
                "grad_kSGSmean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(kSGSmean)
        );

        Info << "\nRetrieving mean modified pressure gradient field grad_p_rghAvg" << endl;
        volVectorField grad_p_rghAvg
        (
            IOobject
            (
                "grad_p_rghAvg",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(p_rghAvg)
        );

        Info << "\nRetrieving mean resolved temperature gradient field grad_TAvg" << endl;
        volVectorField grad_TAvg
        (
            IOobject
            (
                "grad_TAvg",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(TAvg)
        );

        if (!IOobject("grad_UAvg", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nWriting mean resolved velocity gradient field grad_UAvg" << endl;
            grad_UAvg.write();
        }
        else
        {
            Info << "\nMean resolved velocity gradient field grad_UAvg already exists" << endl;
        }

        if (!IOobject("grad_kResolved", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nWriting mean resolved TKE gradient field grad_kResolved" << endl;
            grad_kResolved.write();
        }
        else
        {
            Info << "\nMean resolved TKE gradient field grad_kResolved already exists" << endl;
        }

        if (!IOobject("grad_kSGSmean", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nWriting mean SGS TKE gradient field grad_kSGSmean" << endl;
            grad_kSGSmean.write();
        }
        else
        {
            Info << "\nMean SGS TKE gradient field grad_kSGSmean already exists" << endl;
        }

        if (!IOobject("grad_p_rghAvg", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nWriting mean modified pressure gradient field grad_p_rghAvg" << endl;
            grad_p_rghAvg.write();
        }
        else
        {
            Info << "\nMean modified pressure gradient field grad_p_rghAvg already exists" << endl;
        }

        if (!IOobject("grad_TAvg", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nWriting mean resolved temperature gradient field grad_TAvg" << endl;
            grad_TAvg.write();
        }
        else
        {
            Info << "\nMean resolved temperature gradient field grad_TAvg already exists" << endl;
        }

        // //,
        // //U.boundaryField().types()//, //turbulence->nut()/Prt*fvc::grad(T),
        // //T.boundaryField().types()
        // //const volVectorField q(turbulence->nut()/Prt*fvc::grad(T));

        // // Check availability of tubulence fields

        // /*if (!IOobject("k", runTime.timeName(), mesh).headerOk())
        // {
        //     Info<< "\nWriting turbulence field k" << endl;
        //     k.write();
        // }
        // else
        // {
        //     Info<< "\nTurbulence k field already exists" << endl;
        // }*/

        // /*if (!IOobject("epsilon", runTime.timeName(), mesh).headerOk())
        // {
        //     Info<< "\nWriting turbulence field epsilon" << endl;
        //     epsilon.write();
        // }
        // else
        // {
        //     Info<< "\nTurbulence epsilon field already exists" << endl;
        // }*/

        // if (!IOobject("R", runTime.timeName(), mesh).headerOk())
        // {
        //     Info<< "\nWriting turbulence field R" << endl;
        //     R.write();
        // }
        // else
        // {
        //     Info<< "\nTurbulence R field already exists" << endl;
        // }

        // if (!IOobject("q", runTime.timeName(), mesh).headerOk())
        // {
        //     Info<< "\nWriting heat flux field q" << endl;
        //     q.write();
        // }
        // else
        // {
        //     Info<< "\nHeat flux fiel q field already exists" << endl;
        // }

        // /*if (!IOobject("omega", runTime.timeName(), mesh).headerOk())
        // {
        //     const scalar Cmu = 0.09;

        //     Info<< "creating omega" << endl;
        //     volScalarField omega
        //     (
        //         IOobject
        //         (
        //             "omega",
        //             runTime.timeName(),
        //             mesh
        //         ),
        //         epsilon/(Cmu*k),
        //         epsilon.boundaryField().types()
        //     );
        //     Info<< "\nWriting turbulence field omega" << endl;
        //     omega.write();
        // }
        // else
        // {
        //     Info<< "\nTurbulence omega field already exists" << endl;
        // }*/
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
