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

        Info << "\nRetrieving mean velocity gradient field grad_U" << endl;
        volTensorField grad_U
        (
            IOobject
            (
                "grad_U",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(U)
        );

        Info << "\nRetrieving mean TKE gradient field grad_k" << endl;
        volVectorField grad_k
        (
            IOobject
            (
                "grad_k",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(k)
        );

        Info << "\nRetrieving mean modified pressure gradient field grad_p_rgh" << endl;
        volVectorField grad_p_rgh
        (
            IOobject
            (
                "grad_p_rgh",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(p_rgh)
        );

        Info << "\nRetrieving mean temperature gradient field grad_T" << endl;
        volVectorField grad_T
        (
            IOobject
            (
                "grad_T",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(T)
        );

        if (!IOobject("grad_U", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nWriting mean velocity gradient field grad_U" << endl;
            grad_U.write();
        }
        else
        {
            Info << "\nMean velocity gradient field grad_U already exists" << endl;
        }

        if (!IOobject("grad_k", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nWriting mean TKE gradient field grad_k" << endl;
            grad_k.write();
        }
        else
        {
            Info << "\nMean TKE gradient field grad_k already exists" << endl;
        }

        if (!IOobject("grad_p_rgh", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nWriting mean modified pressure gradient field grad_p_rgh" << endl;
            grad_p_rgh.write();
        }
        else
        {
            Info << "\nMean modified pressure gradient field grad_p_rgh already exists" << endl;
        }

        if (!IOobject("grad_T", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nWriting mean temperature gradient field grad_T" << endl;
            grad_T.write();
        }
        else
        {
            Info << "\nMean temperature gradient field grad_T already exists" << endl;
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
