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
\*---------------------------------------------------------------------------*/

#include "kEpsilonABL.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// LES-RANS bij mixing ratio is default to 0 and cannot be smaller.
// It'll only be overriden if current time t > mix_startTime
scalar mix_ratio = 0.;

defineTypeNameAndDebug(kEpsilonABL, 0);
addToRunTimeSelectionTable(RASModel, kEpsilonABL, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kEpsilonABL::kEpsilonABL
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.03
        )
    ),

    Clambda_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clamda",
            coeffDict_,
            0.075
        )
    ),

    Ceps1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps1",
            coeffDict_,
            1.52
        )
    ),

    Ceps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            coeffDict_,
            1.833
        )
    ),

    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            coeffDict_,
            2.95
        )
    ),

    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            2.95
        )
    ),

    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )
    ),

    lmax_
    (
        "lmax",
        dimLength,
        1.0
    ),

    alphaB_
    (
        IOobject
        (
            "alphaB",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
       ),
        mesh_,
        dimensionedScalar("alphaB",dimless,1.00)
    ),

    Ceps1Star_
    (
        IOobject
        (
            "Ceps1Star",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
       ),
        mesh_,
        dimensionedScalar("Ceps1Star",dimless,1.52)
    ),

    Ceps3_
    (
        IOobject
        (
            "Ceps3",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
       ),
        mesh_,
        dimensionedScalar("Ceps3",dimless,1.0)
    ),

    lm_
    (
        IOobject
        (
            "lm",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
       ),
        mesh_,
        dimensionedScalar("lm",dimLength,1.0)
    ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),

    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateEpsilon("epsilon", mesh_)
    ),

    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    ),

    // // Reynolds stress from LES, read only from time 0
    // uuPrime2_
    // (
    //     IOobject
    //     (
    //         "uuPrime2",
    //         runTime_.timeName(0),  // Only read time 0
    //         mesh_,
    //         IOobject::READ_IF_PRESENT,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh_,
    //     // If not found, it should read as 0 as below
    //     dimensionedSymmTensor
    //     (
    //         "zero",
    //         dimensionedSet(0, 2, -2, 0, 0, 0, 0),
    //         symmTensor (0, 0, 0, 0, 0, 0)
    //     )
    // )

    transportDict_
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    TName_
    (
        coeffDict_.lookupOrDefault<word>("TName","T")
    ),

    T_(U.db().lookupObject<volScalarField>(TName_)),

    g_(U.db().lookupObject<uniformDimensionedVectorField>("g")),

    TRef_(transportDict_.lookup("TRef")),

    Prt_(transportDict_.lookup("Prt")),

    upVec_(vector::zero),

    // LES-RANS bij mixing related parameters
    // Read LES anisotropy tensor field bij
    bij_(U.db().lookupObject<volSymmTensorField>("bij")),
    mix_startTime_(transportDict_.lookupOrDefault<dimensionedScalar>("mix_startTime", 2000)),
    mix_duration_(transportDict_.lookupOrDefault<dimensionedScalar>("mix_duration", 2000)),
    mix_ratio_cap_(transportDict_.lookupOrDefault<dimensionedScalar>("mix_ratio_cap", 0.5)),
    mix_verbose_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "mix_verbose",
            coeffDict_,
            1
        )
    )
{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);

    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();

    upVec_ = -g_.value()/mag(g_.value());

    printCoeffs();

    if (mix_verbose_.value() > 1)
    {
        Info << "LES/ML-RANS bij mixing start time: " << mix_startTime_ << " s" << endl;
        Info << "LES/ML-RANS bij mixing duration: " << mix_duration_ << " s" << endl;
        Info << "LES/ML-RANS bij mix ratio cap: " << mix_ratio_cap_ << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void kEpsilonABL::computeMaxLengthScale()
{
    scalar numerator = 0.0;
    scalar denominator = 0.0;
    forAll(k_, i)
    {
        numerator += sqrt(k_[i]) * (mesh_.C()[i] & upVec_) * mesh_.V()[i];
        denominator += sqrt(k_[i]) * mesh_.V()[i];
    }

    reduce(numerator,sumOp<scalar>());
    reduce(denominator,sumOp<scalar>());

    dimensionedScalar n("n",dimVelocity*dimLength,numerator);
    dimensionedScalar d("d",dimVelocity,max(denominator,1.0E-10));

    lmax_ = Clambda_ * (n/d);
}

// // This is Reynolds stress from k-epsilon model
// // R = 2/3*k*I - nut*(u_i,j + u_j,i) = 2/3*k*I - 2nut*Sij
// // twoSymm returns twice the symmetric part of a tensor, in this case, the tensor is u_i,j, and twoSymm becomes (u_i,j + u_j,i)
// tmp<volSymmTensorField> kEpsilonABL::R() const
// {
//     return tmp<volSymmTensorField>
//     (
//         new volSymmTensorField
//         (
//             IOobject
//             (
//                 "R",
//                 runTime_.timeName(),
//                 mesh_,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
//             k_.boundaryField().types()
//         )
//     );
// }


// LES and RANS blended Reynolds stress with mix_ratio
// If I want to use a new function name I need to declare it in turbulenceModel.H which is a hassel
tmp<volSymmTensorField> kEpsilonABL::R() const
{
    // atio() parses the C-string str interpreting its content as an integral number, which is returned as a value of type int.
    // c_str() gets C-string equivalent of word class timeName()
    // dimensionedScalar tnow_ = atoi(runTime_.timeName().c_str());
    // Get the current time.
    scalar t = runTime_.value();
    // // LES-RANS bij mixing ratio is default to 0 and cannot be smaller.
    // // It'll only be overriden if current time t > mix_startTime
    // scalar mix_ratio = 0.;
    // dimensioned<scalar> tend_ = runTime_.endTime().value();
    // Mixing ratio of LES-RANS bij, capped on request
    if (t > mix_startTime_.value())
    {
        mix_ratio = min((t - mix_startTime_.value())*mix_ratio_cap_.value()/mix_duration_.value(),
        mix_ratio_cap_.value());
    }

    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE  // AUTO_WRITE doesn't work here as R is tmp?
            ),
            (1. - mix_ratio)*((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_))
            + mix_ratio*2.*k_*(bij_ + 1/3.*I),  // Rij_LES = 2k(1/3*I + bij)
            bij_.boundaryField().types()  // Boundary type of R is set to same as bij which are all "calculated"
            // k_.boundaryField().types()
        )
    );
}


// Effective stress tensor incl. laminar stresses
// devReff = -(nu + nut)*(u_i,j + u_j,i - 2/3*u_i,i*delta_ij)
// From Rij_RANS = 2/3*k*I - 2nut*Sij and Rij_LES = 2/3*k*I + 2k*bij,
// we gather -2nut*Sij ~ 2k*bij
// Thus we want to replace eddy-visocity assumption 2nut*Sij,
// and term nut*(u_i,j + u_j,i) -> (1 - mix_ratio)*2nut*Sij + mix_ratio(-2k*bij)
// devReff = -nu(u_i,j + u_j,i - 2/3*u_i,i*delta_ij)
//           - nut(u_i,j + u_j,i) + nut*2/3*u_i,i*delta_ij
//         = -nu(u_i,j + u_j,i - 2/3*u_i,i*delta_ij)
//           - ((1 - mix_ratio)*2nut*Sij + mix_ratio(-2k*bij)) + nut*2/3*u_i,i*delta_ij
// in which nut*2/3*u_i,i should = 0 for incompressible flow due to continuity but we keep it for stability
tmp<volSymmTensorField> kEpsilonABL::devReff() const
{
    // dimensioned<scalar> tnow_ = atoi(runTime_.timeName().c_str());
    scalar t = runTime_.value();
    // scalar mix_ratio = 0.;
    // Mixing ratio of LES-RANS bij, capped on request
    if (t > mix_startTime_.value())
    {
        mix_ratio = min((t - mix_startTime_.value())*mix_ratio_cap_.value()/mix_duration_.value(),
        mix_ratio_cap_.value());
    }

    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
        //    -nuEff()*dev(twoSymm(fvc::grad(U_)))
            - nu()*dev(twoSymm(fvc::grad(U_)))
            - ((1. - mix_ratio)*nut_*twoSymm(fvc::grad(U_)) + mix_ratio*(-2.*k_*bij_))
            + nut_*2/3.*tr(fvc::grad(U_))*I  // tr() is a scalar, thus multiply by I
        )
    );
}


// The source term for the incompressible momentum equation
// divDevReff(u_i) = -div((nu + nut)u_i,j) - div((nu + nut)(u_j,i - (u_j,j/3)*delta_ji))
//                 = -div(nu*u_i,j + nu*u_j,i - nu/3*u_j,j*delta_ji) ...[1]
//                   - div(nut*u_i,j + nut*u_j,i - nut/3*u_j,j*delta_ji) ...[2]
// We leave [1] as it is and incl. blending of LES-RANS bij to [2].
// [2] = -div(2nut*Sij) ...[3]
//       + div(nut/3*u_j,j*delta_ji) ...[4]
// We leave [4] as it is and replace 2nut*Sij in [3] with blended bij, recall 2nut*Sij ~ -2k*bij
// [3] = -div((1 - mix_ratio)*2nut*Sij + mix_ratio*(-2k*bij)) ...[5]
// Summing up, divDevReff(u_i) = [1] + [5] + [4]
//                             = -div(nu*u_i,j + nu*u_j,i - nu/3*u_j,j*delta_ji) ...[1]
//                               - div((1 - mix_ratio)*2nut*Sij)) ...[6]
//                               - div(mix_ratio*(-2k*bij)) ...[7]
//                               + div(nut/3*u_j,j*delta_ji) ...[4]
// Additionally, we're going to make [6] partially implicit, i.e. fvm again as it used to be:
// [6] = -div((1 - mix_ratio)*2nut*u_i,j) ...[8], fvm treatment
//       - div((1 - mix_ratio)*2nut*u_j,i) ...[9]
// Finally, divDevReff(u_i) = [1] + [8] + [9] + [7] + [4]
//                          = -div(nu*u_i,j + nu*u_j,i - nu*1/3*u_j,j*delta_ji)
//                            - div((1 - mix_ratio)*2nut*u_i,j)
//                            - div((1 - mix_ratio)*2nut*u_j,i)
//                            - div(mix_ratio*(-2k*bij))
//                            + div(nut/3*u_j,j*delta_ji)
// Laplacian of u_i results in 3x1 vector, i.e. same rank as u_i
// Note that ideally, due to continuity, u_j,j = 0. But for stability (u_j,j is never actually 0 in simulations), it's kept
tmp<fvVectorMatrix> kEpsilonABL::divDevReff(volVectorField& U) const
{
    // dimensioned<scalar> tnow_ = atoi(runTime_.timeName().c_str());
    scalar t = runTime_.value();
    // scalar mix_ratio = 0.;
    // Mixing ratio of LES-RANS bij, capped on request
    if (t > mix_startTime_.value())
    {
        mix_ratio = min((t - mix_startTime_.value())*mix_ratio_cap_.value()/mix_duration_.value(),
        mix_ratio_cap_.value());
    }

    return
    (
    //   - fvm::laplacian(nuEff(), U)
    //   - fvc::div(nuEff()*dev(T(fvc::grad(U))))
        - fvm::laplacian(nu(), U)  // ...[1]
        - fvc::div(nu()*dev(T(fvc::grad(U))))  // ...[1]
        - fvm::laplacian((1. - mix_ratio)*2.*nut_, U)  // ...[8]
        - fvc::div((1. - mix_ratio)*2.*nut_*T(fvc::grad(U)))  // ...[9]
        // - fvc::div((1. - mix_ratio)*nut_*twoSymm(fvc::grad(U))
        // + mix_ratio*(-2.*k_*bij_))
        - fvc::div(mix_ratio*(-2.*k_*bij_))  // ...[7]
        + fvc::div(nut_/3.*tr(T(fvc::grad(U)))*I)  // ...[4]
    );
}


// TODO: LES-RANS bij mixing for compressible flows
// The source term for the compressible momentum equation
// divDevRhoReff(rho, u_i) = rho*divDevRhoReff
// muEff = rho*(nu + nut)
tmp<fvVectorMatrix> kEpsilonABL::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


bool kEpsilonABL::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        Ceps1_.readIfPresent(coeffDict());
        Ceps2_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void kEpsilonABL::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    // dimensioned<scalar> tnow_ = atoi(runTime_.timeName().c_str());
    scalar t = runTime_.value();
    // scalar mix_ratio = 0.;
    // Mixing ratio of LES-RANS bij, capped on request
    if (t > mix_startTime_.value())
    {
        mix_ratio = min((t - mix_startTime_.value())*mix_ratio_cap_.value()/mix_duration_.value(),
        mix_ratio_cap_.value());
    }

    if (mix_verbose_.value() > 0)
    {
        Info << "Current LES/ML-RANS bij mixing ratio is " << mix_ratio << endl;

        if (mix_verbose_.value() > 1)
        {
            // cmptMin()/cmptMax() takes min/max of all componets in tensor or vector, for each cell
            // Then min()/max() takes min/max of a scalar field
            scalar bij_min = min(cmptMin(bij_));
            scalar bij_max = max(cmptMax(bij_));
            reduce(bij_min, minOp<scalar>());
            reduce(bij_max, maxOp<scalar>());
            Info << "Min LES/ML bij is " << bij_min << "; max is " << bij_max << endl;
        }
    }

    // Update length scale
    lm_ = pow(Cmu_,0.75)*pow(k_,1.5)/epsilon_;

    // Compute maximum length scale
    computeMaxLengthScale();

    // Compute the shear production term. This is where eddy-vicosity approximation comes into play in e-epsilon model
    // G = 2nut*Sij:Sij but also 2nut*Sij:grad(U),
    // where symm(grad(U)) = 0.5(u_i,j + u_j,i) = Sij,
    // and magSqr(Sij) is Sij:Sij
    // Since Rij_LES = 2/3*k*I + 2k*bij and Rij_RANS = 2/3*k*I - 2nut*Sij,
    // 2nut*Sij ~ -2k*bij
    // So with blending, G = ((1 - mix_ratio)*2nut*Sij + mix_ratio*(-2k*bij)):grad(U)
    volScalarField G("kEpsilonABL:G",
    ((1. - mix_ratio)*nut_*twoSymm(fvc::grad(U_))
    + mix_ratio*(-2.*k_*bij_))
    && fvc::grad(U_));  // Double inner dot : is double dimension reduction from 3 x 3 to 1
    if (mix_verbose_.value() > 1)
    {
        // volScalarField G("kEpsilonABL:G", 2.0*nut_*magSqr(symm(fvc::grad(U_))));
        volScalarField G_old = 2.0*nut_*magSqr(symm(fvc::grad(U_)));
        scalar g_min = min(G).value();
        scalar g_max = max(G).value();
        reduce(g_min, minOp<scalar>());
        reduce(g_max, maxOp<scalar>());
        scalar g_avg = G.weightedAverage(mesh_.V()).value();
        Info << "Min G is " << g_min << "; max is " << g_max << "; wighted mean is " << g_avg << endl;
        scalar g_old_min = min(G_old).value();
        scalar g_old_max = max(G_old).value();
        reduce(g_old_min, minOp<scalar>());
        reduce(g_old_max, maxOp<scalar>());
        scalar g_old_avg = G_old.weightedAverage(mesh_.V()).value();
        Info << "Min original G is " << g_old_min << "; max is " << g_old_max << "; wighted mean is " << g_old_avg << endl;
    }

    forAll(G,i)
    {
        if (G[i] == 0.0)
        {
            G[i] = 1.0E-10;
        }
    }
    forAll(G.boundaryField(),b)
    {
        forAll(G.boundaryField()[b],i)
        {
            if (G.boundaryField()[b][i] == 0.0)
            {
                G.boundaryField()[b][i] = 1.0E-10;
            }
        }
    }

    // Compute the buoyancy production term, should be 0 for neutral ABL below inversion layer
    // TODO: this is untouched from LES data but could inject T'T' here too
    volScalarField B("kEpsilonABL:B",(1.0/TRef_)*g_&((nut_/Prt_)*fvc::grad(T_)));

    // Compute the local gradient Richardson number.
    volScalarField Ri = -B/G;

    // Compute alphaB.
    forAll(alphaB_,i)
    {
        if (Ri[i] > 0.0)
        {
            alphaB_[i] = 1.0 - lm_[i]/lmax_.value();
        }
        else
        {
            alphaB_[i] = 1.0 - (1.0 + (Ceps2_.value() - 1.0) / (Ceps2_.value() - Ceps1_.value())) * lm_[i]/lmax_.value();
        }
    }

    // Compute Ceps1Star.
    Ceps1Star_ = Ceps1_ + (Ceps2_ - Ceps1_)*(lm_/lmax_);

    // Compute Ceps3
    Ceps3_ = (Ceps1_ - Ceps2_)*alphaB_ + 1.0;

    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        Ceps1Star_*G*epsilon_/k_
      + Ceps3_*B*epsilon_/k_
      - fvm::Sp(Ceps2_*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

    epsEqn().boundaryManipulate(epsilon_.boundaryField());

    solve(epsEqn);
    bound(epsilon_, epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      + B
      - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity
    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
