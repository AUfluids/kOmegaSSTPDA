/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
    Copyright (C) 2023-2024 M. J. Rincón
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

    This is the openFoam kOmegaSST model corrected with separation factor and secondary flow prediction,
    Progressive data augmented (PDA) 2023, A. Amarloo, M. Rincon

    More details and test cases in following publications:
    M. J. Rincon, M. Reclari, X. I. A. Yang, and M. Abkar,
    "A generalisable data-augmented turbulence model with progressive and interpretable corrections." (2025).
    arXiv preprint arXiv:2503.18568.

    A. Amarloo, M. J. Rincon, M. Reclari, and M. Abkar,
    "Progressive augmentation of RANS models for separated flow prediction
    by CFD-driven surrogate multi-objective optimisation." (2023).
    
    M. J. Rincon, A. Amarloo,  M. Reclari, X.I.A. Yang, and M. Abkar,
    "Progressive augmentation of Reynolds stress tensor models for secondary flow prediction
    by computational fluid dynamics driven surrogate optimisation" (2023).

    The model coefficients  are
    \verbatim
        kOmegaSSTPDABaseCoeffs
        {
            //original coefficients of KOSST
            alphaK1         0.85;
            alphaK2         1.0;
            alphaOmega1     0.5;
            alphaOmega2     0.856;
            beta1           0.075;
            beta2           0.0828;
            betaStar        0.09;
            gamma1          5/9;
            gamma2          0.44;
            a1              0.31;
            b1              1.0;
            c1              10.0;
            F3              no;


            // Separation coefficients
            separationCorrection      true;              // optional - default: true - off: false
            lambda1   1;              // optional - default taken from separationCorrection
            lambda2   1;              // optional - default taken from separationCorrection
            C0                  -2.070;             // optional - default taken from separationCorrection
            C1                  1.119;              // optional - default taken from separationCorrection
            C2                  -0.215;              // optional - default taken from separationCorrection


            // Anisotropy secondary flow coefficients
            anisotropyCorrection       true;              // optional - default: true - off: false
            A0                  -1.584;             // optional - default taken from anisotropyCorrection
            A1                  -0.685;              // optional - default taken from anisotropyCorrection
            A2                  -0.178;              // optional - default taken from anisotropyCorrection


            // Optional decay control
            decayControl    yes;
            kInf            \<far-field k value\>;
            omegaInf        \<far-field omega value\>;
        }
    \endverbatim



\*---------------------------------------------------------------------------*/

#include "kOmegaSSTPDABase.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
void kOmegaSSTPDABase<BasicEddyViscosityModel>::calculateInvariantsAndTensors()
{
    // Get velocity gradient (has dimensions [1/s])
    volTensorField gradU(fvc::grad(this->U_));
    
    // Calculate strain rate and rotation tensors (both have dimensions [1/s])
    Sij_ = symm(gradU);
    Omegaij_ = -0.5*(gradU - gradU.T());
    
    // Calculate magnitude of strain rate (has dimensions [1/s])
    volScalarField S(sqrt(2*magSqr(Sij_)));
    
    // Calculate time scale (has dimensions [s])
    volScalarField tScale(1.0/max(S/a1_ + this->omegaMin_, omega_ + this->omegaMin_));
    
    // Initialize fields with correct dimensions
    forAll(Sij_, cellI)
    {
        scalar tScale_i = tScale[cellI];
        scalar tScale2_i = tScale_i * tScale_i;
        scalar tScale3_i = tScale2_i * tScale_i;
        scalar tScale4_i = tScale3_i * tScale_i;
        
        // Calculate S² and Ω²
        symmTensor S2 = symm(Sij_[cellI] & Sij_[cellI]);
        tensor O2 = Omegaij_[cellI] & Omegaij_[cellI];
        
        // Calculate invariants (all dimensionless)
        I1_[cellI] = tScale2_i * tr(Sij_[cellI] & Sij_[cellI]);
        I2_[cellI] = tScale2_i * tr(Omegaij_[cellI] & Omegaij_[cellI]);
        
        // Calculate tensor bases (all dimensionless)

        // T₂ = S·Ω - Ω·S
        Tij2_[cellI] = tScale2_i * symm((Sij_[cellI] & Omegaij_[cellI]) - (Omegaij_[cellI] & Sij_[cellI]));
        
    }
}

template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTPDABase<BasicEddyViscosityModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTPDABase<BasicEddyViscosityModel>::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTPDABase<BasicEddyViscosityModel>::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTPDABase<BasicEddyViscosityModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class BasicEddyViscosityModel>
void kOmegaSSTPDABase<BasicEddyViscosityModel>::correctNut
(
    const volScalarField& S2
)
{
    // Correct the turbulence viscosity
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}


template<class BasicEddyViscosityModel>
void kOmegaSSTPDABase<BasicEddyViscosityModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}

template<class BasicEddyViscosityModel>
tmp<fvVectorMatrix> kOmegaSSTPDABase<BasicEddyViscosityModel>::divDevReff
(
    volVectorField& U
) const
{
    //Info << "In: nonlinearRST::divDevReff()" << endl;
    return
    (
      - fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)    // linear part
      + fvc::div(dev(2.*this->k_*this->bijDelta_))  // non-linear part
    );
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> kOmegaSSTPDABase<BasicEddyViscosityModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal>
kOmegaSSTPDABase<BasicEddyViscosityModel>::epsilonByk
(
    const volScalarField& F1,
    const volTensorField& gradU
) const
{
    return betaStar_*omega_();
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> kOmegaSSTPDABase<BasicEddyViscosityModel>::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
    return min
    (
        GbyNu0,
        (c1_/a1_)*betaStar_*omega_()
       *max(a1_*omega_(), b1_*F2*sqrt(S2))
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTPDABase<BasicEddyViscosityModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTPDABase<BasicEddyViscosityModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTPDABase<BasicEddyViscosityModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
kOmegaSSTPDABase<BasicEddyViscosityModel>::kOmegaSSTPDABase
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicEddyViscosityModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),


    // Added by Mario

    I1_
    (
        IOobject
        (
            IOobject::groupName("I1", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,                     // Use 'this->mesh_' to initialize the field
        dimensionedScalar("zero", dimless, 0)
    ),
    I2_
    (
        IOobject
        (
            IOobject::groupName("I2", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    I3_
    (
        IOobject
        (
            IOobject::groupName("I3", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    I4_
    (
        IOobject
        (
            IOobject::groupName("I4", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    I5_
    (
        IOobject
        (
            IOobject::groupName("I5", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    Sij_
    (
        IOobject
        (
            IOobject::groupName("Sij", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless/dimTime, symmTensor::zero)
    ),
    Omegaij_
    (
        IOobject
        (
            IOobject::groupName("Omegaij", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedTensor("zero", dimless/dimTime, tensor::zero)
    ),
    Tij1_
    (
        IOobject
        (
            IOobject::groupName("Tij1", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Tij2_
    (
        IOobject
        (
            IOobject::groupName("Tij2", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Tij3_
    (
        IOobject
        (
            IOobject::groupName("Tij3", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Tij4_
    (
        IOobject
        (
            IOobject::groupName("Tij4", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Tij5_
    (
        IOobject
        (
            IOobject::groupName("Tij5", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),

    // Model coefficients
    alphaK1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),

    // Separation correction coefficients
    separationCorrection_
    (
        Switch::getOrAddToDict
        (
            "separationCorrection",
            this->coeffDict_,
            true
        )
    ),
    lambda1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "lambda1",
            this->coeffDict_,
            18.622
        )
    ),
    lambda2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "lambda2",
            this->coeffDict_,
            4.698
        )
    ),
    C0_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C0",
            this->coeffDict_,
            -2.070
        )
    ),
    C1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.119
        )
    ),
    C2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2",
            this->coeffDict_,
            -0.215
        )
    ),

    // Secondary flow coefficients
    anisotropyCorrection_
    (
        Switch::getOrAddToDict
        (
            "anisotropyCorrection",
            this->coeffDict_,
            true
        )
    ),
    anisotropyRelaxation_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "secondary_relaxation",
            this->coeffDict_,
            0.5
        )
    ),
    A0_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A0",
            this->coeffDict_,
            -1.584
        )
    ),
    A1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A1",
            this->coeffDict_,
            -0.685
        )
    ),
    A2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A2",
            this->coeffDict_,
            -0.178
        )
    ),

    // Fields
    bijDelta_
    (
        IOobject
        (
            IOobject::groupName("bijDelta", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        0*dev(symm(fvc::grad(U)))/(omega_ + this->omegaMin_)
    ),
    Rij_
    (
        IOobject
        (
            IOobject::groupName("Rij", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        0.0*(((2.0/3.0)*I)*k_ - this->nut_*twoSymm(fvc::grad(this->U_)))
    ),
    F3_
    (
        Switch::getOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),
    y_(wallDist::New(this->mesh_).y()),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    separationFactor_
    (
        IOobject
        (
            IOobject::groupName("separationFactor", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        0*k_/k_
    ),

    // Decay control
    decayControl_
    (
        Switch::getOrAddToDict
        (
            "decayControl",
            this->coeffDict_,
            false
        )
    ),
    kInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kInf",
            this->coeffDict_,
            k_.dimensions(),
            0
        )
    ),
    omegaInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "omegaInf",
            this->coeffDict_,
            omega_.dimensions(),
            0
        )
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    setDecayControl(this->coeffDict_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
void kOmegaSSTPDABase<BasicEddyViscosityModel>::setDecayControl
(
    const dictionary& dict
)
{
    decayControl_.readIfPresent("decayControl", dict);

    if (decayControl_)
    {
        kInf_.read(dict);
        omegaInf_.read(dict);

        Info<< "    Employing decay control with kInf:" << kInf_
            << " and omegaInf:" << omegaInf_ << endl;
    }
    else
    {
        kInf_.value() = 0;
        omegaInf_.value() = 0;
    }
}


template<class BasicEddyViscosityModel>
bool kOmegaSSTPDABase<BasicEddyViscosityModel>::read()
{
    if (BasicEddyViscosityModel::read())
    {
        // Model coefficients
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        // Separation correction coefficients
        separationCorrection_.readIfPresent("separationCorrection", this->coeffDict());
        lambda1_.readIfPresent(this->coeffDict());
        lambda2_.readIfPresent(this->coeffDict());
        C0_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());

        // Secondary flow coefficients
        anisotropyCorrection_.readIfPresent("anisotropyCorrection", this->coeffDict());
        anisotropyRelaxation_.readIfPresent(this->coeffDict());
        A0_.readIfPresent(this->coeffDict());
        A1_.readIfPresent(this->coeffDict());
        A2_.readIfPresent(this->coeffDict());

        // Decay control
        setDecayControl(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicEddyViscosityModel>
void kOmegaSSTPDABase<BasicEddyViscosityModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    BasicEddyViscosityModel::correct();

    volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));

    // Calculate strain and rotation tensors
    volTensorField gradU(fvc::grad(U));
    volSymmTensorField Sij(symm(gradU));
    volTensorField Oij(-0.5*(gradU - gradU.T()));
    volScalarField S(sqrt(2*magSqr(symm(fvc::grad(U)))));

    // Calculate time scales
    volScalarField tScale(1./max(S/a1_ + this->omegaMin_, omega_ + this->omegaMin_));
    volScalarField tScale2(tScale*tScale);

    // Calculate invariants
    volScalarField I1_(tScale2 * tr(Sij & Sij));
    volScalarField I2_(tScale2 * tr(Oij & Oij));
    volSymmTensorField Tij2_(tScale2 * symm((Sij & Oij) - (Oij & Sij)));

    // Calculate invariants cell by cell
    forAll(Sij, CellI)
    {
        I1_[CellI] = tScale2[CellI] * tr(Sij[CellI] & Sij[CellI]);
        I2_[CellI] = tScale2[CellI] * tr(Oij[CellI] & Oij[CellI]);
        Tij2_[CellI] = tScale2[CellI] * symm((Sij[CellI] & Oij[CellI]) - (Oij[CellI] & Sij[CellI]));
    }

    dimensionedScalar nutMin("nutMin", dimensionSet(0, 2, -1, 0, 0, 0, 0), 1e-9);

    // Calculate separation factor
    volScalarField alpha_S(I1_ * 0.0);
    if (separationCorrection_)
    {
        alpha_S = C0_
                 + C1_*(I1_ - 2.86797085e-02) / 1.96630250e-02
                 + C2_*(I2_ - -1.21140076e-02) / 1.83587958e-02;  // z-score standardisation
    }

    separationFactor_ = pow(
        max
        (
            min
            (
                (scalar(1)-pow(nut*omega_/k_, lambda1_.value())),
                scalar(1)
            ),
            scalar(0)
        ), lambda2_.value())*(alpha_S);

    // Calculate secondary flow factor
    volScalarField alpha_A(I1_ * 0.0);
    if (anisotropyCorrection_)
    {
        alpha_A = A0_
                 + A1_*(I1_ - 4.13641572e-02) / 9.70441569e-03
                 + A2_*(I2_ - -4.13023579e-02) / 9.75952414e-03;  // z-score standardisation
    }

    // Update Reynolds stress tensor
    bijDelta_ = bijDelta_ + ((nut*omega_/(k_ + this->kMin_))*(alpha_A*Tij2_) - bijDelta_)*anisotropyRelaxation_;
    Rij_ = ((2.0/3.0)*I)*k_ - 2.0*nut*Sij + 2*k_*bijDelta_;

    // Calculate production terms
    volSymmTensorField dAij(2*k_*bijDelta_);
    volSymmTensorField P(-twoSymm(dAij & gradU));
    volScalarField Pk_bijDelta_(0.5*tr(P));

    // Calculate G/nu
    volScalarField::Internal GbyNu0
    (
        this->type() + ":GbyNu",
        ((tgradU() && dev(twoSymm(tgradU()))) + Pk_bijDelta_/(nut + nutMin))
    );

    volScalarField::Internal G(this->GName(), nut*GbyNu0);

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    // Solve omega equation
    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        GbyNu0 = GbyNu(GbyNu0, F23(), S2());
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            alpha()*rho()*gamma*(GbyNu0 + GbyNu0*separationFactor_)
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + alpha()*rho()*beta*sqr(omegaInf_)
          + Qsas(S2(), gamma, beta)
          + omegaSource()
          + fvOptions(alpha, rho, omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }

    // Solve k equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        alpha()*rho()*Pk(G)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilonByk(F1, tgradU()), k_)
      + alpha()*rho()*betaStar_*omegaInf_*kInf_
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    tgradU.clear();

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(S2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
