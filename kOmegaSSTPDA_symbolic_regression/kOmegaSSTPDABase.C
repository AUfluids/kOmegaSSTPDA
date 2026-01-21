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
            separationRelaxation     1.0;              // optional - default: 1.0 - relaxation factor for separation factor


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

Note for Developers
    The lnInclude directory is automatically created by wmake when building
    the library. It should contain only header files (.H), not implementation
    files (.C). The lnInclude directory allows other libraries (e.g., the
    compressible version) to include these headers. To regenerate lnInclude
    automatically, simply run: wmake libso

\*---------------------------------------------------------------------------*/

#include "kOmegaSSTPDABase.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

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
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
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
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), symmTensor::zero)
    ),
    Omegaij_
    (
        IOobject
        (
            IOobject::groupName("Omegaij", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedTensor("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), tensor::zero)
    ),
    Tij2_
    (
        IOobject
        (
            IOobject::groupName("Tij2", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Tij6_
    (
        IOobject
        (
            IOobject::groupName("Tij6", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Tij7_
    (
        IOobject
        (
            IOobject::groupName("Tij7", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Tij8_
    (
        IOobject
        (
            IOobject::groupName("Tij8", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Tij9_
    (
        IOobject
        (
            IOobject::groupName("Tij9", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Tij10_
    (
        IOobject
        (
            IOobject::groupName("Tij10", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
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
            15.175753
        )
    ),
    lambda2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "lambda2",
            this->coeffDict_,
            18.535054
        )
    ),
    C0_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C0",
            this->coeffDict_,
            0.0
        )
    ),
    C1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1",
            this->coeffDict_,
            0.0
        )
    ),
    C2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2",
            this->coeffDict_,
            0.0
        )
    ),
    C3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0.0
        )
    ),
    C4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C4",
            this->coeffDict_,
            0.0
        )
    ),
    C5_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C5",
            this->coeffDict_,
            0.0
        )
    ),
    separationRelaxation_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "separationRelaxation",
            this->coeffDict_,
            1.0
        )
    ),
    alpha_S_
    (
        IOobject
        (
            IOobject::groupName("alpha_S_", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
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
            1.0
        )
    ),
    A0_2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A0_2",
            this->coeffDict_,
            0.0
        )
    ),
    A1_2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A1_2",
            this->coeffDict_,
            0.0
        )
    ),
    A2_2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A2_2",
            this->coeffDict_,
            0.0
        )
    ),
    A3_2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A3_2",
            this->coeffDict_,
            0.0
        )
    ),
    A4_2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A4_2",
            this->coeffDict_,
            0.0
        )
    ),
    A5_2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A5_2",
            this->coeffDict_,
            0.0
        )
    ),
    A0_3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A0_3",
            this->coeffDict_,
            0.0
        )
    ),
    A1_3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A1_3",
            this->coeffDict_,
            0.0
        )
    ),
    A2_3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A2_3",
            this->coeffDict_,
            0.0
        )
    ),
    A3_3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A3_3",
            this->coeffDict_,
            0.0
        )
    ),
    A4_3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A4_3",
            this->coeffDict_,
            0.0
        )
    ),
    A5_3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A5_3",
            this->coeffDict_,
            0.0
        )
    ),
    A0_4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A0_4",
            this->coeffDict_,
            0.0
        )
    ),
    A1_4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A1_4",
            this->coeffDict_,
            0.0
        )
    ),
    A2_4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A2_4",
            this->coeffDict_,
            0.0
        )
    ),
    A3_4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A3_4",
            this->coeffDict_,
            0.0
        )
    ),
    A4_4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A4_4",
            this->coeffDict_,
            0.0
        )
    ),
    A5_4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A5_4",
            this->coeffDict_,
            0.0
        )
    ),
    A0_5_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A0_5",
            this->coeffDict_,
            0.0
        )
    ),
    A1_5_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A1_5",
            this->coeffDict_,
            0.0
        )
    ),
    A2_5_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A2_5",
            this->coeffDict_,
            0.0
        )
    ),
    A3_5_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A3_5",
            this->coeffDict_,
            0.0
        )
    ),
    A4_5_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A4_5",
            this->coeffDict_,
            0.0
        )
    ),
    A5_5_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A5_5",
            this->coeffDict_,
            0.0
        )
    ),
    A0_6_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A0_6",
            this->coeffDict_,
            0.0
        )
    ),
    A1_6_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A1_6",
            this->coeffDict_,
            0.0
        )
    ),
    A2_6_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A2_6",
            this->coeffDict_,
            0.0
        )
    ),
    A3_6_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A3_6",
            this->coeffDict_,
            0.0
        )
    ),
    A4_6_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A4_6",
            this->coeffDict_,
            0.0
        )
    ),
    A5_6_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A5_6",
            this->coeffDict_,
            0.0
        )
    ),
    A0_7_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A0_7",
            this->coeffDict_,
            0.0
        )
    ),
    A1_7_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A1_7",
            this->coeffDict_,
            0.0
        )
    ),
    A2_7_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A2_7",
            this->coeffDict_,
            0.0
        )
    ),
    A3_7_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A3_7",
            this->coeffDict_,
            0.0
        )
    ),
    A4_7_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A4_7",
            this->coeffDict_,
            0.0
        )
    ),
    A5_7_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A5_7",
            this->coeffDict_,
            0.0
        )
    ),
    A0_8_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A0_8",
            this->coeffDict_,
            0.0
        )
    ),
    A1_8_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A1_8",
            this->coeffDict_,
            0.0
        )
    ),
    A2_8_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A2_8",
            this->coeffDict_,
            0.0
        )
    ),
    A3_8_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A3_8",
            this->coeffDict_,
            0.0
        )
    ),
    A4_8_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A4_8",
            this->coeffDict_,
            0.0
        )
    ),
    A5_8_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A5_8",
            this->coeffDict_,
            0.0
        )
    ),
    A0_9_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A0_9",
            this->coeffDict_,
            0.0
        )
    ),
    A1_9_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A1_9",
            this->coeffDict_,
            0.0
        )
    ),
    A2_9_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A2_9",
            this->coeffDict_,
            0.0
        )
    ),
    A3_9_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A3_9",
            this->coeffDict_,
            0.0
        )
    ),
    A4_9_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A4_9",
            this->coeffDict_,
            0.0
        )
    ),
    A5_9_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A5_9",
            this->coeffDict_,
            0.0
        )
    ),
    A0_10_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A0_10",
            this->coeffDict_,
            0.0
        )
    ),
    A1_10_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A1_10",
            this->coeffDict_,
            0.0
        )
    ),
    A2_10_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A2_10",
            this->coeffDict_,
            0.0
        )
    ),
    A3_10_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A3_10",
            this->coeffDict_,
            0.0
        )
    ),
    A4_10_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A4_10",
            this->coeffDict_,
            0.0
        )
    ),
    A5_10_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A5_10",
            this->coeffDict_,
            0.0
        )
    ),
    alpha_A_2_
    (
        IOobject
        (
            IOobject::groupName("alpha_A_2", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    alpha_A_3_
    (
        IOobject
        (
            IOobject::groupName("alpha_A_3", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    alpha_A_4_
    (
        IOobject
        (
            IOobject::groupName("alpha_A_4", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    alpha_A_5_
    (
        IOobject
        (
            IOobject::groupName("alpha_A_5", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    alpha_A_6_
    (
        IOobject
        (
            IOobject::groupName("alpha_A_6", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    alpha_A_7_
    (
        IOobject
        (
            IOobject::groupName("alpha_A_7", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    alpha_A_8_
    (
        IOobject
        (
            IOobject::groupName("alpha_A_8", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    alpha_A_9_
    (
        IOobject
        (
            IOobject::groupName("alpha_A_9", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    alpha_A_10_
    (
        IOobject
        (
            IOobject::groupName("alpha_A_10", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
    ),
    anisotropyFactor_
    (
        IOobject
        (
            IOobject::groupName("anisotropyFactor", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),

    // Fields
    bijDelta_
    (
        IOobject
        (
            IOobject::groupName("bijDelta", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    Rij_
    (
        IOobject
        (
            IOobject::groupName("Rij", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimensionSet(0, 2, -2, 0, 0, 0, 0), symmTensor::zero)
    ),
    bij_
    (
        IOobject
        (
            IOobject::groupName("bij", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            writePDAFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
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
    writePDAFields_
    (
        Switch::getOrAddToDict
        (
            "writeInvariantsAndTensors",
            this->coeffDict_,
            true
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
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0)
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
    ),
    // Z-score standardisation constants for separation correction
    I1_mean_separation
    (
        dimensionedScalar("I1_mean_separation", dimless, 0.029745472322525918)
    ),
    I1_std_separation
    (
        dimensionedScalar("I1_std_separation", dimless, 0.01781867158784395)
    ),
    I2_mean_separation
    (
        dimensionedScalar("I2_mean_separation", dimless, -0.024867181279038093)
    ),
    I2_std_separation
    (
        dimensionedScalar("I2_std_separation", dimless, 0.01800275771769106)
    ),
    I3_mean_separation
    (
        dimensionedScalar("I3_mean_separation", dimless, 7.935733464624863e-06)
    ),
    I3_std_separation
    (
        dimensionedScalar("I3_std_separation", dimless, 0.00010183240848778372)
    ),
    I4_mean_separation
    (
        dimensionedScalar("I4_mean_separation", dimless, -2.3092425370425096e-07)
    ),
    I4_std_separation
    (
        dimensionedScalar("I4_std_separation", dimless, 2.7771796282956458e-05)
    ),
    I5_mean_separation
    (
        dimensionedScalar("I5_mean_separation", dimless, -0.0004981212641134251)
    ),
    I5_std_separation
    (
        dimensionedScalar("I5_std_separation", dimless, 0.0004679466893069297)
    ),
    // Z-score standardisation constants for anisotropy correction
    I1_mean_anisotropy
    (
        dimensionedScalar("I1_mean_anisotropy", dimless, 0.03679851253346419)
    ),
    I1_std_anisotropy
    (
        dimensionedScalar("I1_std_anisotropy", dimless, 0.016109597689079515)
    ),
    I2_mean_anisotropy
    (
        dimensionedScalar("I2_mean_anisotropy", dimless, -0.03681156960681101)
    ),
    I2_std_anisotropy
    (
        dimensionedScalar("I2_std_anisotropy", dimless, 0.01612305604796805)
    ),
    I3_mean_anisotropy
    (
        dimensionedScalar("I3_mean_anisotropy", dimless, 0.0)
    ),
    I3_std_anisotropy
    (
        dimensionedScalar("I3_std_anisotropy", dimless, 1.0)
    ),
    I4_mean_anisotropy
    (
        dimensionedScalar("I4_mean_anisotropy", dimless, 7.606158379651854e-05)
    ),
    I4_std_anisotropy
    (
        dimensionedScalar("I4_std_anisotropy", dimless, 0.00014946062409414713)
    ),
    I5_mean_anisotropy
    (
        dimensionedScalar("I5_mean_anisotropy", dimless, -0.0008195433737371881)
    ),
    I5_std_anisotropy
    (
        dimensionedScalar("I5_std_anisotropy", dimless, 0.00043533017048806246)
    ),
    // Flags to track which invariants/tensors are active (initialized in constructor body)
    useI1_(false),
    useI2_(false),
    useI3_(false),
    useI4_(false),
    useI5_(false),
    useTij2_(false),
    useTij3_(false),
    useSymbolicRegression_(false),
    separationExpressionStr_(""),
    anisotropyExpressionStrs_(),
    separationExpressionParser_(nullptr),
    anisotropyExpressionParsers_(),
    debugSymbolicRegression_(false),
    useTij4_(false),
    useTij5_(false),
    useTij6_(false),
    useTij7_(false),
    useTij8_(false),
    useTij9_(false),
    useTij10_(false)
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    setDecayControl(this->coeffDict_);
    
    // Read symbolic regression configuration from constructor
    // (in case read() is not called during initialization)
    useSymbolicRegression_.readIfPresent("useSymbolicRegression", this->coeffDict_);
    debugSymbolicRegression_.readIfPresent("debugSymbolicRegression", this->coeffDict_);
    
    Pout<< "========================================" << endl;
    Pout<< "[CONSTRUCTOR] SYMBOLIC REGRESSION INIT" << endl;
    Pout<< "========================================" << endl;
    Pout<< "    Symbolic regression enabled: " << useSymbolicRegression_ << endl;
    Pout<< "    Debug symbolic regression: " << debugSymbolicRegression_ << endl;
    
    // Read expression dictionaries if symbolic regression is enabled
    if (useSymbolicRegression_)
    {
        Pout<< "    [CONSTRUCTOR] Reading symbolic regression expressions..." << endl;
        
        // Read separation expression
        if (this->coeffDict_.found("separationExpressionDict"))
        {
            fileName exprFile(this->coeffDict_.lookup("separationExpressionDict"));
            // Remove "constant/" prefix if present, as runTime_.constant() already provides it
            if (exprFile.find("constant/") == 0)
            {
                exprFile = exprFile.substr(9); // Remove "constant/" (9 characters)
            }
            Pout<< "    [CONSTRUCTOR] Attempting to read separation expression from: constant/" << exprFile << endl;
            IOdictionary exprDict
            (
                IOobject
                (
                    exprFile,
                    this->runTime_.constant(),
                    this->mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            
            if (exprDict.found("expression"))
            {
                separationExpressionStr_ = exprDict.lookupOrDefault<string>("expression", "");
                Pout<< "    [CONSTRUCTOR] Loaded separation expression:" << nl
                    << "        " << separationExpressionStr_ << endl;
                
                // Create expression parser
                separationExpressionParser_ = new expressionParser(exprDict);
                
                if (separationExpressionParser_->isValid())
                {
                    Pout<< "    [CONSTRUCTOR] Separation expression compiled successfully!" << endl;
                }
                else
                {
                    Pout<< "    [CONSTRUCTOR] ERROR: Failed to compile separation expression!" << endl;
                    WarningInFunction
                        << "Failed to compile separation expression. Falling back to hardcoded coefficients." << endl;
                }
            }
            else
            {
                WarningInFunction
                    << "Expression dictionary " << exprFile
                    << " does not contain 'expression' key" << endl;
            }
        }
        else
        {
            Pout<< "    [CONSTRUCTOR] WARNING: separationExpressionDict not found in turbulenceProperties" << endl;
        }
        
        // Read consolidated anisotropy expressions
        if (this->coeffDict_.found("anisotropyExpressionsDict"))
        {
            fileName exprFile(this->coeffDict_.lookup("anisotropyExpressionsDict"));
            // Remove "constant/" prefix if present, as runTime_.constant() already provides it
            if (exprFile.find("constant/") == 0)
            {
                exprFile = exprFile.substr(9); // Remove "constant/" (9 characters)
            }
            Pout<< "    [CONSTRUCTOR] Attempting to read anisotropy expressions from: constant/" << exprFile << endl;
            IOdictionary exprDict
            (
                IOobject
                (
                    exprFile,
                    this->runTime_.constant(),
                    this->mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            
            if (exprDict.found("tensors"))
            {
                const dictionary& tensorsDict = exprDict.subDict("tensors");
                
                // Read each tensor expression
                for (label i = 2; i <= 10; ++i)
                {
                    word tensorName("Tij" + name(i));
                    if (tensorsDict.found(tensorName))
                    {
                        const dictionary& tensorDict = tensorsDict.subDict(tensorName);
                        if (tensorDict.found("expression"))
                        {
                            string exprStr = tensorDict.lookupOrDefault<string>("expression", "");
                            word key = name(i);  // Convert label to word for HashTable key
                            anisotropyExpressionStrs_.insert(key, exprStr);
                            Pout<< "    [CONSTRUCTOR] Loaded " << tensorName
                                << " expression:" << nl
                                << "        " << exprStr << endl;
                            
                            // Create expression parser for this tensor
                            // Create parser from the full anisotropy dictionary (to get variables)
                            dictionary mergedDict(exprDict);
                            mergedDict.merge(tensorDict);
                            
                            expressionParser* parser = new expressionParser(mergedDict);
                            anisotropyExpressionParsers_.insert(key, parser);
                            
                            if (parser->isValid())
                            {
                                Pout<< "    [CONSTRUCTOR] " << tensorName 
                                    << " expression compiled successfully!" << endl;
                            }
                            else
                            {
                                Pout<< "    [CONSTRUCTOR] ERROR: Failed to compile " << tensorName << " expression!" << endl;
                                WarningInFunction
                                    << "Failed to compile " << tensorName 
                                    << " expression. Falling back to hardcoded coefficients." << endl;
                            }
                        }
                    }
                }
            }
            else
            {
                WarningInFunction
                    << "Anisotropy expressions dictionary " << exprFile
                    << " does not contain 'tensors' section" << endl;
            }
        }
        else
        {
            Pout<< "    [CONSTRUCTOR] WARNING: anisotropyExpressionsDict not found in turbulenceProperties" << endl;
        }
    }
    else
    {
        Pout<< "    [CONSTRUCTOR] Symbolic regression is disabled. Using hardcoded coefficients." << endl;
    }
    
    // Check which coefficients are present and set flags accordingly
    checkActiveCoefficients();
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
void kOmegaSSTPDABase<BasicEddyViscosityModel>::checkActiveCoefficients()
{
    // Check separation coefficients (C1-C5) for invariants
    // Note: Since getOrAddToDict adds keys with default 0.0, we check if values are non-zero
    // If all coefficients for an invariant are zero (or not present), the invariant won't be calculated
    useI1_ = (mag(C1_.value()) > SMALL);
    useI2_ = (mag(C2_.value()) > SMALL);
    useI3_ = (mag(C3_.value()) > SMALL);
    useI4_ = (mag(C4_.value()) > SMALL);
    useI5_ = (mag(C5_.value()) > SMALL);
    
    // Check anisotropy coefficients for tensors
    // Set flag to true if any coefficient for this tensor is non-zero
    useTij2_ = (mag(A0_2_.value()) > SMALL || mag(A1_2_.value()) > SMALL ||
                mag(A2_2_.value()) > SMALL || mag(A3_2_.value()) > SMALL ||
                mag(A4_2_.value()) > SMALL || mag(A5_2_.value()) > SMALL);
    useTij3_ = (mag(A0_3_.value()) > SMALL || mag(A1_3_.value()) > SMALL ||
                mag(A2_3_.value()) > SMALL || mag(A3_3_.value()) > SMALL ||
                mag(A4_3_.value()) > SMALL || mag(A5_3_.value()) > SMALL);
    useTij4_ = (mag(A0_4_.value()) > SMALL || mag(A1_4_.value()) > SMALL ||
                mag(A2_4_.value()) > SMALL || mag(A3_4_.value()) > SMALL ||
                mag(A4_4_.value()) > SMALL || mag(A5_4_.value()) > SMALL);
    useTij5_ = (mag(A0_5_.value()) > SMALL || mag(A1_5_.value()) > SMALL ||
                mag(A2_5_.value()) > SMALL || mag(A3_5_.value()) > SMALL ||
                mag(A4_5_.value()) > SMALL || mag(A5_5_.value()) > SMALL);
    useTij6_ = (mag(A0_6_.value()) > SMALL || mag(A1_6_.value()) > SMALL ||
                mag(A2_6_.value()) > SMALL || mag(A3_6_.value()) > SMALL ||
                mag(A4_6_.value()) > SMALL || mag(A5_6_.value()) > SMALL);
    useTij7_ = (mag(A0_7_.value()) > SMALL || mag(A1_7_.value()) > SMALL ||
                mag(A2_7_.value()) > SMALL || mag(A3_7_.value()) > SMALL ||
                mag(A4_7_.value()) > SMALL || mag(A5_7_.value()) > SMALL);
    useTij8_ = (mag(A0_8_.value()) > SMALL || mag(A1_8_.value()) > SMALL ||
                mag(A2_8_.value()) > SMALL || mag(A3_8_.value()) > SMALL ||
                mag(A4_8_.value()) > SMALL || mag(A5_8_.value()) > SMALL);
    useTij9_ = (mag(A0_9_.value()) > SMALL || mag(A1_9_.value()) > SMALL ||
                mag(A2_9_.value()) > SMALL || mag(A3_9_.value()) > SMALL ||
                mag(A4_9_.value()) > SMALL || mag(A5_9_.value()) > SMALL);
    useTij10_ = (mag(A0_10_.value()) > SMALL || mag(A1_10_.value()) > SMALL ||
                 mag(A2_10_.value()) > SMALL || mag(A3_10_.value()) > SMALL ||
                 mag(A4_10_.value()) > SMALL || mag(A5_10_.value()) > SMALL);
}


template<class BasicEddyViscosityModel>
void kOmegaSSTPDABase<BasicEddyViscosityModel>::limitBijDelta()
{
    // Limit bijDelta to ensure physically valid Reynolds stress tensor
    // The constraint is that Rij diagonal components must be non-negative:
    // Rij = (2/3)*k*I - 2*nut*Sij + 2*k*bijDelta
    // For diagonal components: Rxx, Ryy, Rzz >= 0
    
    // Conservative bounds for bijDelta diagonal components
    // These ensure Rij diagonal stays positive even in extreme cases
    const dimensionedScalar bijMin("bijMin", dimless, -1.0/3.0);  // Tight lower bound of -1/3
    const dimensionedScalar bijMax("bijMax", dimless, 2.0/3.0);   // Tight upper bound of 2/3
    
    // Limit diagonal components (xx, yy, zz) of bijDelta
    forAll(bijDelta_, celli)
    {
        symmTensor& bij = bijDelta_[celli];
        
        // Limit diagonal components
        bij.xx() = max(min(bij.xx(), bijMax.value()), bijMin.value());
        bij.yy() = max(min(bij.yy(), bijMax.value()), bijMin.value());
        bij.zz() = max(min(bij.zz(), bijMax.value()), bijMin.value());
        
        // Optionally limit off-diagonal components to prevent extreme values
        // This helps maintain numerical stability
        const scalar offDiagMax = 0.3;
        bij.xy() = max(min(bij.xy(), offDiagMax), -offDiagMax);
        bij.xz() = max(min(bij.xz(), offDiagMax), -offDiagMax);
        bij.yz() = max(min(bij.yz(), offDiagMax), -offDiagMax);
    }
    
    // Correct boundary conditions
    bijDelta_.correctBoundaryConditions();
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
        writePDAFields_.readIfPresent("writeInvariantsAndTensors", this->coeffDict());
        
        // Symbolic regression configuration
        useSymbolicRegression_.readIfPresent("useSymbolicRegression", this->coeffDict());
        debugSymbolicRegression_.readIfPresent("debugSymbolicRegression", this->coeffDict());
        
        Pout<< "========================================" << endl;
        Pout<< "SYMBOLIC REGRESSION CONFIGURATION" << endl;
        Pout<< "========================================" << endl;
        Pout<< "    Symbolic regression enabled: " << useSymbolicRegression_ << endl;
        Pout<< "    Debug symbolic regression: " << debugSymbolicRegression_ << endl;
        Pout<< "========================================" << endl;
        
        // Read expression dictionaries if symbolic regression is enabled
        if (useSymbolicRegression_)
        {
            // Read separation expression
            if (this->coeffDict().found("separationExpressionDict"))
            {
                fileName exprFile(this->coeffDict().lookup("separationExpressionDict"));
                // Remove "constant/" prefix if present, as runTime_.constant() already provides it
                if (exprFile.find("constant/") == 0)
                {
                    exprFile = exprFile.substr(9); // Remove "constant/" (9 characters)
                }
                Pout<< "    [SR] Attempting to read separation expression from: constant/" << exprFile << endl;
                IOdictionary exprDict
                (
                    IOobject
                    (
                        exprFile,
                        this->runTime_.constant(),
                        this->mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );
                
                if (exprDict.found("expression"))
                {
                    separationExpressionStr_ = exprDict.lookupOrDefault<string>("expression", "");
                    Pout<< "    [SR] Loaded separation expression:" << nl
                        << "        " << separationExpressionStr_ << endl;
                    
                    // Create expression parser
                    if (separationExpressionParser_)
                    {
                        delete separationExpressionParser_;
                    }
                    separationExpressionParser_ = new expressionParser(exprDict);
                    
                    if (separationExpressionParser_->isValid())
                    {
                        Pout<< "    [SR] Separation expression compiled successfully!" << endl;
                    }
                    else
                    {
                        Pout<< "    [SR] ERROR: Failed to compile separation expression!" << endl;
                        WarningInFunction
                            << "Failed to compile separation expression. Falling back to hardcoded coefficients." << endl;
                    }
                }
                else
                {
                    WarningInFunction
                        << "Expression dictionary " << exprFile
                        << " does not contain 'expression' key" << endl;
                }
            }
            
            // Read consolidated anisotropy expressions
            if (this->coeffDict().found("anisotropyExpressionsDict"))
            {
                fileName exprFile(this->coeffDict().lookup("anisotropyExpressionsDict"));
                // Remove "constant/" prefix if present, as runTime_.constant() already provides it
                if (exprFile.find("constant/") == 0)
                {
                    exprFile = exprFile.substr(9); // Remove "constant/" (9 characters)
                }
                Pout<< "    [SR] Attempting to read anisotropy expressions from: constant/" << exprFile << endl;
                IOdictionary exprDict
                (
                    IOobject
                    (
                        exprFile,
                        this->runTime_.constant(),
                        this->mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );
                
                if (exprDict.found("tensors"))
                {
                    const dictionary& tensorsDict = exprDict.subDict("tensors");
                    
                    // Read each tensor expression
                    for (label i = 2; i <= 10; ++i)
                    {
                        word tensorName("Tij" + name(i));
                        if (tensorsDict.found(tensorName))
                        {
                            const dictionary& tensorDict = tensorsDict.subDict(tensorName);
                            if (tensorDict.found("expression"))
                            {
                            string exprStr = tensorDict.lookupOrDefault<string>("expression", "");
                            word key = name(i);  // Convert label to word for HashTable key
                            anisotropyExpressionStrs_.insert(key, exprStr);
                            Pout<< "    [SR] Loaded " << tensorName
                                << " expression:" << nl
                                << "        " << exprStr << endl;
                                
                                // Create expression parser for this tensor
                                if (anisotropyExpressionParsers_.found(key))
                                {
                                    delete anisotropyExpressionParsers_[key];
                                }
                                
                                // Create parser from the full anisotropy dictionary (to get variables)
                                // We need to merge the tensor-specific expression with the variables
                                dictionary mergedDict(exprDict);
                                mergedDict.merge(tensorDict);
                                
                                expressionParser* parser = new expressionParser(mergedDict);
                                anisotropyExpressionParsers_.insert(key, parser);
                                
                            if (parser->isValid())
                            {
                                Pout<< "    [SR] " << tensorName 
                                    << " expression compiled successfully!" << endl;
                            }
                            else
                            {
                                Pout<< "    [SR] ERROR: Failed to compile " << tensorName << " expression!" << endl;
                                WarningInFunction
                                    << "Failed to compile " << tensorName 
                                    << " expression. Falling back to hardcoded coefficients." << endl;
                            }
                            }
                        }
                    }
                }
                else
                {
                    WarningInFunction
                        << "Anisotropy expressions dictionary " << exprFile
                        << " does not contain 'tensors' section" << endl;
                }
            }
            else
            {
                Info<< "    Warning: anisotropyExpressionsDict not found in turbulenceProperties" << endl;
            }
        }
        else
        {
            Info<< "    Symbolic regression is disabled. Using hardcoded coefficients." << endl;
        }

        // Separation correction coefficients
        separationCorrection_.readIfPresent("separationCorrection", this->coeffDict());
        lambda1_.readIfPresent(this->coeffDict());
        lambda2_.readIfPresent(this->coeffDict());
        C0_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        C4_.readIfPresent(this->coeffDict());
        C5_.readIfPresent(this->coeffDict());
        separationRelaxation_.readIfPresent(this->coeffDict());

        // Secondary flow coefficients
        anisotropyCorrection_.readIfPresent("anisotropyCorrection", this->coeffDict());
        anisotropyRelaxation_.readIfPresent(this->coeffDict());
        A0_2_.readIfPresent(this->coeffDict());
        A1_2_.readIfPresent(this->coeffDict());
        A2_2_.readIfPresent(this->coeffDict());
        A3_2_.readIfPresent(this->coeffDict());
        A4_2_.readIfPresent(this->coeffDict());
        A5_2_.readIfPresent(this->coeffDict());
        A0_3_.readIfPresent(this->coeffDict());
        A1_3_.readIfPresent(this->coeffDict());
        A2_3_.readIfPresent(this->coeffDict());
        A3_3_.readIfPresent(this->coeffDict());
        A4_3_.readIfPresent(this->coeffDict());
        A5_3_.readIfPresent(this->coeffDict());
        A0_4_.readIfPresent(this->coeffDict());
        A1_4_.readIfPresent(this->coeffDict());
        A2_4_.readIfPresent(this->coeffDict());
        A3_4_.readIfPresent(this->coeffDict());
        A4_4_.readIfPresent(this->coeffDict());
        A5_4_.readIfPresent(this->coeffDict());
        A0_5_.readIfPresent(this->coeffDict());
        A1_5_.readIfPresent(this->coeffDict());
        A2_5_.readIfPresent(this->coeffDict());
        A3_5_.readIfPresent(this->coeffDict());
        A4_5_.readIfPresent(this->coeffDict());
        A5_5_.readIfPresent(this->coeffDict());
        A0_6_.readIfPresent(this->coeffDict());
        A1_6_.readIfPresent(this->coeffDict());
        A2_6_.readIfPresent(this->coeffDict());
        A3_6_.readIfPresent(this->coeffDict());
        A4_6_.readIfPresent(this->coeffDict());
        A5_6_.readIfPresent(this->coeffDict());
        A0_7_.readIfPresent(this->coeffDict());
        A1_7_.readIfPresent(this->coeffDict());
        A2_7_.readIfPresent(this->coeffDict());
        A3_7_.readIfPresent(this->coeffDict());
        A4_7_.readIfPresent(this->coeffDict());
        A5_7_.readIfPresent(this->coeffDict());
        A0_8_.readIfPresent(this->coeffDict());
        A1_8_.readIfPresent(this->coeffDict());
        A2_8_.readIfPresent(this->coeffDict());
        A3_8_.readIfPresent(this->coeffDict());
        A4_8_.readIfPresent(this->coeffDict());
        A5_8_.readIfPresent(this->coeffDict());
        A0_9_.readIfPresent(this->coeffDict());
        A1_9_.readIfPresent(this->coeffDict());
        A2_9_.readIfPresent(this->coeffDict());
        A3_9_.readIfPresent(this->coeffDict());
        A4_9_.readIfPresent(this->coeffDict());
        A5_9_.readIfPresent(this->coeffDict());
        A0_10_.readIfPresent(this->coeffDict());
        A1_10_.readIfPresent(this->coeffDict());
        A2_10_.readIfPresent(this->coeffDict());
        A3_10_.readIfPresent(this->coeffDict());
        A4_10_.readIfPresent(this->coeffDict());
        A5_10_.readIfPresent(this->coeffDict());

        // Decay control
        setDecayControl(this->coeffDict());
        
        // Check which coefficients are active and update flags
        checkActiveCoefficients();

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
    Sij_ = (symm(gradU));
    Omegaij_ = -0.5*(gradU - gradU.T());
    volScalarField S(sqrt(2*magSqr(symm(fvc::grad(U)))));

    // Calculate turbulent time scales tau
    volScalarField tauScale(1./max(S/a1_ + this->omegaMin_, omega_ + this->omegaMin_));
    volScalarField tauScale2(tauScale*tauScale);
    volScalarField tauScale3(tauScale2*tauScale);
    volScalarField tauScale4(tauScale2*tauScale2);
    volScalarField tauScale5(tauScale2*tauScale3);

    // Dimensionless strain/rotation tensors
    tmp<volSymmTensorField> tSdim(tauScale*Sij_);
    const volSymmTensorField& Sdim = tSdim();
    tmp<volTensorField> tOdim(tauScale*Omegaij_);
    const volTensorField& Odim = tOdim();

    // Calculate invariants (normalised by omega via tauScale powers) - only if coefficients are present
    if (useI1_ || separationCorrection_)
    {
        I1_ = (tauScale2 * tr(Sij_ & Sij_));
    }
    if (useI2_ || separationCorrection_)
    {
        I2_ = (tauScale2 * tr(Omegaij_ & Omegaij_));
    }
    if (useI3_ || separationCorrection_)
    {
        I3_ = (tauScale3 * tr(Sij_ & Sij_ & Sij_));
    }
    if (useI4_ || separationCorrection_)
    {
        I4_ = (tauScale3 * tr(Omegaij_ & Omegaij_ & Sij_));
    }
    if (useI5_ || separationCorrection_)
    {
        I5_ = (tauScale4 * tr(Omegaij_ & Omegaij_ & Sij_ & Sij_));
    }

    // Calculate base tensors (dimensionless using Sdim/Odim) - only if coefficients are present
    if (useTij2_)
    {
        Tij2_ = symm((Sdim & Odim) - (Odim & Sdim));
    }
    if (useTij3_)
    {
        Tij3_ = symm(Sdim & Sdim - (scalar(1.0)/3.0)*tr(Sdim & Sdim)*tensor::I);
    }
    if (useTij4_)
    {
        Tij4_ = symm(Odim & Odim - (scalar(1.0)/3.0)*tr(Odim & Odim)*tensor::I);
    }
    if (useTij5_)
    {
        Tij5_ = symm((Odim & Sdim & Sdim) - (Sdim & Sdim & Odim));
    }
    if (useTij6_)
    {
        Tij6_ = symm((Odim & Odim & Sdim) + (Sdim & Odim & Odim) - (scalar(2.0)/3.0)*tr(Sdim & Odim & Odim)*tensor::I);
    }
    if (useTij7_)
    {
        Tij7_ = symm((Odim & Sdim & Odim & Odim) - (Odim & Odim & Sdim & Odim));
    }
    if (useTij8_)
    {
        Tij8_ = symm((Odim & Sdim & Sdim & Sdim) - (Sdim & Sdim & Odim & Sdim));
    }
    if (useTij9_)
    {
        Tij9_ = symm((Odim & Odim & Sdim & Sdim) + (Sdim & Sdim & Odim & Odim) - (scalar(2.0)/3.0)*tr(Sdim & Sdim & Odim & Odim)*tensor::I);
    }
    if (useTij10_)
    {
        Tij10_ = symm((Odim & Sdim & Sdim & Sdim & Odim) - (Odim & Odim & Sdim & Sdim & Odim));
    }

    // Calculate invariants and base tensors cell by cell (only for active fields)
    forAll(Sij_, CellI)
    {
        if (useI1_ || separationCorrection_)
        {
            I1_[CellI] = tauScale2[CellI] * tr(Sij_[CellI] & Sij_[CellI]);
        }
        if (useI2_ || separationCorrection_)
        {
            I2_[CellI] = tauScale2[CellI] * tr(Omegaij_[CellI] & Omegaij_[CellI]);
        }
        if (useI3_ || separationCorrection_)
        {
            I3_[CellI] = tauScale3[CellI] * tr(Sij_[CellI] & Sij_[CellI] & Sij_[CellI]);
        }
        if (useI4_ || separationCorrection_)
        {
            I4_[CellI] = tauScale3[CellI] * tr(Omegaij_[CellI] & Omegaij_[CellI] & Sij_[CellI]);
        }
        if (useI5_ || separationCorrection_)
        {
            I5_[CellI] = tauScale4[CellI] * tr(Omegaij_[CellI] & Omegaij_[CellI] & Sij_[CellI] & Sij_[CellI]);
        }

        // Base tensors (dimensionless via tauScale powers)
        if (useTij2_)
        {
            Tij2_[CellI] = symm((Sdim[CellI] & Odim[CellI]) - (Odim[CellI] & Sdim[CellI]));
        }
        if (useTij3_)
        {
            Tij3_[CellI] = symm(Sdim[CellI] & Sdim[CellI] - (scalar(1.0)/3.0)*tr(Sdim[CellI] & Sdim[CellI])*tensor::I);
        }
        if (useTij4_)
        {
            Tij4_[CellI] = symm(Odim[CellI] & Odim[CellI] - (scalar(1.0)/3.0)*tr(Odim[CellI] & Odim[CellI])*tensor::I);
        }
        if (useTij5_)
        {
            Tij5_[CellI] = symm((Odim[CellI] & Sdim[CellI] & Sdim[CellI]) - (Sdim[CellI] & Sdim[CellI] & Odim[CellI]));
        }
        if (useTij6_)
        {
            Tij6_[CellI] = symm((Odim[CellI] & Odim[CellI] & Sdim[CellI]) + (Sdim[CellI] & Odim[CellI] & Odim[CellI]) - (scalar(2.0)/3.0)*tr(Sdim[CellI] & Odim[CellI] & Odim[CellI])*tensor::I);
        }
        if (useTij7_)
        {
            Tij7_[CellI] = symm((Odim[CellI] & Sdim[CellI] & Odim[CellI] & Odim[CellI]) - (Odim[CellI] & Odim[CellI] & Sdim[CellI] & Odim[CellI]));
        }
        if (useTij8_)
        {
            Tij8_[CellI] = symm((Odim[CellI] & Sdim[CellI] & Sdim[CellI] & Sdim[CellI]) - (Sdim[CellI] & Sdim[CellI] & Odim[CellI] & Sdim[CellI]));
        }
        if (useTij9_)
        {
            Tij9_[CellI] = symm((Odim[CellI] & Odim[CellI] & Sdim[CellI] & Sdim[CellI]) + (Sdim[CellI] & Sdim[CellI] & Odim[CellI] & Odim[CellI]) - (scalar(2.0)/3.0)*tr(Sdim[CellI] & Sdim[CellI] & Odim[CellI] & Odim[CellI])*tensor::I);
        }
        if (useTij10_)
        {
            Tij10_[CellI] = symm((Odim[CellI] & Sdim[CellI] & Sdim[CellI] & Sdim[CellI] & Odim[CellI]) - (Odim[CellI] & Odim[CellI] & Sdim[CellI] & Sdim[CellI] & Odim[CellI]));
        }
    }

    dimensionedScalar nutMin("nutMin", dimensionSet(0, 2, -1, 0, 0, 0, 0), 1e-9);

    // Calculate separation factor
    // alpha_S = C0 + C1*I1 + C2*I2 + C3*I3 + C4*I4 + C5*I5
    // Only include terms for invariants that have non-zero coefficients
    if (separationCorrection_)
    {
        if (useSymbolicRegression_ && separationExpressionParser_ && separationExpressionParser_->isValid())
        {
            // Evaluate symbolic regression expression
            // Create normalized invariant fields for expression evaluation
            volScalarField I1_norm((I1_ - I1_mean_separation)/I1_std_separation);
            volScalarField I2_norm((I2_ - I2_mean_separation)/I2_std_separation);
            volScalarField I3_norm((I3_ - I3_mean_separation)/I3_std_separation);
            volScalarField I4_norm((I4_ - I4_mean_separation)/I4_std_separation);
            volScalarField I5_norm((I5_ - I5_mean_separation)/I5_std_separation);
            
            // Register field variables with parser
            separationExpressionParser_->registerFieldVariable("I1", I1_norm);
            separationExpressionParser_->registerFieldVariable("I2", I2_norm);
            separationExpressionParser_->registerFieldVariable("I3", I3_norm);
            separationExpressionParser_->registerFieldVariable("I4", I4_norm);
            separationExpressionParser_->registerFieldVariable("I5", I5_norm);
            
            // Evaluate expression
            separationExpressionParser_->evaluateField(alpha_S_);
        }
        else
        {
            // Fall back to hardcoded calculation
            alpha_S_ = C0_;
            if (useI1_)
            {
                alpha_S_ += C1_*(I1_ - I1_mean_separation.value())/I1_std_separation.value();
            }
            if (useI2_)
            {
                alpha_S_ += C2_*(I2_ - I2_mean_separation.value())/I2_std_separation.value();
            }
            if (useI3_)
            {
                alpha_S_ += C3_*(I3_ - I3_mean_separation.value())/I3_std_separation.value();
            }
            if (useI4_)
            {
                alpha_S_ += C4_*(I4_ - I4_mean_separation.value())/I4_std_separation.value();
            }
            if (useI5_)
            {
                alpha_S_ += C5_*(I5_ - I5_mean_separation.value())/I5_std_separation.value();
            }
        }
        
        // Debug output
        static label debugCounter = 0;
        debugCounter++;
        if (debugSymbolicRegression_ && (debugCounter % 100 == 0))
        {
            Info<< nl << "    ========================================" << nl
                << "    [Debug] Symbolic Regression - Separation" << nl
                << "    ========================================" << nl
                << "    Iteration: " << debugCounter << nl
                << "    alpha_S statistics:" << nl
                << "        min: " << min(alpha_S_).value() << nl
                << "        max: " << max(alpha_S_).value() << nl
                << "        mean: " << average(alpha_S_).value() << nl
                << "        method: " 
                << (useSymbolicRegression_ && separationExpressionParser_ && separationExpressionParser_->isValid()
                    ? "symbolic regression" : "hardcoded coefficients") << nl;
            if (useSymbolicRegression_ && separationExpressionParser_ && separationExpressionParser_->isValid())
            {
                Info<< "        expression: " << separationExpressionStr_ << nl;
            }
            Info<< "    ========================================" << nl << endl;
        }
    }

    // Calculate target separation factor (before relaxation)
    volScalarField separationFactorTarget
    (
        pow
        (
            max
            (
                min
                (
                    (scalar(1)-pow(nut*omega_/k_, lambda1_.value())),
                    scalar(1)
                ),
                scalar(0)
            ), 
            lambda2_.value()
        )*(alpha_S_)
    );
    
    // Apply relaxation for robustness: new = old + (target - old) * relaxation
    separationFactor_ = separationFactor_ 
                      + (separationFactorTarget - separationFactor_)*separationRelaxation_;

    // Calculate secondary flow factors
    // alpha_A_i = A0_i + A1_i*I1 + A2_i*I2 + A3_i*I3 + A4_i*I4 + A5_i*I5
    // anisotropyFactor = sum of all active alpha_A_i*Tij_i
    if (anisotropyCorrection_)
    {
        if (useTij2_)
        {
            word key2 = name(2);  // Convert label to word for HashTable key
            if (useSymbolicRegression_ && anisotropyExpressionParsers_.found(key2) 
                && anisotropyExpressionParsers_[key2] && anisotropyExpressionParsers_[key2]->isValid())
            {
                // Evaluate symbolic regression expression
                // Create normalized invariant fields
                volScalarField I1_norm((I1_ - I1_mean_anisotropy)/I1_std_anisotropy);
                volScalarField I2_norm((I2_ - I2_mean_anisotropy)/I2_std_anisotropy);
                volScalarField I3_norm((I3_ - I3_mean_anisotropy)/I3_std_anisotropy);
                volScalarField I4_norm((I4_ - I4_mean_anisotropy)/I4_std_anisotropy);
                volScalarField I5_norm((I5_ - I5_mean_anisotropy)/I5_std_anisotropy);
                
                // Register field variables
                expressionParser* parser = anisotropyExpressionParsers_[key2];
                parser->registerFieldVariable("I1", I1_norm);
                parser->registerFieldVariable("I2", I2_norm);
                parser->registerFieldVariable("I3", I3_norm);
                parser->registerFieldVariable("I4", I4_norm);
                parser->registerFieldVariable("I5", I5_norm);
                
                // Evaluate expression
                parser->evaluateField(alpha_A_2_);
            }
            else
            {
                // Fall back to hardcoded calculation
                alpha_A_2_ = A0_2_;
                if (useI1_) alpha_A_2_ += A1_2_*(I1_ - I1_mean_anisotropy.value())/I1_std_anisotropy.value();
                if (useI2_) alpha_A_2_ += A2_2_*(I2_ - I2_mean_anisotropy.value())/I2_std_anisotropy.value();
                if (useI3_) alpha_A_2_ += A3_2_*(I3_ - I3_mean_anisotropy.value())/I3_std_anisotropy.value();
                if (useI4_) alpha_A_2_ += A4_2_*(I4_ - I4_mean_anisotropy.value())/I4_std_anisotropy.value();
                if (useI5_) alpha_A_2_ += A5_2_*(I5_ - I5_mean_anisotropy.value())/I5_std_anisotropy.value();
            }
            
            static label debugCounterTij2 = 0;
            debugCounterTij2++;
            if (debugSymbolicRegression_ && (debugCounterTij2 % 100 == 0))
            {
                Pout<< nl << "========================================" << nl
                    << "[DEBUG] Symbolic Regression - Tij2" << nl
                    << "========================================" << nl
                    << "    Iteration: " << debugCounterTij2 << nl
                    << "    alpha_A_2 statistics:" << nl
                    << "        min: " << min(alpha_A_2_).value() << nl
                    << "        max: " << max(alpha_A_2_).value() << nl
                    << "        mean: " << average(alpha_A_2_).value() << nl
                    << "        method: " 
                    << (useSymbolicRegression_ && anisotropyExpressionParsers_.found(key2) 
                        && anisotropyExpressionParsers_[key2] && anisotropyExpressionParsers_[key2]->isValid()
                        ? "symbolic regression" : "hardcoded coefficients") << nl;
                if (useSymbolicRegression_ && anisotropyExpressionParsers_.found(key2) 
                    && anisotropyExpressionParsers_[key2] && anisotropyExpressionParsers_[key2]->isValid())
                {
                    Pout<< "        expression: " << anisotropyExpressionStrs_[key2] << nl;
                }
                Pout<< "========================================" << nl << endl;
            }
        }
        if (useTij3_)
        {
            alpha_A_3_ = A0_3_;
            if (useI1_) alpha_A_3_ += A1_3_*(I1_ - I1_mean_anisotropy.value())/I1_std_anisotropy.value();
            if (useI2_) alpha_A_3_ += A2_3_*(I2_ - I2_mean_anisotropy.value())/I2_std_anisotropy.value();
            if (useI3_) alpha_A_3_ += A3_3_*(I3_ - I3_mean_anisotropy.value())/I3_std_anisotropy.value();
            if (useI4_) alpha_A_3_ += A4_3_*(I4_ - I4_mean_anisotropy.value())/I4_std_anisotropy.value();
            if (useI5_) alpha_A_3_ += A5_3_*(I5_ - I5_mean_anisotropy.value())/I5_std_anisotropy.value();
        }
        if (useTij4_)
        {
            alpha_A_4_ = A0_4_;
            if (useI1_) alpha_A_4_ += A1_4_*(I1_ - I1_mean_anisotropy.value())/I1_std_anisotropy.value();
            if (useI2_) alpha_A_4_ += A2_4_*(I2_ - I2_mean_anisotropy.value())/I2_std_anisotropy.value();
            if (useI3_) alpha_A_4_ += A3_4_*(I3_ - I3_mean_anisotropy.value())/I3_std_anisotropy.value();
            if (useI4_) alpha_A_4_ += A4_4_*(I4_ - I4_mean_anisotropy.value())/I4_std_anisotropy.value();
            if (useI5_) alpha_A_4_ += A5_4_*(I5_ - I5_mean_anisotropy.value())/I5_std_anisotropy.value();
        }
        if (useTij5_)
        {
            alpha_A_5_ = A0_5_;
            if (useI1_) alpha_A_5_ += A1_5_*(I1_ - I1_mean_anisotropy.value())/I1_std_anisotropy.value();
            if (useI2_) alpha_A_5_ += A2_5_*(I2_ - I2_mean_anisotropy.value())/I2_std_anisotropy.value();
            if (useI3_) alpha_A_5_ += A3_5_*(I3_ - I3_mean_anisotropy.value())/I3_std_anisotropy.value();
            if (useI4_) alpha_A_5_ += A4_5_*(I4_ - I4_mean_anisotropy.value())/I4_std_anisotropy.value();
            if (useI5_) alpha_A_5_ += A5_5_*(I5_ - I5_mean_anisotropy.value())/I5_std_anisotropy.value();
        }
        if (useTij6_)
        {
            alpha_A_6_ = A0_6_;
            if (useI1_) alpha_A_6_ += A1_6_*(I1_ - I1_mean_anisotropy.value())/I1_std_anisotropy.value();
            if (useI2_) alpha_A_6_ += A2_6_*(I2_ - I2_mean_anisotropy.value())/I2_std_anisotropy.value();
            if (useI3_) alpha_A_6_ += A3_6_*(I3_ - I3_mean_anisotropy.value())/I3_std_anisotropy.value();
            if (useI4_) alpha_A_6_ += A4_6_*(I4_ - I4_mean_anisotropy.value())/I4_std_anisotropy.value();
            if (useI5_) alpha_A_6_ += A5_6_*(I5_ - I5_mean_anisotropy.value())/I5_std_anisotropy.value();
        }
        if (useTij7_)
        {
            alpha_A_7_ = A0_7_;
            if (useI1_) alpha_A_7_ += A1_7_*(I1_ - I1_mean_anisotropy.value())/I1_std_anisotropy.value();
            if (useI2_) alpha_A_7_ += A2_7_*(I2_ - I2_mean_anisotropy.value())/I2_std_anisotropy.value();
            if (useI3_) alpha_A_7_ += A3_7_*(I3_ - I3_mean_anisotropy.value())/I3_std_anisotropy.value();
            if (useI4_) alpha_A_7_ += A4_7_*(I4_ - I4_mean_anisotropy.value())/I4_std_anisotropy.value();
            if (useI5_) alpha_A_7_ += A5_7_*(I5_ - I5_mean_anisotropy.value())/I5_std_anisotropy.value();
        }
        if (useTij8_)
        {
            alpha_A_8_ = A0_8_;
            if (useI1_) alpha_A_8_ += A1_8_*(I1_ - I1_mean_anisotropy.value())/I1_std_anisotropy.value();
            if (useI2_) alpha_A_8_ += A2_8_*(I2_ - I2_mean_anisotropy.value())/I2_std_anisotropy.value();
            if (useI3_) alpha_A_8_ += A3_8_*(I3_ - I3_mean_anisotropy.value())/I3_std_anisotropy.value();
            if (useI4_) alpha_A_8_ += A4_8_*(I4_ - I4_mean_anisotropy.value())/I4_std_anisotropy.value();
            if (useI5_) alpha_A_8_ += A5_8_*(I5_ - I5_mean_anisotropy.value())/I5_std_anisotropy.value();
        }
        if (useTij9_)
        {
            alpha_A_9_ = A0_9_;
            if (useI1_) alpha_A_9_ += A1_9_*(I1_ - I1_mean_anisotropy.value())/I1_std_anisotropy.value();
            if (useI2_) alpha_A_9_ += A2_9_*(I2_ - I2_mean_anisotropy.value())/I2_std_anisotropy.value();
            if (useI3_) alpha_A_9_ += A3_9_*(I3_ - I3_mean_anisotropy.value())/I3_std_anisotropy.value();
            if (useI4_) alpha_A_9_ += A4_9_*(I4_ - I4_mean_anisotropy.value())/I4_std_anisotropy.value();
            if (useI5_) alpha_A_9_ += A5_9_*(I5_ - I5_mean_anisotropy.value())/I5_std_anisotropy.value();
        }
        if (useTij10_)
        {
            alpha_A_10_ = A0_10_;
            if (useI1_) alpha_A_10_ += A1_10_*(I1_ - I1_mean_anisotropy.value())/I1_std_anisotropy.value();
            if (useI2_) alpha_A_10_ += A2_10_*(I2_ - I2_mean_anisotropy.value())/I2_std_anisotropy.value();
            if (useI3_) alpha_A_10_ += A3_10_*(I3_ - I3_mean_anisotropy.value())/I3_std_anisotropy.value();
            if (useI4_) alpha_A_10_ += A4_10_*(I4_ - I4_mean_anisotropy.value())/I4_std_anisotropy.value();
            if (useI5_) alpha_A_10_ += A5_10_*(I5_ - I5_mean_anisotropy.value())/I5_std_anisotropy.value();
        }
    }

    // Build anisotropyFactor from all active tensors
    anisotropyFactor_ = dimensionedSymmTensor("zero", dimless, symmTensor::zero);
    if (useTij2_) anisotropyFactor_ += alpha_A_2_*Tij2_;
    if (useTij3_) anisotropyFactor_ += alpha_A_3_*Tij3_;
    if (useTij4_) anisotropyFactor_ += alpha_A_4_*Tij4_;
    if (useTij5_) anisotropyFactor_ += alpha_A_5_*Tij5_;
    if (useTij6_) anisotropyFactor_ += alpha_A_6_*Tij6_;
    if (useTij7_) anisotropyFactor_ += alpha_A_7_*Tij7_;
    if (useTij8_) anisotropyFactor_ += alpha_A_8_*Tij8_;
    if (useTij9_) anisotropyFactor_ += alpha_A_9_*Tij9_;
    if (useTij10_) anisotropyFactor_ += alpha_A_10_*Tij10_;
    
    // Debug output for anisotropy factor
    static label anisotropyDebugCounter = 0;
    if (debugSymbolicRegression_ && anisotropyCorrection_ && (anisotropyDebugCounter++ % 100 == 0))
    {
        Info<< "    [Debug] anisotropyFactor statistics (iteration " << anisotropyDebugCounter << "):" << nl
            << "        min trace: " << min(tr(anisotropyFactor_)).value() << nl
            << "        max trace: " << max(tr(anisotropyFactor_)).value() << nl
            << "        mean trace: " << average(tr(anisotropyFactor_)).value() << nl
            << "        Active tensors: ";
        if (useTij2_) Info<< "Tij2 ";
        if (useTij3_) Info<< "Tij3 ";
        if (useTij4_) Info<< "Tij4 ";
        if (useTij5_) Info<< "Tij5 ";
        if (useTij6_) Info<< "Tij6 ";
        if (useTij7_) Info<< "Tij7 ";
        if (useTij8_) Info<< "Tij8 ";
        if (useTij9_) Info<< "Tij9 ";
        if (useTij10_) Info<< "Tij10 ";
        Info<< endl;
    }
    // Update Reynolds stress tensor
    bijDelta_ = bijDelta_ + ((nut*omega_/(k_ + this->kMin_))*anisotropyFactor_ - bijDelta_)*anisotropyRelaxation_;
    
    // Limit bijDelta to ensure physically valid Reynolds stress tensor
    limitBijDelta();
    
    // Calculate total and updated Reynolds stress tensor
    Rij_ = ((2.0/3.0)*I)*k_ - 2.0*nut*Sij_ + 2*k_*bijDelta_;
    
    // Calculate normalized Reynolds stress anisotropy tensor
    // bij = (Rij - (2/3)*k*delta_ij) / (2*k)  (dimensionless)
    const volScalarField kPlusMax(k_ + this->kMin_);
    bij_ = (Rij_ - (2.0/3.0)*k_*I) / (2.0*kPlusMax);

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
            alpha()*rho()*gamma*GbyNu0*(scalar(1) + separationFactor_())
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
