/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 M. J. Rincón
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

#include "calcPDAFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(calcPDAFields, 0);

addToRunTimeSelectionTable
(
    functionObject,
    calcPDAFields,
    dictionary
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

calcPDAFields::calcPDAFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    epsilon_(1e-10),
    normalise_(true),
    useMean_(false),
    UName_("U")
{
    read(dict);
}


bool calcPDAFields::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("epsilon", epsilon_);
    dict.readIfPresent("normalise", normalise_);
    dict.readIfPresent("useMean", useMean_);
    dict.readIfPresent("U", UName_);

    return true;
}


bool calcPDAFields::execute()
{
    return true;
}


bool calcPDAFields::write()
{
    // Determine field names based on useMean flag
    word UFieldName = useMean_ ? "UMean" : (UName_.empty() ? "U" : UName_);
    word kFieldName = useMean_ ? "kMean" : "k";
    word omegaFieldName = useMean_ ? "omegaMean" : "omega";
    
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UFieldName);

    // Calculate strain and rotation tensors
    volTensorField gradU(fvc::grad(U));
    volSymmTensorField Sij(symm(gradU));
    volTensorField Omegaij(-0.5*(gradU - gradU.T()));
    volScalarField S(sqrt(2*magSqr(symm(gradU))));

    // Calculate time scale (always with time dimensions, like kOmegaSSTPDABase)
    // tauScale always has time dimensions [0 0 1 0 0 0 0]
    // The normalise flag controls whether invariants/tensors are normalized by it
    volScalarField tauScale
    (
        IOobject
        (
            "tauScale",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("one", dimensionSet(0, 0, 1, 0, 0, 0, 0), 1.0)  // Always time dimensions
    );
    
    // Try to find omega field (instantaneous or mean based on useMean flag)
    const volScalarField* omegaField = nullptr;
    
    if (mesh_.foundObject<volScalarField>(omegaFieldName))
    {
        omegaField = &mesh_.lookupObject<volScalarField>(omegaFieldName);
    }
    
    if (omegaField)
    {
        // Use the same tauScale calculation as kOmegaSSTPDABase
        // tauScale = 1.0 / max(S/a1 + omegaMin, omega + omegaMin)
        // Create directly from expression to ensure dimensions are correctly inferred
        const dimensionedScalar a1("a1", dimless, 0.31);
        const dimensionedScalar omegaMin
        (
            "omegaMin",
            omegaField->dimensions(),
            1.0e-15
        );
        // Create tauScale directly from expression (like kOmegaSSTPDABase)
        // This ensures dimensions are correctly inferred: 1.0 (dimless) / [0 0 -1 0 0 0 0] = [0 0 1 0 0 0 0]
        tauScale = 1./max(S/a1 + omegaMin, (*omegaField) + omegaMin);
    }
    else if (normalise_)
    {
        WarningInFunction
            << "Omega field '" << omegaFieldName << "' not found. "
            << "Cannot calculate tauScale for normalisation. "
            << "tauScale will remain unity (1.0 s)." << endl;
    }

    // Calculate time scale powers
    volScalarField tauScale2(tauScale*tauScale);
    volScalarField tauScale3(tauScale2*tauScale);
    volScalarField tauScale4(tauScale2*tauScale2);
    volScalarField tauScale5(tauScale2*tauScale3);

    // Dimensionless strain/rotation tensors (only if normalise is true)
    // Otherwise use dimensional Sij and Omegaij directly
    volSymmTensorField Sdim
    (
        IOobject
        (
            "Sdim",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        normalise_ ? (tauScale*Sij)() : Sij
    );
    volTensorField Odim
    (
        IOobject
        (
            "Odim",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        normalise_ ? (tauScale*Omegaij)() : Omegaij
    );

    // Calculate invariants
    // Normalise by tauScale powers only if normalise flag is true
    volScalarField I1
    (
        IOobject
        (
            "I1",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        normalise_ ? (tauScale2 * tr(Sij & Sij)) : tr(Sij & Sij)
    );
    I1.checkOut();  // Remove from object registry to prevent automatic writing

    volScalarField I2
    (
        IOobject
        (
            "I2",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        normalise_ ? (tauScale2 * tr(Omegaij & Omegaij)) : tr(Omegaij & Omegaij)
    );
    I2.checkOut();  // Remove from object registry to prevent automatic writing

    volScalarField I3
    (
        IOobject
        (
            "I3",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        normalise_ ? (tauScale3 * tr(Sij & Sij & Sij)) : tr(Sij & Sij & Sij)
    );
    I3.checkOut();

    volScalarField I4
    (
        IOobject
        (
            "I4",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        normalise_ ? (tauScale3 * tr(Omegaij & Omegaij & Sij)) : tr(Omegaij & Omegaij & Sij)
    );
    I4.checkOut();

    volScalarField I5
    (
        IOobject
        (
            "I5",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        normalise_ ? (tauScale4 * tr(Omegaij & Omegaij & Sij & Sij)) : tr(Omegaij & Omegaij & Sij & Sij)
    );
    I5.checkOut();

    // Calculate base tensors
    volSymmTensorField Tij2
    (
        IOobject
        (
            "Tij2",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        symm((Sdim & Odim) - (Odim & Sdim))
    );
    Tij2.checkOut();

    volSymmTensorField Tij3
    (
        IOobject
        (
            "Tij3",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        symm((Sdim & Sdim) - (scalar(1.0)/3.0)*tr(Sdim & Sdim)*tensor::I)
    );
    Tij3.checkOut();

    volSymmTensorField Tij4
    (
        IOobject
        (
            "Tij4",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        symm((Odim & Odim) - (scalar(1.0)/3.0)*tr(Odim & Odim)*tensor::I)
    );
    Tij4.checkOut();

    volSymmTensorField Tij5
    (
        IOobject
        (
            "Tij5",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        symm((Odim & Sdim & Sdim) - (Sdim & Sdim & Odim))
    );
    Tij5.checkOut();

    volSymmTensorField Tij6
    (
        IOobject
        (
            "Tij6",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        symm((Odim & Odim & Sdim) + (Sdim & Odim & Odim) 
             - (scalar(2.0)/3.0)*tr(Sdim & Odim & Odim)*tensor::I)
    );
    Tij6.checkOut();

    volSymmTensorField Tij7
    (
        IOobject
        (
            "Tij7",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        symm((Odim & Sdim & Odim & Odim) - (Odim & Odim & Sdim & Odim))
    );
    Tij7.checkOut();

    volSymmTensorField Tij8
    (
        IOobject
        (
            "Tij8",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        symm((Odim & Sdim & Sdim & Sdim) - (Sdim & Sdim & Odim & Sdim))
    );
    Tij8.checkOut();

    volSymmTensorField Tij9
    (
        IOobject
        (
            "Tij9",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        symm((Odim & Odim & Sdim & Sdim) + (Sdim & Sdim & Odim & Odim)
             - (scalar(2.0)/3.0)*tr(Sdim & Sdim & Odim & Odim)*tensor::I)
    );
    Tij9.checkOut();

    volSymmTensorField Tij10
    (
        IOobject
        (
            "Tij10",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        symm((Odim & Sdim & Sdim & Sdim & Odim) 
             - (Odim & Odim & Sdim & Sdim & Odim))
    );
    Tij10.checkOut();

    // Output strain-rate and rotation tensors
    volSymmTensorField SijField
    (
        IOobject
        (
            "Sij",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sij
    );
    SijField.checkOut();

    volTensorField OmegaijField
    (
        IOobject
        (
            "Omegaij",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Omegaij
    );
    OmegaijField.checkOut();

    // Retrieve turbulent kinetic energy if available (k or kMean based on useMean flag)
    const volScalarField* kPtr = nullptr;
    if (mesh_.foundObject<volScalarField>(kFieldName))
    {
        kPtr = &mesh_.lookupObject<volScalarField>(kFieldName);
    }

    // Retrieve turbulent viscosity if available (fallback to zero)
    const volScalarField* nutPtr =
        mesh_.foundObject<volScalarField>("nut")
            ? &mesh_.lookupObject<volScalarField>("nut")
            : nullptr;

    // Calculate total Reynolds stress tensor Rij
    volSymmTensorField Rij
    (
        IOobject
        (
            "Rij",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor("zero", dimensionSet(0, 2, -2, 0, 0, 0, 0), symmTensor::zero)
    );

    if (kPtr && nutPtr)
    {
        const volScalarField& k = *kPtr;
        Rij = ((2.0/3.0)*k)*symmTensor::I - 2.0*(*nutPtr)*Sij;
    }
    else if (kPtr)
    {
        const volScalarField& k = *kPtr;
        Rij = ((2.0/3.0)*k)*symmTensor::I;
    }
    Rij.checkOut();

    // Calculate normalised Reynolds stress anisotropy tensor bij
    // bij = (Rij - (2/3)*k*delta_ij) / (2*k)  (dimensionless)
    volSymmTensorField bij
    (
        IOobject
        (
            "bij",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    if (kPtr)
    {
        const volScalarField& k = *kPtr;
        const dimensionedScalar kEps("kEps", k.dimensions(), epsilon_);
        const volScalarField kPlusMax(k + kEps);
        bij = (Rij - ((2.0/3.0)*k)*symmTensor::I) / (2.0*kPlusMax);
    }
    bij.checkOut();

    // Write fields
    Info<< "    functionObjects::" << name() << " " << name() 
        << " writing fields:" << endl;
    
    Info<< "        Invariants: I1, I2, I3, I4, I5" << endl;
    I1.write();
    I2.write();
    I3.write();
    I4.write();
    I5.write();
    
    Info<< "        Base tensors: Tij2, Tij3, Tij4, Tij5, Tij6, Tij7, Tij8, Tij9, Tij10" << endl;
    Tij2.write();
    Tij3.write();
    Tij4.write();
    Tij5.write();
    Tij6.write();
    Tij7.write();
    Tij8.write();
    Tij9.write();
    Tij10.write();
    
    Info<< "        Strain/rotation tensors: Sij, Omegaij" << endl;
    SijField.write();
    OmegaijField.write();
    
    Info<< "        Reynolds stress: Rij, bij" << endl;
    Rij.write();
    bij.write();
    
    Info<< "    functionObjects::" << name() << " " << name() 
        << " wrote " << 18 << " fields" << endl;

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// ************************************************************************* //

