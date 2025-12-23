/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 M. J. Rincón
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
    normalisationFieldName_("omega"),
    kName_("k"),
    epsilonName_("epsilon"),
    useKEpsilon_(false)
{
    read(dict);
}


bool calcPDAFields::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("epsilon", epsilon_);
    dict.readIfPresent("normalise", normalise_);
    dict.readIfPresent("normalisationField", normalisationFieldName_);
    dict.readIfPresent("k", kName_);
    dict.readIfPresent("epsilon", epsilonName_);
    dict.readIfPresent("useKEpsilon", useKEpsilon_);

    return true;
}


bool calcPDAFields::execute()
{
    return true;
}


bool calcPDAFields::write()
{
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    // Calculate strain and rotation tensors
    volTensorField gradU(fvc::grad(U));
    volSymmTensorField Sij(symm(gradU));
    volTensorField Omegaij(-0.5*(gradU - gradU.T()));
    volScalarField S(sqrt(2*magSqr(symm(gradU))));

    // Calculate time scale for normalisation
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
        dimensionedScalar("one", dimless, 1.0)
    );
    
    if (normalise_)
    {
        // Try to find the normalisation field
        if (mesh_.foundObject<volScalarField>(normalisationFieldName_))
        {
            const volScalarField& normField = 
                mesh_.lookupObject<volScalarField>(normalisationFieldName_);
            tauScale = 1.0/(normField + epsilon_);
        }
        else if (normalisationFieldName_ == "omega")
        {
            // Default: use strain rate if omega not available
            tauScale = 1.0/(S + epsilon_);
        }
        else if (normalisationFieldName_ == "epsilon")
        {
            // Use epsilon if available
            if (mesh_.foundObject<volScalarField>("epsilon"))
            {
                const volScalarField& epsilon = 
                    mesh_.lookupObject<volScalarField>("epsilon");
                tauScale = 1.0/(epsilon + epsilon_);
            }
            else
            {
                tauScale = 1.0/(S + epsilon_);
            }
        }
        else if (normalisationFieldName_ == "S" || normalisationFieldName_ == "strain")
        {
            // Use strain rate magnitude
            tauScale = 1.0/(S + epsilon_);
        }
        else if (useKEpsilon_ && mesh_.foundObject<volScalarField>(kName_))
        {
            // Use k/epsilon if available
            const volScalarField& k = mesh_.lookupObject<volScalarField>(kName_);
            if (mesh_.foundObject<volScalarField>(epsilonName_))
            {
                const volScalarField& epsilon = 
                    mesh_.lookupObject<volScalarField>(epsilonName_);
                tauScale = k/(epsilon + epsilon_);
            }
            else
            {
                tauScale = 1.0/(S + epsilon_);
            }
        }
        else
        {
            // Default: use strain rate
            tauScale = 1.0/(S + epsilon_);
        }
    }
    // If normalisation is disabled, tauScale remains unity

    // Calculate time scale powers
    volScalarField tauScale2(tauScale*tauScale);
    volScalarField tauScale3(tauScale2*tauScale);
    volScalarField tauScale4(tauScale2*tauScale2);
    volScalarField tauScale5(tauScale2*tauScale3);

    // Dimensionless strain/rotation tensors (or dimensional if not normalised)
    volSymmTensorField Sdim(tauScale*Sij);
    volTensorField Odim(tauScale*Omegaij);

    // Calculate invariants
    volScalarField I1
    (
        IOobject
        (
            "I1",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        tauScale2 * tr(Sij & Sij)
    );

    volScalarField I2
    (
        IOobject
        (
            "I2",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        tauScale2 * tr(Omegaij & Omegaij)
    );

    volScalarField I3
    (
        IOobject
        (
            "I3",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        tauScale3 * tr(Sij & Sij & Sij)
    );

    volScalarField I4
    (
        IOobject
        (
            "I4",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        tauScale3 * tr(Omegaij & Omegaij & Sij)
    );

    volScalarField I5
    (
        IOobject
        (
            "I5",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        tauScale4 * tr(Omegaij & Omegaij & Sij & Sij)
    );

    // Calculate base tensors
    volSymmTensorField Tij2
    (
        IOobject
        (
            "Tij2",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        symm((Sdim & Odim) - (Odim & Sdim))
    );

    volSymmTensorField Tij3
    (
        IOobject
        (
            "Tij3",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        symm((Sdim & Sdim) - (scalar(1.0)/3.0)*tr(Sdim & Sdim)*tensor::I)
    );

    volSymmTensorField Tij4
    (
        IOobject
        (
            "Tij4",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        symm((Odim & Odim) - (scalar(1.0)/3.0)*tr(Odim & Odim)*tensor::I)
    );

    volSymmTensorField Tij5
    (
        IOobject
        (
            "Tij5",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        symm((Odim & Sdim & Sdim) - (Sdim & Sdim & Odim))
    );

    volSymmTensorField Tij6
    (
        IOobject
        (
            "Tij6",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        symm((Odim & Odim & Sdim) + (Sdim & Odim & Odim) 
             - (scalar(2.0)/3.0)*tr(Sdim & Odim & Odim)*tensor::I)
    );

    volSymmTensorField Tij7
    (
        IOobject
        (
            "Tij7",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        symm((Odim & Sdim & Odim & Odim) - (Odim & Odim & Sdim & Odim))
    );

    volSymmTensorField Tij8
    (
        IOobject
        (
            "Tij8",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        symm((Odim & Sdim & Sdim & Sdim) - (Sdim & Sdim & Odim & Sdim))
    );

    volSymmTensorField Tij9
    (
        IOobject
        (
            "Tij9",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        symm((Odim & Odim & Sdim & Sdim) + (Sdim & Sdim & Odim & Odim)
             - (scalar(2.0)/3.0)*tr(Sdim & Sdim & Odim & Odim)*tensor::I)
    );

    volSymmTensorField Tij10
    (
        IOobject
        (
            "Tij10",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        symm((Odim & Sdim & Sdim & Sdim & Odim) 
             - (Odim & Odim & Sdim & Sdim & Odim))
    );

    // Output strain-rate and rotation tensors
    volSymmTensorField SijField
    (
        IOobject
        (
            "Sij",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Sij
    );

    volTensorField OmegaijField
    (
        IOobject
        (
            "Omegaij",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Omegaij
    );

    // Write fields
    I1.write();
    I2.write();
    I3.write();
    I4.write();
    I5.write();
    Tij2.write();
    Tij3.write();
    Tij4.write();
    Tij5.write();
    Tij6.write();
    Tij7.write();
    Tij8.write();
    Tij9.write();
    Tij10.write();
    SijField.write();
    OmegaijField.write();

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// ************************************************************************* //

