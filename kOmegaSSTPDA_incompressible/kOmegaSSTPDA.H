/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Class
    Foam::RASModels::kOmegaSSTPDA

Group
    grpRASTurbulence

Description
    Implementation of the k-omega-SST-PDA turbulence model for
    incompressible and compressible flows.

    The model includes:
    - Separation factor correction
    - Secondary flow prediction
    - Progressive data augmentation (PDA)

    References:
    \verbatim
        Rincon, M. J., Reclari, M., Yang, X. I. A., & Abkar, M. (2025).
        A generalisable data-augmented turbulence model with progressive 
        and interpretable corrections.
        arXiv preprint arXiv:2503.18568.

        Amarloo, A., Rincon, M. J., Reclari, M., & Abkar, M. (2023).
        Progressive augmentation of RANS models for separated flow prediction
        by CFD-driven surrogate multi-objective optimisation.

        Rincon, M. J., Amarloo, A., Reclari, M., Yang, X. I. A., & Abkar, M. (2023).
        Progressive augmentation of Reynolds stress tensor models for secondary flow prediction
        by computational fluid dynamics driven surrogate optimisation.
    \endverbatim

    Light wrapper around base class, see kOmegaSSTPDABase.

SourceFiles
    kOmegaSSTPDA.C

SeeAlso
    kOmegaSSTPDABase.H

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTPDA_H
#define kOmegaSSTPDA_H

#include "kOmegaSSTPDABase.H"
#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class kOmegaSSTPDA Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaSSTPDA
:
    public kOmegaSSTPDABase<eddyViscosity<RASModel<BasicTurbulenceModel>>>
{
    // Private Member Functions

        //- No copy construct
        kOmegaSSTPDA(const kOmegaSSTPDA&) = delete;

        //- No copy assignment
        void operator=(const kOmegaSSTPDA&) = delete;


protected:

    // Protected Member Functions

        virtual void correctNut(const volScalarField& S2);
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTPDA");


    // Constructors

        //- Construct from components
        kOmegaSSTPDA
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kOmegaSSTPDA() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmegaSSTPDA.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
