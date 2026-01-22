/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 M. J. Rincón
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
    Foam::functionObjects::calcPDAFields

Description
    Calculates PDA turbulence model invariants (I1-I5) and base tensors (Tij2-Tij10)
    from velocity field. Suitable for LES simulations.

    Fields calculated:
    - I1, I2, I3, I4, I5: Invariants
    - Tij2, Tij3, Tij4, Tij5, Tij6, Tij7, Tij8, Tij9, Tij10: Base tensors

    The time scale is calculated from the strain rate magnitude:
    \f[
        \tau = \frac{1}{S + \epsilon}
    \f]
    where S is the strain rate magnitude and \f$\epsilon\f$ is a small value
    to prevent division by zero.

SourceFiles
    calcPDAFields.C

\*---------------------------------------------------------------------------*/

#ifndef calcPDAFields_H
#define calcPDAFields_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class calcPDAFields Declaration
\*---------------------------------------------------------------------------*/

class calcPDAFields
:
    public fvMeshFunctionObject
{
    // Private Member Functions

        //- No copy construct
        calcPDAFields(const calcPDAFields&) = delete;

        //- No copy assignment
        void operator=(const calcPDAFields&) = delete;


protected:

    // Protected data

        //- Small value to prevent division by zero
        scalar epsilon_;

        //- Normalisation control
        bool normalise_;
        
        //- Use mean fields (UMean, kMean, omegaMean) or instantaneous fields (U, k, omega)
        bool useMean_;
        
        //- Field name for velocity (used for gradients)
        word UName_;


public:

    //- Runtime type information
    TypeName("calcPDAFields");


    // Constructors

        //- Construct from Time and dictionary
        calcPDAFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~calcPDAFields() = default;


    // Member Functions

        //- Read the calcPDAFields data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the calcPDAFields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

