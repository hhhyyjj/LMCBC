/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::pressureOutletLMCBCFvPatchField

Group
    grpOutletBoundaryConditions

Description
    This boundary condition provides an pressureOutletLMCBC outflow condition, based on
    solving DDt(W, field) = 0 at the boundary where \c W is the wave velocity
    and \c field is the field to which this boundary condition is applied.

    The standard (Euler, backward, CrankNicolson, localEuler) time schemes are
    supported.  Additionally an optional mechanism to relax the value at
    the boundary to a specified far-field value is provided which is
    switched on by specifying the relaxation length-scale \c lInf and the
    far-field value \c fieldInf.

    The flow/wave speed \c (w) at the outlet is provided by the virtual function
    advectionSpeed() the default implementation of which requires the name of
    the flux field \c (phi) and optionally the density \c (rho) if the
    mass-flux rather than the volumetric-flux is given.

    The flow/wave speed at the outlet can be changed by deriving a specialised
    BC from this class and over-riding advectionSpeed()  e.g. in
    waveTransmissiveFvPatchField the advectionSpeed() calculates and returns
    the flow-speed plus the acoustic wave speed creating an acoustic wave
    transmissive boundary condition.

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
        fieldInf     | value of field beyond patch | no      |
        lInf         | distance beyond patch for \c fieldInf | no |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            pressureOutletLMCBC;
        phi             phi;
    }
    \endverbatim

Note
    If \c lInf is specified, \c fieldInf will be required; \c rho is only
    required in the case of a mass-based flux.

SourceFiles
    pressureOutletLMCBCFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureOutletLMCBCFvPatchField_H
#define pressureOutletLMCBCFvPatchField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class pressureOutletLMCBCFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class pressureOutletLMCBCFvPatchField
:
    public mixedFvPatchField<Type>
{
protected:

    // Private data

        //- Name of the flux transporting the field
        word phiName_;

        word axis_;

        //- acoustic speed
        scalar cInf_;

        //- density
        scalar rho_;


public:

    //- Runtime type information
    TypeName("pressureOutletLMCBC");


    // Constructors

        //- Construct from patch and internal field
        pressureOutletLMCBCFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureOutletLMCBCFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given pressureOutletLMCBCFvPatchField
        //  onto a new patch
        pressureOutletLMCBCFvPatchField
        (
            const pressureOutletLMCBCFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pressureOutletLMCBCFvPatchField
        (
            const pressureOutletLMCBCFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new pressureOutletLMCBCFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pressureOutletLMCBCFvPatchField
        (
            const pressureOutletLMCBCFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new pressureOutletLMCBCFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

        // Evaluation functions

            //- Calculate and return the advection speed at the boundary
            virtual tmp<scalarField> advectionSpeed() const;
            virtual tmp<scalarField> advectionSpeed2() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pressureOutletLMCBCFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
