/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "pressureOutletLMCBCFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::pressureOutletLMCBCFvPatchField<Type>::pressureOutletLMCBCFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_("phi"),
    axis_("axis"),
    cInf_(340.0),
    rho_(1.225)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::pressureOutletLMCBCFvPatchField<Type>::pressureOutletLMCBCFvPatchField
(
    const pressureOutletLMCBCFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    axis_(ptf.axis_),
    cInf_(ptf.cInf_),
    rho_(ptf.rho_)
{}


template<class Type>
Foam::pressureOutletLMCBCFvPatchField<Type>::pressureOutletLMCBCFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    axis_(dict.get<word>("axis")),
    cInf_(dict.get<scalar>("cInf")),
    rho_(dict.get<scalar>("rho"))
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }

    this->refValue() = *this;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;

}


template<class Type>
Foam::pressureOutletLMCBCFvPatchField<Type>::pressureOutletLMCBCFvPatchField
(
    const pressureOutletLMCBCFvPatchField& ptpsf
)
:
    mixedFvPatchField<Type>(ptpsf),
    phiName_(ptpsf.phiName_),
    axis_(ptpsf.axis_),
    cInf_(ptpsf.cInf_),
    rho_(ptpsf.rho_)
{}


template<class Type>
Foam::pressureOutletLMCBCFvPatchField<Type>::pressureOutletLMCBCFvPatchField
(
    const pressureOutletLMCBCFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptpsf, iF),
    phiName_(ptpsf.phiName_),
    axis_(ptpsf.axis_),
    cInf_(ptpsf.cInf_),
    rho_(ptpsf.rho_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::pressureOutletLMCBCFvPatchField<Type>::advectionSpeed() const
{
    const fvPatchVectorField& uMean_f =
        this->patch().template lookupPatchField<volVectorField, vector>("uMean");
   
    Foam::vectorField n = this->patch().nf();
    Foam::scalarField uMean_n =  uMean_f & n;
    const scalarField cP(this->patch().size(), cInf_);
    return uMean_n + cP;
}

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::pressureOutletLMCBCFvPatchField<Type>::advectionSpeed2() const
{
    const fvPatchVectorField& uMean_f =
        this->patch().template lookupPatchField<volVectorField, vector>("uMean");

    Foam::vectorField n = this->patch().nf();
    Foam::scalarField uMean_n = uMean_f & n;
    const scalarField cP(this->patch().size(), cInf_);
    return cP- uMean_n;
}

template<class Type>
void Foam::pressureOutletLMCBCFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvMesh& mesh = this->internalField().mesh();
    scalar deltaT = this->db().time().deltaTValue();
    const GeometricField<Type, fvPatchField, volMesh>& field =
        this->db().objectRegistry::template
        lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            this->internalField().name()
        );
    label patchi = this->patch().index();
    const scalarField cP(this->patch().size(), cInf_);
    const scalarField rhoP(this->patch().size(), rho_);

    const scalarField lambda4(advectionSpeed());
    const scalarField lambda1(advectionSpeed2());
    const scalarField alpha(lambda4 * deltaT * this->patch().deltaCoeffs());
    this->valueFraction() = 1.0 / (1.0 + alpha);
    this->refValue() = field.oldTime().boundaryField()[patchi];
    Foam::vectorField n = this->patch().nf();

    const fvPatchVectorField& uMean_f =
        this->patch().template lookupPatchField<volVectorField, vector>("uMean");
    const fvPatchVectorField& uAcoustic_f =
        this->patch().template lookupPatchField<volVectorField, vector>("uAcoustic");

    Foam::scalarField RefGrad(this->patch().size());
    const scalarField OneField(this->patch().size(), 1);
    Foam::vectorField uAcoustic_f_snGrad = OneField * uAcoustic_f.snGrad();

    const GeometricField<vector, fvPatchField, volMesh>& uAcoustic_c =
        this->db().objectRegistry::template
        lookupObject<GeometricField<vector, fvPatchField, volMesh>>
        ("uAcoustic");

    if ("x" == axis_)
    {
        forAll(RefGrad, facei)
        {
            RefGrad[facei] =
                -0.5 * (uMean_f[facei][1] * cInf_ * rho_ / lambda4[facei]) * uAcoustic_f_snGrad[facei][1]
                - 0.5 * (uMean_f[facei][2] * cInf_ * rho_ / lambda4[facei]) * uAcoustic_f_snGrad[facei][2]
                - 0.5 * (uMean_f[facei][1] * cInf_ * rho_ / lambda1[facei]) * uAcoustic_f_snGrad[facei][1]
                - 0.5 * (uMean_f[facei][2] * cInf_ * rho_ / lambda1[facei]) * uAcoustic_f_snGrad[facei][2];
        }
    }
    else if ("y" == axis_)
    {
        forAll(RefGrad, facei)
        {
            RefGrad[facei] =
                -0.5 * (uMean_f[facei][0] * cInf_ * rho_ / lambda4[facei]) * uAcoustic_f_snGrad[facei][0]
                - 0.5 * (uMean_f[facei][2] * cInf_ * rho_ / lambda4[facei]) * uAcoustic_f_snGrad[facei][2]
                - 0.5 * (uMean_f[facei][0] * cInf_ * rho_ / lambda1[facei]) * uAcoustic_f_snGrad[facei][0]
                - 0.5 * (uMean_f[facei][2] * cInf_ * rho_ / lambda1[facei]) * uAcoustic_f_snGrad[facei][2];
        }
    }
    else
    {
        forAll(RefGrad, facei)
        {
            RefGrad[facei] =
                -0.5 * (uMean_f[facei][0] * cInf_ * rho_ / lambda4[facei]) * uAcoustic_f_snGrad[facei][0]
                - 0.5 * (uMean_f[facei][1] * cInf_ * rho_ / lambda4[facei]) * uAcoustic_f_snGrad[facei][1]
                - 0.5 * (uMean_f[facei][0] * cInf_ * rho_ / lambda1[facei]) * uAcoustic_f_snGrad[facei][0]
                - 0.5 * (uMean_f[facei][1] * cInf_ * rho_ / lambda1[facei]) * uAcoustic_f_snGrad[facei][1];
        }
    }

    this->refGrad() = Type(pTraits<Type>::one) * RefGrad;
    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::pressureOutletLMCBCFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntry("axis", axis_);
    os.writeEntry("cInf", cInf_);
    os.writeEntry("rho", rho_);

    this->writeEntry("value", os);
}


// ************************************************************************* //
