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

#include "velocityOutletLMCBCFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityOutletLMCBCFvPatchVectorField::velocityOutletLMCBCFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    cInf_(0.0),
    rho_(0.0),
    fieldInf_(Zero),
    axis_("axis")
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = Zero;
}


Foam::velocityOutletLMCBCFvPatchVectorField::velocityOutletLMCBCFvPatchVectorField
(
    const velocityOutletLMCBCFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    cInf_(ptf.cInf_),
    rho_(ptf.rho_),
    fieldInf_(ptf.fieldInf_),
    axis_(ptf.axis_)
{}


Foam::velocityOutletLMCBCFvPatchVectorField::velocityOutletLMCBCFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    cInf_(dict.get<scalar>("cInf")),
    rho_(dict.get<scalar>("rho")),
    axis_(dict.get<word>("axis"))
{
    if (dict.found("value"))
    {
        fvPatchVectorField::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchVectorField::operator=(this->patchInternalField());
    }

    this->refValue() = *this;
    this->refGrad() = Zero;
    this->valueFraction() = Zero;

    dict.lookup("fieldInf") >> fieldInf_;
}

Foam::velocityOutletLMCBCFvPatchVectorField::velocityOutletLMCBCFvPatchVectorField
(
    const velocityOutletLMCBCFvPatchVectorField& ptpsf
)
:
    directionMixedFvPatchVectorField(ptpsf),
    cInf_(ptpsf.cInf_),
    rho_(ptpsf.rho_),
    fieldInf_(ptpsf.fieldInf_),
    axis_(ptpsf.axis_)
{}

Foam::velocityOutletLMCBCFvPatchVectorField::velocityOutletLMCBCFvPatchVectorField
(
    const velocityOutletLMCBCFvPatchVectorField& ptpsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptpsf, iF),
    cInf_(ptpsf.cInf_),
    rho_(ptpsf.rho_),
    fieldInf_(ptpsf.fieldInf_),
    axis_(ptpsf.axis_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::velocityOutletLMCBCFvPatchVectorField::advectionSpeed() const
{
    const fvPatchVectorField& uMean_f =
        this->patch().template lookupPatchField<volVectorField, vector>("uMean");

    Foam::vectorField n = this->patch().nf();
    Foam::scalarField uMean_n = uMean_f & n;
    return uMean_n + tmp<scalarField>::New(this->size(), cInf_);
}

Foam::tmp<Foam::scalarField>
Foam::velocityOutletLMCBCFvPatchVectorField::advectionSpeed2() const
{
    const fvPatchVectorField& uMean_f =
        this->patch().template lookupPatchField<volVectorField, vector>("uMean");

    Foam::vectorField n = this->patch().nf();
    Foam::scalarField uMean_n = uMean_f & n;
    return tmp<scalarField>::New(this->size(), cInf_) - uMean_n;
}


void Foam::velocityOutletLMCBCFvPatchVectorField::updateCoeffs()
{
	if (this->updated())
	{
		return;
	}

	const fvMesh& mesh = this->internalField().mesh();

	word ddtScheme
		(
		 mesh.ddtScheme(this->internalField().name())
		);
	scalar deltaT = this->db().time().deltaTValue();

	const GeometricField<vector, fvPatchField, volMesh>& field =
		this->db().objectRegistry::template
		lookupObject<GeometricField<vector, fvPatchField, volMesh>>
		(
		 this->internalField().name()
		);

    const scalarField cP(this->patch().size(), cInf_);
    const scalarField rhoP(this->patch().size(), rho_);
    const  scalarField lambda4(advectionSpeed());
    const  scalarField lambda1(advectionSpeed2());
    
    const fvPatchScalarField& pAcoustic_f =
        this->patch().template lookupPatchField<volScalarField, scalar>("pAcoustic");

	label patchi = this->patch().index();
    symmTensorField Lund(symmTensorField(this->patch().size()));
    scalarField deltaCoeffs = this->patch().deltaCoeffs();

    if ("x" == axis_)
    {
        forAll(Lund, facei)
        {
            symmTensor& lws = Lund[facei];
            lws.zz() = 1.0 / (1.0 + lambda4[facei] * deltaT * deltaCoeffs[facei]);
            lws.xy() = 0.0;
            lws.xz() = 0.0;
            lws.yy() = 1.0 / (1.0 + lambda4[facei] * deltaT * deltaCoeffs[facei]);
            lws.yz() = 0.0;
            lws.xx() = 1.0 / (1.0 + lambda4[facei] * deltaT * deltaCoeffs[facei]);
        }
    }

    else if ("y" == axis_)
    {
        forAll(Lund, facei)
        {
            symmTensor& lws = Lund[facei];
            lws.zz() = 1.0 / (1.0 + lambda4[facei] * deltaT * deltaCoeffs[facei]);
            lws.xy() = 0.0;
            lws.xz() = 0.0;
            lws.yy() = 1.0 / (1.0 + lambda4[facei] * deltaT * deltaCoeffs[facei]);
            lws.yz() = 0.0;
            lws.xx() = 1.0 / (1.0 + lambda4[facei] * deltaT * deltaCoeffs[facei]);
        }
    }

    else
    {
        forAll(Lund, facei)
        {
            symmTensor& lws = Lund[facei];
            lws.zz() = 1.0 / (1.0 + lambda4[facei] * deltaT * deltaCoeffs[facei]);
            lws.xy() = 0.0;
            lws.xz() = 0.0;
            lws.yy() = 1.0 / (1.0 + lambda4[facei] * deltaT * deltaCoeffs[facei]);
            lws.yz() = 0.0;
            lws.xx() = 1.0 / (1.0 + lambda4[facei] * deltaT * deltaCoeffs[facei]);
        }
    }

    this->valueFraction() = Lund;

    this->refValue() = field.oldTime().boundaryField()[patchi];
    const fvPatchVectorField& uMean_f =
        this->patch().template lookupPatchField<volVectorField, vector>("uMean");
    const fvPatchVectorField& uAcoustic_f =
        this->patch().template lookupPatchField<volVectorField, vector>("uAcoustic");

    Foam::vectorField n = this->patch().nf();
    Foam::scalarField RefGrad(this->patch().size());
    const scalarField OneField(this->patch().size(), 1);
    Foam::vectorField uAcoustic_f_snGrad = OneField * uAcoustic_f.snGrad();

    if ("x" == axis_)
    {
        forAll(RefGrad, facei)
        {
            RefGrad[facei] =
                -0.5 * (uMean_f[facei][1] / lambda4[facei]) * uAcoustic_f_snGrad[facei][1]
                - 0.5 * (uMean_f[facei][2] / lambda4[facei]) * uAcoustic_f_snGrad[facei][2]
                + 0.5 * (uMean_f[facei][1] / lambda1[facei]) * uAcoustic_f_snGrad[facei][1]
                + 0.5 * (uMean_f[facei][2] / lambda1[facei]) * uAcoustic_f_snGrad[facei][2];
        }
    }

    else if ("y" == axis_)
    {
        forAll(RefGrad, facei)
        {
            RefGrad[facei] =
                -0.5 * (uMean_f[facei][0] / lambda4[facei]) * uAcoustic_f_snGrad[facei][0]
                - 0.5 * (uMean_f[facei][2] / lambda4[facei]) * uAcoustic_f_snGrad[facei][2]
                + 0.5 * (uMean_f[facei][0] / lambda1[facei]) * uAcoustic_f_snGrad[facei][0]
                + 0.5 * (uMean_f[facei][2] / lambda1[facei]) * uAcoustic_f_snGrad[facei][2];
        }
    }

    else
    {
        forAll(RefGrad, facei)
        {
            RefGrad[facei] =
                -0.5 * (uMean_f[facei][0] / lambda4[facei]) * uAcoustic_f_snGrad[facei][0]
                - 0.5 * (uMean_f[facei][1] / lambda4[facei]) * uAcoustic_f_snGrad[facei][1]
                + 0.5 * (uMean_f[facei][0] / lambda1[facei]) * uAcoustic_f_snGrad[facei][0]
                + 0.5 * (uMean_f[facei][1] / lambda1[facei]) * uAcoustic_f_snGrad[facei][1];
        }
    }

    this->refGrad() = RefGrad * fieldInf_;

    directionMixedFvPatchVectorField::updateCoeffs();
}

void Foam::velocityOutletLMCBCFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    os.writeEntry("cInf", cInf_);
    os.writeEntry("rho", rho_);
    os.writeEntry("fieldInf", fieldInf_);
    os.writeEntry("axis", axis_);

    this->writeEntry("value", os);
}

// ************************************************************************* //
namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        velocityOutletLMCBCFvPatchVectorField
    );
}