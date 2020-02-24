/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "timeVaryingInletOutletEpsilonFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeVaryingInletOutletEpsilonFvPatchScalarField::
timeVaryingInletOutletEpsilonFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    phiName_("phi"),
    zDir_(Zero),
    kappa_(0.41),
    Cmu_(0.09),
    Zref_(0),
    Uref_(),
    z0_(0),
    zGround_(0),
    Ustar_(0)    
{
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::timeVaryingInletOutletEpsilonFvPatchScalarField::
timeVaryingInletOutletEpsilonFvPatchScalarField
(
    const timeVaryingInletOutletEpsilonFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    zDir_(ptf.zDir_),
    kappa_(ptf.kappa_),
    Cmu_(ptf.Cmu_),
    Zref_(ptf.Zref_),
    Uref_(ptf.Uref_,false),
    z0_(mapper(ptf.z0_)),
    zGround_(mapper(ptf.zGround_)),
    Ustar_(mapper(ptf.Ustar_))
{}


Foam::timeVaryingInletOutletEpsilonFvPatchScalarField::
timeVaryingInletOutletEpsilonFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    zDir_(dict.lookup("zDir")),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    Zref_(readScalar(dict.lookup("Zref"))),
    Uref_(Function1<scalar>::New("Uref",dict)),
    z0_("z0", dict, p.size()),
    zGround_("zGround", dict, p.size()),
    Ustar_(p.size())
{
    /*fvPatchScalarField::operator=(vectorField("value", dict, p.size()));*/

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        mixedFvPatchScalarField::evaluate();
    }

    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::timeVaryingInletOutletEpsilonFvPatchScalarField::
timeVaryingInletOutletEpsilonFvPatchScalarField
(
    const timeVaryingInletOutletEpsilonFvPatchScalarField& pivpvf
)
:
    mixedFvPatchScalarField(pivpvf),
    phiName_(pivpvf.phiName_),
    zDir_(pivpvf.zDir_),
    kappa_(pivpvf.kappa_),
    Cmu_(pivpvf.Cmu_),
    Zref_(pivpvf.Zref_),
    Uref_(pivpvf.Uref_,false),
    z0_(pivpvf.z0_),
    zGround_(pivpvf.zGround_),
    Ustar_(pivpvf.Ustar_)
{}


Foam::timeVaryingInletOutletEpsilonFvPatchScalarField::
timeVaryingInletOutletEpsilonFvPatchScalarField
(
    const timeVaryingInletOutletEpsilonFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    zDir_(pivpvf.zDir_),
    kappa_(pivpvf.kappa_),
    Cmu_(pivpvf.Cmu_),
    Zref_(pivpvf.Zref_),
    Uref_(pivpvf.Uref_,false),
    z0_(pivpvf.z0_),
    zGround_(pivpvf.zGround_),
    Ustar_(pivpvf.Ustar_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeVaryingInletOutletEpsilonFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    m(z0_,z0_);
    m(zGround_,zGround_);
    m(Ustar_,Ustar_);
}


void Foam::timeVaryingInletOutletEpsilonFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const timeVaryingInletOutletEpsilonFvPatchScalarField& tiptf =
        refCast<const timeVaryingInletOutletEpsilonFvPatchScalarField>
        (ptf);

    z0_.rmap(tiptf.z0_, addr);
    zGround_.rmap(tiptf.zGround_, addr);
    Ustar_.rmap(tiptf.Ustar_, addr);
}


void Foam::timeVaryingInletOutletEpsilonFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    scalar Uref = Uref_->value(db().time().timeOutputValue());
    Ustar_ = kappa_*Uref/(log((Zref_ + z0_)/z0_));

    const vectorField& p = patch().Cf();
    refValue() = (pow3(Ustar_)/(kappa_*((p & zDir_) - zGround_ +z0_)));

    //pos0 is a function that returns 1 if phip is positive or 0 
    valueFraction() = 1.0 - pos0(phip);

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::timeVaryingInletOutletEpsilonFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    Uref_->writeData(os);
    writeEntry(os, "z0", z0_) ;
    os.writeKeyword("zDir")
        << zDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa")
        << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cmu")
        << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("Zref")
        << Zref_ << token::END_STATEMENT << nl;
    writeEntry(os, "zGround", zGround_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::timeVaryingInletOutletEpsilonFvPatchScalarField::operator=
(
    const fvPatchField<scalar>& pvf
)
{
    fvPatchField<scalar>::operator=
    (
        valueFraction()*refValue()
      + (1 - valueFraction())*pvf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        timeVaryingInletOutletEpsilonFvPatchScalarField
    );
}

// ************************************************************************* //
