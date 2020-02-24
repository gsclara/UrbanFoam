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

#include "timeVaryingInletOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "inletOutletFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeVaryingInletOutletVelocityFvPatchVectorField::
timeVaryingInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    phiName_("phi"),
    flowDir_(),
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


Foam::timeVaryingInletOutletVelocityFvPatchVectorField::
timeVaryingInletOutletVelocityFvPatchVectorField
(
    const timeVaryingInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    flowDir_(ptf.flowDir_,false),
    zDir_(ptf.zDir_),
    kappa_(ptf.kappa_),
    Cmu_(ptf.Cmu_),
    Zref_(ptf.Zref_),
    Uref_(ptf.Uref_,false),
    z0_(mapper(ptf.z0_)),
    zGround_(mapper(ptf.zGround_)),
    Ustar_(mapper(ptf.Ustar_))
{}


Foam::timeVaryingInletOutletVelocityFvPatchVectorField::
timeVaryingInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    flowDir_(Function1<vector>::New("flowDir",dict)),
    zDir_(dict.lookup("zDir")),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    Zref_(readScalar(dict.lookup("Zref"))),
    Uref_(Function1<scalar>::New("Uref",dict)),
    z0_("z0", dict, p.size()),
    zGround_("zGround", dict, p.size()),
    Ustar_(p.size())
{
    /*fvPatchVectorField::operator=(vectorField("value", dict, p.size()));*/
    vector flowDir = flowDir_->value(db().time().timeOutputValue()); 
    if (mag(flowDir) < small || mag(zDir_) < small)
    {   
        FatalErrorInFunction
            << "magnitude of n or z must be greater than zero"
            << abort(FatalError);
    }

    
    // Ensure direction vectors are normalized
    flowDir /= mag(flowDir);
    zDir_ /= mag(zDir_);

    // Initialise with the value entry or internal field
    if (dict.found("value"))
    {
        fvPatchVectorField::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
	mixedFvPatchVectorField::evaluate(); 
    }
    
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::timeVaryingInletOutletVelocityFvPatchVectorField::
timeVaryingInletOutletVelocityFvPatchVectorField
(
    const timeVaryingInletOutletVelocityFvPatchVectorField& pivpvf
)
:
    mixedFvPatchVectorField(pivpvf),
    phiName_(pivpvf.phiName_),
    flowDir_(pivpvf.flowDir_,false),
    zDir_(pivpvf.zDir_),
    kappa_(pivpvf.kappa_),
    Cmu_(pivpvf.Cmu_),
    Zref_(pivpvf.Zref_),
    Uref_(pivpvf.Uref_,false),
    z0_(pivpvf.z0_),
    zGround_(pivpvf.zGround_),
    Ustar_(pivpvf.Ustar_)
{}


Foam::timeVaryingInletOutletVelocityFvPatchVectorField::
timeVaryingInletOutletVelocityFvPatchVectorField
(
    const timeVaryingInletOutletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    flowDir_(pivpvf.flowDir_,false),
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

void Foam::timeVaryingInletOutletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchVectorField::autoMap(m);
    m(z0_,z0_);
    m(zGround_,zGround_);
    m(Ustar_,Ustar_);
}


void Foam::timeVaryingInletOutletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    mixedFvPatchVectorField::rmap(ptf, addr);

    const timeVaryingInletOutletVelocityFvPatchVectorField& tiptf =
        refCast<const timeVaryingInletOutletVelocityFvPatchVectorField>
        (ptf);

    z0_.rmap(tiptf.z0_, addr);
    zGround_.rmap(tiptf.zGround_, addr);
    Ustar_.rmap(tiptf.Ustar_, addr);
}


void Foam::timeVaryingInletOutletVelocityFvPatchVectorField::updateCoeffs()
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
    vector flowDir = flowDir_->value(db().time().timeOutputValue());
    Ustar_ = kappa_*Uref/(log((Zref_ + z0_)/z0_));

    const vectorField& p = patch().Cf();
    scalarField Un
    (
        (Ustar_/kappa_)
       *log(((zDir_ & p) - zGround_ + z0_)/z0_)
    );

    refValue() = (flowDir*Un);

    //pos0 is a function that returns 1 if phip is positive or 0 
    valueFraction() = 1.0 - pos0(phip);

    mixedFvPatchVectorField::updateCoeffs();
}


void Foam::timeVaryingInletOutletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    Uref_->writeData(os);
    flowDir_->writeData(os);
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
    writeEntry(os, "value",*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::timeVaryingInletOutletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    vector flowDir = flowDir_->value(db().time().timeOutputValue());
    fvPatchField<vector>::operator=
    (
        valueFraction()*(flowDir*(flowDir & pvf))
      + (1 - valueFraction())*pvf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        timeVaryingInletOutletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
