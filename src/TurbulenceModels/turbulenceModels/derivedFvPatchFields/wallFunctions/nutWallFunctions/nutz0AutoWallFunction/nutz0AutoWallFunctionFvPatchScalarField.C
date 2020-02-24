/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "nutz0AutoWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "surfaceFields.H"
//#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutz0AutoWallFunctionFvPatchScalarField::calcNut() const
{
	const surfaceScalarField& z0 =
            db().lookupObject<surfaceScalarField>(z0Name_);
    	
	const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalarField> tnutw(new scalarField(*this));
    scalarField& nutw = tnutw.ref();

    forAll(nutw, facei)
    {
        label celli = patch().faceCells()[facei];

        scalar uStar = Cmu25*sqrt(k[celli]);
        scalar yPlus = uStar*y[facei]/nuw[facei];

        scalar Edash = (y[facei] + z0[facei])/z0[facei];

	// Modified by CGS on July,2019
	scalar yPlusPrime = uStar*(y[facei] + z0[facei])/nuw[facei];

        nutw[facei] =
            nuw[facei]*(yPlus*kappa_/log(max(Edash, 1+1e-4)) - 1);

        if (debug)
        {
            Info<< "yPlus = " << yPlus
                << ", Edash = " << Edash
                << ", nutw = " << nutw[facei]
                << endl;
        }
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutz0AutoWallFunctionFvPatchScalarField::
nutz0AutoWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF),
    z0Name_("z0")
//    z0_(p.size(), 0.0)
{}


nutz0AutoWallFunctionFvPatchScalarField::
nutz0AutoWallFunctionFvPatchScalarField
(
    const nutz0AutoWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    z0Name_(ptf.z0Name_)
//    z0_(mapper(ptf.z0_))
{}


nutz0AutoWallFunctionFvPatchScalarField::
nutz0AutoWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict),
    z0Name_(dict.lookupOrDefault<word>("z0", "z0"))
//    z0_("z0", dict, p.size())
{}


nutz0AutoWallFunctionFvPatchScalarField::
nutz0AutoWallFunctionFvPatchScalarField
(
    const nutz0AutoWallFunctionFvPatchScalarField& rwfpsf
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf),
    z0Name_(rwfpsf.z0Name_)
//    z0_(rwfpsf.z0_)
{}


nutz0AutoWallFunctionFvPatchScalarField::
nutz0AutoWallFunctionFvPatchScalarField
(
    const nutz0AutoWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf, iF),
    z0Name_(rwfpsf.z0Name_)
//    z0_(rwfpsf.z0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*void nutz0AutoWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    nutkWallFunctionFvPatchScalarField::autoMap(m);
    m(z0_, z0_);
}


void nutz0AutoWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    nutkWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const nutz0AutoWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const nutz0AutoWallFunctionFvPatchScalarField>(ptf);

    z0_.rmap(nrwfpsf.z0_, addr);
}
*/

void Foam::nutz0AutoWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& z0 =
            db().lookupObject<surfaceScalarField>(z0Name_);
}

void nutz0AutoWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntryIfDifferent<word>(os, "z0", "z0", z0Name_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutz0AutoWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
