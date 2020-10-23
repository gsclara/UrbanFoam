/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "nutkz0WallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutkz0WallFunctionFvPatchScalarField::nut() const
{
    const label patchi = patch().index();

    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
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

        scalar Edash = (y[facei] + z0_[facei])/z0_[facei];

        // Modified by CGS on October,2020
        scalar yPlusPrime = uStar*(y[facei] + z0_[facei])/nuw[facei];

        nutw[facei] =
            nuw[facei]*(yPlus*kappa_/log(max(Edash, 1+1e-4)) - 1);
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutkz0WallFunctionFvPatchScalarField::nutkz0WallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF),
    z0_(p.size(), 0.0)
{}


nutkz0WallFunctionFvPatchScalarField::nutkz0WallFunctionFvPatchScalarField
(
    const nutkz0WallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    z0_(mapper(ptf.z0_))
{}


nutkz0WallFunctionFvPatchScalarField::nutkz0WallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict),
    z0_("z0", dict, p.size())
{}


nutkz0WallFunctionFvPatchScalarField::nutkz0WallFunctionFvPatchScalarField
(
    const nutkz0WallFunctionFvPatchScalarField& rwfpsf
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf),
    z0_(rwfpsf.z0_)
{}


nutkz0WallFunctionFvPatchScalarField::nutkz0WallFunctionFvPatchScalarField
(
    const nutkz0WallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf, iF),
    z0_(rwfpsf.z0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//tmp<scalarField> nutkz0WallFunctionFvPatchScalarField::yPlus() const
//{
//    const label patchi = patch().index();
//
//    const momentumTransportModel& turbModel =
//        db().lookupObject<momentumTransportModel>
//        (
//            IOobject::groupName
//            (
//                momentumTransportModel::typeName,
//                internalField().group()
//            )
//        );
//
//    const scalarField& y = turbModel.y()[patchi];
//
//    const tmp<volScalarField> tk = turbModel.k();
//    const volScalarField& k = tk();
//    tmp<scalarField> kwc = k.boundaryField()[patchi].patchInternalField();
//    const tmp<scalarField> tnuw = turbModel.nu(patchi);
//    const scalarField& nuw = tnuw();
//
//    return pow025(Cmu_)*y*sqrt(kwc)/nuw;
//}
// ip-try
void nutkz0WallFunctionFvPatchScalarField::autoMap
(
     const fvPatchFieldMapper& m
     )
{
    nutkWallFunctionFvPatchScalarField::autoMap(m);
    m(z0_, z0_);
}


void nutkz0WallFunctionFvPatchScalarField::rmap
(
     const fvPatchScalarField& ptf,
         const labelList& addr
         )
{
    nutkWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const nutkz0WallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const nutkz0WallFunctionFvPatchScalarField>(ptf);

    z0_.rmap(nrwfpsf.z0_, addr);
}

void nutkz0WallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry(os, "z0", z0_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutkz0WallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
