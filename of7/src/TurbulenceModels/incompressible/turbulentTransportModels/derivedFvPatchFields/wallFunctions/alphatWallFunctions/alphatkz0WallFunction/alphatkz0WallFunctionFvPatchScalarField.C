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

#include "alphatkz0WallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar alphatkz0WallFunctionFvPatchScalarField::tolerance_ = 0.01;
label alphatkz0WallFunctionFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatkz0WallFunctionFvPatchScalarField::
alphatkz0WallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    z0_terrain_(0.0),
    Prt_(0.85)
{}


alphatkz0WallFunctionFvPatchScalarField::
alphatkz0WallFunctionFvPatchScalarField
(
    const alphatkz0WallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    z0_terrain_(ptf.z0_terrain_),
    Prt_(ptf.Prt_)
{}


alphatkz0WallFunctionFvPatchScalarField::
alphatkz0WallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    z0_terrain_(readScalar(dict.lookup("z0_terrain"))),
    Prt_(readScalar(dict.lookup("Prt")))  // force read to avoid ambiguity
{}


alphatkz0WallFunctionFvPatchScalarField::
alphatkz0WallFunctionFvPatchScalarField
(
    const alphatkz0WallFunctionFvPatchScalarField& wfpsf
)
:
    fixedValueFvPatchScalarField(wfpsf),
    z0_terrain_(wfpsf.z0_terrain_),
    Prt_(wfpsf.Prt_)
{}


alphatkz0WallFunctionFvPatchScalarField::
alphatkz0WallFunctionFvPatchScalarField
(
    const alphatkz0WallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wfpsf, iF),
    z0_terrain_(wfpsf.z0_terrain_),
    Prt_(wfpsf.Prt_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatkz0WallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    // Retrieve turbulence properties from model

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);

    const scalar Cmu25 = pow(nutw.Cmu(), 0.25);
    const scalarField& y = turbModel.y()[patchi];
    const tmp<volScalarField> tnu = turbModel.nu();
    const volScalarField& nu = tnu();
    const scalarField& nuw = nu.boundaryField()[patchi];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const IOdictionary& transportProperties =
        db().lookupObject<IOdictionary>("transportProperties");

    // Molecular Prandtl number
    const scalar Pr
    (
        dimensionedScalar
        (
            "Pr",
            dimless,
            transportProperties.lookup("Pr")
        ).value()
    );

    // Populate boundary values
    scalarField& alphatw = *this;
    forAll(alphatw, facei)
    {
        const label celli = patch().faceCells()[facei];

        const scalar uStar = Cmu25*sqrt(k[celli]);

        const scalar Edash = (y[facei]+z0_terrain_)/(z0_terrain_+ 1e-4);

        // Update turbulent thermal condctivity    
        alphatw[facei] =
            uStar*nutw.kappa()*y[facei]/(Prt_*log(max(Edash, 1 + 1e-4)))
          + nuw[facei]/Pr;
    }

    // lower bound values to avoid unrealistic
    // negative temperatures on the ground
    alphatw = max(alphatw, scalar(0.01));

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void alphatkz0WallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry(os, "z0_terrain", z0_terrain_);
    writeEntry(os, "Prt", Prt_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatkz0WallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
