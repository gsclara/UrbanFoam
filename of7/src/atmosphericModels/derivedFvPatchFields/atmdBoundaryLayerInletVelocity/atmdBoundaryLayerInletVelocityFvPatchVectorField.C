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

#include "atmdBoundaryLayerInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmdBoundaryLayerInletVelocityFvPatchVectorField::
atmdBoundaryLayerInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchVectorField(p, iF),
    atmdBoundaryLayer(p.patch())
{}


atmdBoundaryLayerInletVelocityFvPatchVectorField::
atmdBoundaryLayerInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchVectorField(p, iF),
    atmdBoundaryLayer(patch().Cf(),  dict, p.patch())
{
    phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    refValue() = U(patch().Cf());
    refGrad() = Zero;
    valueFraction() = 1;

    if (dict.found("value"))
    {
        vectorField::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        vectorField::operator=(refValue());
    }
}


atmdBoundaryLayerInletVelocityFvPatchVectorField::
atmdBoundaryLayerInletVelocityFvPatchVectorField
(
    const atmdBoundaryLayerInletVelocityFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchVectorField(pvf, p, iF, mapper),
    atmdBoundaryLayer(pvf, p, mapper)
{}


atmdBoundaryLayerInletVelocityFvPatchVectorField::
atmdBoundaryLayerInletVelocityFvPatchVectorField
(
    const atmdBoundaryLayerInletVelocityFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchVectorField(pvf, iF),
    atmdBoundaryLayer(pvf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmdBoundaryLayerInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    inletOutletFvPatchVectorField::autoMap(m);
    atmdBoundaryLayer::autoMap(m);
}


void atmdBoundaryLayerInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& pvf,
    const labelList& addr
)
{
    inletOutletFvPatchVectorField::rmap(pvf, addr);

    const atmdBoundaryLayerInletVelocityFvPatchVectorField& blpvf =
        refCast<const atmdBoundaryLayerInletVelocityFvPatchVectorField>(pvf);

    atmdBoundaryLayer::rmap(blpvf, addr);
}


void atmdBoundaryLayerInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    atmdBoundaryLayer::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    atmdBoundaryLayerInletVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
