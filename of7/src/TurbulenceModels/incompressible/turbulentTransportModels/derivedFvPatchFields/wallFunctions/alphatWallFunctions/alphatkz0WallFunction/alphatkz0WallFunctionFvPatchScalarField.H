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

Class
    Foam::incompressible::alphatkz0WallFunctionFvPatchScalarField

Description
    This boundary condition provides a kinematic turbulent thermal conductivity
    for using wall functions, using the Jayatilleke 'P' function.

Usage
    \table
        Property     | Description                           | Required    | Default value
        Prt          | turbulent Prandtl number              | no         | 0.85
        Cmu          | model coefficient                     | no          | 0.09
        kappa        | Von Karman constant                   | no          | 0.41
        z0           | Surface roughness lengt               | yes         | N/A
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            alphatkz0WallFunction;
    }
    \endverbatim

Note
    The units of kinematic turbulent thermal conductivity are [m2/s]

See also
    Foam::fixedValueFvPatchField

SourceFiles
    alphatkz0WallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphatkz0WallFunctionFvPatchScalarField_H
#define alphatkz0WallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
       Class alphatkz0WallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatkz0WallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
protected:

    // Protected data

        //- Surface roughness length
        scalar z0_terrain_;

        //- Turbulent Prandtl number
        scalar Prt_;

        //- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;


        // Solution parameters

            static scalar tolerance_;
            static label maxIters_;

public:

    //- Runtime type information
    TypeName("alphatkz0WallFunction");


    // Constructors

        //- Construct from patch and internal field
        alphatkz0WallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphatkz0WallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  alphatkz0WallFunctionFvPatchScalarField
        //  onto a new patch
        alphatkz0WallFunctionFvPatchScalarField
        (
            const alphatkz0WallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        alphatkz0WallFunctionFvPatchScalarField
        (
            const alphatkz0WallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new alphatkz0WallFunctionFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        alphatkz0WallFunctionFvPatchScalarField
        (
            const alphatkz0WallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new alphatkz0WallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
