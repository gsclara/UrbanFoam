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

Class
    Foam::timeVaryingInletOutletVelocityFvPatchVectorField

Description
    This velocity inlet/outlet boundary condition is applied to pressure
    boundaries where the pressure is specified.  A zero-gradient condition is
    applied for outflow (as defined by the flux); for inflow, the velocity
    is obtained from the flux with the specified inlet direction.

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
        inletDirection | inlet direction per patch face | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            pressureDirectedInletOutletVelocity;
        phi             phi;
        rho             rho;
        inletDirection  uniform (1 0 0);
        value           uniform 0;
    }
    \endverbatim

Note
    Sign conventions:
    - positive flux (out of domain): apply zero-gradient condition
    - negative flux (into of domain): derive from the flux with specified
      direction

See also
    Foam::mixedFvPatchVectorField

SourceFiles
    timeVaryingInletOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingInletOutletVelocityFvPatchVectorField_H
#define timeVaryingInletOutletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "Function1.H"
#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class timeVaryingInletOutletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingInletOutletVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    // Private data

        //- Flux field name
        word phiName_;

	//- Flow direction
        autoPtr<Function1<vector>> flowDir_;

        //- Direction of the z-coordinate
        vector zDir_;

        //- Von Karman constant
        const scalar kappa_;

        //- Turbulent viscosity coefficient
        const scalar Cmu_;

        //- Reference height
        const scalar Zref_;

        //- Reference velocity
        autoPtr<Function1<scalar>> Uref_;


        //- Surface roughness height
        scalarField z0_;

        //- Minimum coordinate value in z direction
        scalarField zGround_;

        //- Friction velocity
        scalarField Ustar_;



public:

    //- Runtime type information
    TypeName("timeVaryingInletOutletVelocity");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  timeVaryingInletOutletVelocityFvPatchVectorField
        //  onto a new patch
        timeVaryingInletOutletVelocityFvPatchVectorField
        (
            const timeVaryingInletOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingInletOutletVelocityFvPatchVectorField
        (
            const timeVaryingInletOutletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingInletOutletVelocityFvPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingInletOutletVelocityFvPatchVectorField
        (
            const timeVaryingInletOutletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingInletOutletVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Attributes

            //- Return true: this patch field is altered by assignment
            virtual bool assignable() const
            {
                return true;
            }

        // Access

            //- Return the name of phi
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& phiName()
            {
                return phiName_;
            }

            //- Return z-direction
            const vector& zDir() const
            {
                return zDir_;
            }

            //- Return z-direction
            const scalarField& Ustar() const
            {
                return Ustar_;
            }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
