/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::functionObjects::firstCellHeight

Description
    Evaluates and outputs turbulence y+ for models. Values written to
    time directories as field 'firstCellHeight'.

    Example of function object specification:
    \verbatim
    firstCellHeight1
    {
        type        firstCellHeight;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

Usage
    \table
        Property | Description                | Required   | Default value
        type     | type name: firstCellHeight           | yes        |
    \endtable

Note
    Writing field 'firstCellHeight' is done by default, but it can be overridden by
	    defining an empty \c objects list. For details see writeLocalObjects.

	See also
	    Foam::functionObject
	    Foam::functionObjects::fvMeshFunctionObject
	    Foam::functionObjects::logFiles
	    Foam::functionObjects::writeLocalObjects
	    Foam::functionObjects::timeControl

	SourceFiles
	    firstCellHeight.C

	\*---------------------------------------------------------------------------*/

	#ifndef functionObjects_firstCellHeight_H
	#define functionObjects_firstCellHeight_H

	#include "fvMeshFunctionObject.H"
	#include "logFiles.H"
	#include "writeLocalObjects.H"
	#include "volFieldsFwd.H"

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	namespace Foam
	{

	// Forward declaration of classes
	class momentumTransportModel;

	namespace functionObjects
	{

	/*---------------------------------------------------------------------------*\
				  Class firstCellHeight Declaration
	\*---------------------------------------------------------------------------*/

	class firstCellHeight
	:
	    public fvMeshFunctionObject,
	    public logFiles,
	    public writeLocalObjects
	{
	    // Private Member Functions

		//- File header information
		virtual void writeFileHeader(const label i);

		//- Calculate y+
		tmp<volScalarField> calcZfirst
		(
            const momentumTransportModel& turbModel
        );


public:

    //- Runtime type information
    TypeName("firstCellHeight");


    // Constructors

        //- Construct from Time and dictionary
        firstCellHeight
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        firstCellHeight(const firstCellHeight&) = delete;


    //- Destructor
    virtual ~firstCellHeight();


    // Member Functions

        //- Read the firstCellHeight data
        virtual bool read(const dictionary&);

        //- Calculate the firstCellHeight field
        virtual bool execute();

        //- Write the firstCellHeight field
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const firstCellHeight&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
