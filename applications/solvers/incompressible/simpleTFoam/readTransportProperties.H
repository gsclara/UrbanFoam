//    singlePhaseTransportModel transportProperties(U, phi);
Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


// Laminar Prandtl number
dimensionedScalar Pr(transportProperties.lookup("Pr"));

// Turbulent Prandtl number
dimensionedScalar Prt(transportProperties.lookup("Prt"));
