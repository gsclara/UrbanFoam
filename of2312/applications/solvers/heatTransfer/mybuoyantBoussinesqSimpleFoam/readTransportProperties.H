singlePhaseTransportModel laminarTransport(U, phi);

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport
);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);

// Laminar Schmidt number
dimensionedScalar Sc_("Sc_", dimless, laminarTransport);

// Turbulent Schmidt number
dimensionedScalar Sct_("Sct_", dimless, laminarTransport);

// Contaminant mass diffusivity
dimensionedScalar D_ca("D_ca", dimArea/dimTime, laminarTransport);

// Vapour mass diffusivity
dimensionedScalar D_wa("D_wa", dimArea/dimTime, laminarTransport);
