// Passive scalar transport
{
    volScalarField Dc
    (
        IOobject
        (
            "Dc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence->nut()/Sct_+D_ca
    );
    Dc.correctBoundaryConditions();

    fvScalarMatrix w_cEqn
    (
        fvm::ddt(w_c)
      + fvm::div(phi, w_c)
      - fvm::laplacian(Dc, w_c)
     ==
        fvOptions(w_c)
    );

    w_cEqn.relax();
    fvOptions.constrain(w_cEqn);
    w_cEqn.solve();
    fvOptions.correct(w_c);

//    w_c==max(w_c, 0.0*w_c);
}
