// Passive scalar transport - water vapour
{
    volScalarField Dw
    (
        IOobject
        (
            "Dw",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence->nut()/Sct_+D_wa
    );
    Dw.correctBoundaryConditions();

    fvScalarMatrix w_wEqn
    (
        fvm::ddt(rho, w_w)
      + fvm::div(phi, w_w)
      - fvm::laplacian(Dw*rho, w_w)
     ==
        fvOptions(rho, w_w)
    );

    w_wEqn.relax();
    fvOptions.constrain(w_wEqn);
    w_wEqn.solve();
    fvOptions.correct(w_w);

//    w_w==max(w_w, 0.0*w_w);
}
