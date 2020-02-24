    // Passive scalar equation for temperature field
    volScalarField alphat("alphat", turbulence->nut()/Prt+turbulence->nu()/Pr);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphat, T)
     ==
        fvOptions(T)
    );

    TEqn.relax();

    TEqn.solve();
