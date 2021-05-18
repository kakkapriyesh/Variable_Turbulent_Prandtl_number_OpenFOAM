{  Info<< "reading Prtt above alpha"<<endl;

   Info<<turbulence->Prtt()<<endl;
   // Info<< "reading nut"<<endl;
    //Info<<turbulence->nut()<<endl;
    alphat = max(turbulence->nut()/(turbulence->Prtt()+dimensionedScalar("small", 10^-30)),(-1.0*turbulence->nu()/Pr));
  //Info<<"reading alphat" <<endl;
  //Info<< alphat <<endl;
    alphat.correctBoundaryConditions();
 // Info<<"boundary condition alphat" <<endl;
 //  Info<< alphat <<endl;
    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);
 //   Info<< "readinng alphaEff\n" <<endl;
 // Info<< alphaEff <<endl;
Info<<"starting T eqn \n" <<endl;
//Info<< T <<endl;
    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        radiation->ST(rhoCpRef, T)
      + fvOptions(T)
    );
//Info<<"done T eqn \n" <<endl;
    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    radiation->correct();

    fvOptions.correct(T);

    rhok = 1.0 - beta*(T - TRef);
}
