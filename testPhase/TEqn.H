{
  volScalarField::Internal CeSource = rhoSmear * (alpha1()*alpha1() +
                                        alpha2() * alpha2());
  volScalarField::Internal HvSource = rhoDotSharp *(alpha1()/mixture.thermo1().Cv()()+alpha2()/mixture.thermo2().Cv()())* phaseChange.latentHeat(T).ref();

  fvScalarMatrix TEqn
  (
        fvm::ddt(rho, T) + fvm::div(rhoPhi, T) - fvm::Sp(contErr, T)
      - fvm::laplacian(turbulence.alphaEff(), T)
      + (
             mixture.totalInternalEnergy()
          ?
            fvc::div(fvc::absolute(phi, U), p)()() // - contErr/rho*p
          + (fvc::ddt(rho, K) + fvc::div(rhoPhi, K))()()
          - (U()&(fvModels.source(rho, U)&U)()) - contErr*K
          :
            p*fvc::div(fvc::absolute(phi, U))()()
        )
       *(
           alpha1()/mixture.thermo1().Cv()()
         + alpha2()/mixture.thermo2().Cv()()
        )
     ==
      fvm::SuSp(CeSource,T)+HvSource
  );

  TEqn.relax();

  fvConstraints.constrain(TEqn);

  TEqn.solve();

  fvConstraints.constrain(T);

  mixture.correctThermo();
  mixture.correct();
}
