// tmp <volVectorField> relvel = Uc-uParticle;

fvVectorMatrix UcEqn
(
    fvm::ddt(alphac, Uc) + fvm::div(alphaPhic, Uc)
  - fvm::Sp(fvc::ddt(alphac) + fvc::div(alphaPhic), Uc)
  + continuousPhaseTurbulence->divDevRhoReff(Uc)
 ==
    fvm::Sp(uCoeff, Uc) + fvc::Sp(uSourceDrag, uParticle)
);

UcEqn.relax();

volScalarField rAUc(1.0/UcEqn.A());
surfaceScalarField rAUcf("Dp", fvc::interpolate(rAUc));

surfaceScalarField phicForces
(
   fvc::flux(rAUc*uSource) + rAUcf*(g & mesh.Sf())
);

if (pimple.momentumPredictor())
{
    solve
    (
        UcEqn
     ==
        fvc::reconstruct
        (
            phicForces/rAUcf - fvc::snGrad(p)*mesh.magSf()
        )
    );
}
