{
    U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);

    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
