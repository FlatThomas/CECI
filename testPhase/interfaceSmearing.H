//Correct Phasechange properties
Reconstruction.correct();
rhoDotSharp.ref()=phaseChange.RhoDotSharp(alpha1)[0];

//Initialize Volume Integral
scalar volIntegral = 0.0;

// Loop Over Cells to get Value
forAll(rhoDotSharp, CellI)
{
    volIntegral += rhoDotSharp[CellI] * mesh.V()[CellI];
}

/*---------------------------------------------------*\
              Smear Sharp Source Term
    \*---------------------------------------------------*/

// Dummy Initializations, Fix Later
dimensionedScalar DDT = phaseChange.lookup("DDT");

// Initialize smeared field, copied from sharp field
 rhoSmear=rhoDotSharp;
// Create PDE
Foam::fvScalarMatrix smearEqn(
    fvm::Sp(scalar(1), rhoSmear) - fvm::laplacian(DDT, rhoSmear) == rhoDotSharp

);

// Solve!
smearEqn.solve();

/*---------------------------------------------------*\
              Crop Smeared Field
\*---------------------------------------------------*/
const scalar fCut = .001;
volScalarField::Internal leftHeavy(heavySide(-1 * alpha1(), fCut).ref());
volScalarField::Internal rightHeavy(heavySide(alpha1(), -1 + fCut).ref());

forAll(mesh.C(), CellI)
{
    if (alpha1()[CellI] < 1 - fCut && alpha1()[CellI] > fCut)
    {
        rhoSmear[CellI] = 0;
    }
}


/*---------------------------------------------------*\
              Solve For Coefficients
    \*---------------------------------------------------*/

dimensionedScalar heavyVolL(dimless, 0.0);
dimensionedScalar heavyVolV(dimless, 0.0);

forAll(mesh.C(), CellI)
{
    heavyVolL += rightHeavy[CellI] * rhoSmear[CellI] * mesh.V()[CellI];
    heavyVolV += leftHeavy[CellI] * rhoSmear[CellI] * mesh.V()[CellI];
}
Info << "heavyVolL=" << heavyVolL << endl;
Info << "heavyVolV=" << heavyVolV << endl;
dimensionedScalar Nl("Nl", dimless, 0.0);
dimensionedScalar Nv("Nv", dimless, 0.0);
if(heavyVolL.value()>1e-99)
{
    Nl= volIntegral / heavyVolL;
}
if(heavyVolV.value()>1e-99)
{
    Nv= volIntegral / heavyVolV;
}

/*---------------------------------------------------*\
        Calculate Final Source Term Distribution
         \*---------------------------------------------------*/

forAll(mesh.C(), CellI)
{
    rhoSmear[CellI] = Nv.value() * leftHeavy[CellI] * rhoSmear[CellI] - Nl.value() * rightHeavy[CellI] * rhoSmear[CellI];
}

scalar continuityTest = 0;
scalar continuityTest2 = 0;
// Validate Continuity
forAll(rhoSmear, CellI)
{
    continuityTest += rhoSmear[CellI] * mesh.V()[CellI] * leftHeavy[CellI];
    continuityTest2 += rhoSmear[CellI] * mesh.V()[CellI] * rightHeavy[CellI];
}

Info << "Continuity Value: " << continuityTest << endl;
Info << "Sharp Continuity: " << volIntegral << endl;
Info << "Continuity2 Value: " << continuityTest2 << endl<<endl;
