//Variable Declarations
rhoDotSharp.ref()=phaseChange.Salpha(alpha1)[0].ref();
dimensionedScalar volIntegral("volIntegral", dimensionSet(1, -3, -1, 0, 0, 0,0),0.0);

  
//Loop Over Cells to get Value
forAll(rhoDotSharp, CellI)
{
    volIntegral += rhoDotSharp[CellI] * mesh.V()[CellI];
}

/*---------------------------------------------------*\
              Smear Sharp Source Term
    \*---------------------------------------------------*/

//Dummy Initializations, Fix Later
dimensionedScalar DDT=phaseChange.lookup("DDT");

//Initialize smeared field, copied from sharp field
rhoSmear.ref()=rhoDotSharp;
    //Create PDE
    Foam::fvScalarMatrix smearEqn(
        fvm::Sp(scalar(1),rhoSmear)-fvm::laplacian(DDT, rhoSmear) == rhoDotSharp

    );

//Solve!
solve(smearEqn);

/*---------------------------------------------------*\
              Crop Smeared Field
\*---------------------------------------------------*/
const scalar fCut=.001; 

forAll(mesh.C(), CellI)
{
    if (alpha1()[CellI] < 1 - fCut || alpha1()[CellI] > fCut)
    {
        rhoSmear[CellI] = 0;
    }
}
/*---------------------------------------------------*\
              Solve For Coefficients
    \*---------------------------------------------------*/

dimensionedScalar heavyVolL = 0;
dimensionedScalar heavyVolV = 0;

forAll(mesh.C(), CellI)
{
    heavyVolL += heavySide(alpha1()[CellI] - 1 + fCut,fCut) * rhoSmear[CellI] * mesh.V()[CellI];
    heavyVolV += heavySide(fCut - alpha1()[CellI],fCut) * rhoSmear[CellI] * mesh.V()[CellI];
}

dimensionedScalar Nl = volIntegral / heavyVolL;
dimensionedScalar Nv = volIntegral / heavyVolV;

/*---------------------------------------------------*\
        Calculate Final Source Term Distribution
    \*---------------------------------------------------*/

forAll(mesh.C(), CellI)
{
    rhoSmear[CellI] = Nv.value() * heavySide(fCut - alpha1()[CellI],fCut).value() * rhoSmear[CellI] - Nl.value() * heavySide(alpha1()[CellI] - 1 + fCut,fCut).value() * rhoSmear[CellI];
}


