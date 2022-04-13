//Variable Declarations
Info<<"Declaring Variables"<<endl;
tmp<volScalarField::Internal> trhodotSharp=phaseChange.Salpha(alpha1)[0];
rhoDotSharp.ref()=trhodotSharp.ref();
scalar volIntegral=0.0;

//Loop Over Cells to get Value
forAll(rhoDotSharp, CellI)
{
    volIntegral += rhoDotSharp[CellI] * mesh.V()[CellI];
}
Info<<"volIntegral="<<volIntegral<<endl;

/*---------------------------------------------------*\
              Smear Sharp Source Term
    \*---------------------------------------------------*/
Info<<"You've made it this far (1)!"<<endl;

//Dummy Initializations, Fix Later
dimensionedScalar DDT=phaseChange.lookup("DDT");

//Initialize smeared field, copied from sharp field
//rhoSmear.ref()=rhoDotSharp;
    //Create PDE
    Foam::fvScalarMatrix smearEqn(
        fvm::Sp(scalar(1),rhoSmear)-fvm::laplacian(DDT, rhoSmear) == rhoDotSharp

    );

//Solve!
smearEqn.solve();

/*---------------------------------------------------*\
              Crop Smeared Field
\*---------------------------------------------------*/
const scalar fCut=.001; 
volScalarField::Internal leftHeavy(heavySide(-1*alpha1(), fCut).ref());
volScalarField::Internal rightHeavy(heavySide(alpha1(),-1+fCut).ref());



forAll(mesh.C(), CellI)
{
    if (alpha1()[CellI] < 1 - fCut && alpha1()[CellI] > fCut)
    {
        rhoSmear[CellI] = 0;
    }
}


Info<<"You've made it this far (2)!"<<endl;


/*---------------------------------------------------*\
              Solve For Coefficients
    \*---------------------------------------------------*/

dimensionedScalar heavyVolL(dimless,0.0);
dimensionedScalar heavyVolV(dimless,0.0);

forAll(mesh.C(), CellI)
{
    heavyVolL += rightHeavy[CellI] * rhoSmear[CellI] * mesh.V()[CellI];
    heavyVolV += leftHeavy[CellI] * rhoSmear[CellI] * mesh.V()[CellI];
}
Info<<"You've made it this far (3)!"<<endl;
Info<<"heavyVolL="<<heavyVolL<<endl;
Info<<"heavyVolV="<<heavyVolV<<endl;

dimensionedScalar Nl = volIntegral / heavyVolL;
dimensionedScalar Nv = volIntegral / heavyVolV;

/*---------------------------------------------------*\
        Calculate Final Source Term Distribution
    \*---------------------------------------------------*/

forAll(mesh.C(), CellI)
{
    rhoSmear[CellI] = Nv.value() * leftHeavy[CellI] * rhoSmear[CellI]
                    - Nl.value() * rightHeavy[CellI]* rhoSmear[CellI];
}

rhoSmear.write();
scalar continuityTest=0;
scalar continuityTest2=0;
//Validate Continuity
forAll(rhoSmear,CellI)
{
    continuityTest+=rhoSmear[CellI]*mesh.V()[CellI]*leftHeavy[CellI];
    continuityTest2+=rhoSmear[CellI]*mesh.V()[CellI]*rightHeavy[CellI];
}

Info<<"Continuity Value: "<<continuityTest<<endl;
Info<<"Sharp Continuity: "<<volIntegral<<endl;
Info<<"Continuity2 Value: "<<continuityTest2<<endl;

