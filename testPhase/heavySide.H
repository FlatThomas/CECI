//Basic Heavy Side Function
inline tmp<volScalarField::Internal> heavySide
       (const volScalarField::Internal &A, const scalar &fCut) 
    
{
    //Tmp Class Wrapper
    tmp<volScalarField::Internal> tmpH(
            new volScalarField::Internal(
                IOobject(
                    "H",
                    A.time().timeName(),
                    A.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE),
                A.mesh(),
                dimensionedScalar("H",dimless,0)));
    volScalarField::Internal &H=tmpH.ref();

    //Write Cell Values
       forAll(A,CellI)
       {
           if (A[CellI]+fCut>0)
                {
                    H[CellI] = 1;
                }
       } 
         
   return tmpH;
};



