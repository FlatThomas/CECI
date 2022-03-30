inline dimensionedScalar heavySide(scalar A, const scalar &fCut)
        {
            dimensionedScalar H("H",dimless,0);
            if (A < 1 - fCut || A > fCut)
            {
                H = 0;
            }
            else
            {
                H = 1;
            }
            return H;
        };



