function Psi_xN = PointWiseCentralDifference(PsiHandle,x,Eps)
    
    if (nargin < 3) || isempty(Eps)
        Eps = DefaultDerivativeDelta();
    end
    
    PsipmEps = PsiHandle([x+Eps,x-Eps]);    
    Psi_xN   = (PsipmEps(:,1) - PsipmEps(:,2))./(2*Eps);
    
end