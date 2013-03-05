function Psi_xN = PointWiseCentralDifference(PsiHandle,x,Eps)
    
    if (nargin < 3) || isempty(Eps)
        Eps = DefaultDerivativeDelta();
    end
    
    PsipmEps = PsiHandle([x+Eps,x-Eps]);    
    Psi_xN   = (PsipmEps(:,1) - PsipmEps(:,2))./(2*Eps);
%     PsipmEps = PsiHandle([(1+Eps)*x,(1-Eps)*x]);  
%     Psi_xN   = (PsipmEps(:,1) - PsipmEps(:,2))./(2*Eps*x);
    
end