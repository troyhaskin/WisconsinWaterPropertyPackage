function HelmDeriv = HelmholtzIdealGas_tt(~,tau)
    
    persistent n_o gamma_o
    
    if isempty(n_o)
        [n_o,gamma_o] = HelmholtzIdealGas_Coefficients()    ;
    end
    
    HelmDeriv   = - n_o(3) * tau.^(-2)  ;
    SumErr      = 0                     ;
    
    for k = 4:8
        %   The default equation from the IAPWS-95 document is this:
        %       Part = - n_o(k) * gamma_o(k)^2 * exp(-gamma_o(k).*tau)./(1-exp(-gamma_o(k).*tau)).^2;
        %   However, the equivalent reformulation below is used since it is less 
        %   computationally intensive.
        %
        Part    = 0.5 * n_o(k) * gamma_o(k)^2 ./(1 - cosh(gamma_o(k) * tau));
        [HelmDeriv,SumErr] = KahanSum(HelmDeriv,Part,SumErr);
    end
end