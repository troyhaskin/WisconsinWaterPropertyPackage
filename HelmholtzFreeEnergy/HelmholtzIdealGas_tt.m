function HelmDeriv = HelmholtzIdealGas_tt(~,tau)
    
    persistent n_o gamma_o
    
    if isempty(n_o)
        [n_o,gamma_o] = HelmholtzIdealGas_Coefficients()    ;
    end
    
    HelmDeriv   = - n_o(3) * tau.^(-2)  ;
    SumErr      = 0                     ;
    
    for k = 4:8
        Part    = exp(-gamma_o(k).*tau);
        Part    = - n_o(k) * gamma_o(k)^2 * Part./(1-Part).^2;
        [HelmDeriv,SumErr] = KahanSum(HelmDeriv,Part,SumErr);
    end
end