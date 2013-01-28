function HelmDeriv = HelmholtzIdealGas_t(~,tau)
    
    persistent n_o gamma_o
    
    if isempty(n_o)
        [n_o,gamma_o] = HelmholtzIdealGas_Coefficients()    ;
    end
    
    HelmDeriv   = n_o(2)    ;
    SumErr      = 0         ;
    
    
    [HelmDeriv,SumErr]  = KahanSum(HelmDeriv,n_o(3) ./tau,SumErr);
    
    for k = 4:8
        Part    = n_o(k) * gamma_o(k) * ((1-exp(-gamma_o(k).*tau)).^(-1) - 1);
        [HelmDeriv,SumErr] = KahanSum(HelmDeriv,Part,SumErr);
    end
end