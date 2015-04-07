function HelmDeriv = HelmholtzIdealGas_t(~,tau)
    
    [n_o,gamma_o] = Coefficients_HelmholtzIdealGas()    ;
    
    HelmDeriv   = n_o(2)    ;
    SumErr      = 0         ;
    
    
    [HelmDeriv,SumErr]  = KahanSum(HelmDeriv,n_o(3) ./tau,SumErr);
    
    for k = 4:8
        %   The default equation from the IAPWS-95 document is this:
        %       Part    = n_o(k) * gamma_o(k) * ((1-exp(-gamma_o(k).*tau)).^(-1) - 1);
        %   However, the equivalent reformulation below is used since it is less 
        %   computationally intensive.
        %
        Part    = n_o(k) * gamma_o(k) ./ (exp(gamma_o(k)*tau) - 1);
        [HelmDeriv,SumErr] = KahanSum(HelmDeriv,Part,SumErr);
    end
end