function Helm = HelmholtzIdealGas(delta,tau)
    
    [n_o,gamma_o] = Coefficients_HelmholtzIdealGas()    ;
    
    Helm    = log(delta)    ;
    SumErr  = 0             ;
    
    [Helm,SumErr]   = KahanSum(Helm,n_o(1)           ,SumErr);
    [Helm,SumErr]   = KahanSum(Helm,n_o(2) * tau     ,SumErr);
    [Helm,SumErr]   = KahanSum(Helm,n_o(3) * log(tau),SumErr);
        
    for k = 4:8
        Part   = n_o(k) * log(1- exp(-gamma_o(k)*tau));
        [Helm,SumErr] = KahanSum(Helm,Part,SumErr);
    end

end