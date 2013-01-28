function Helm = HelmholtzIdealGas(delta,tau)
    
    persistent n_o gamma_o
    
    if isempty(n_o)
        [n_o,gamma_o] = HelmholtzIdealGas_Coefficients()    ;
    end
    
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