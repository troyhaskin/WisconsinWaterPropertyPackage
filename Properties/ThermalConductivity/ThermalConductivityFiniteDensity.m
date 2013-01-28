function lam1 = ThermalConductivityFiniteDensity(delta,tau)
    
    L        = ThermalConductivityFiniteDensityConstants();
    PowerSum = 0;
    SumErrk  = 0;
    
    for k = 1:5
        c       = 0 ;
        SumErrm = 0 ;
        for m = 1:6
            Part        = L(k,m) * (delta - 1).^(m-1);
            [c,SumErrm] = KahanSum(c,Part,SumErrm);
        end
        
        Part               = c .* (tau - 1).^(k-1);
        [PowerSum,SumErrk] = KahanSum(PowerSum,Part,SumErrk);
    end

    lam1 = exp(delta.*PowerSum);
    
end