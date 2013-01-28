function mu1 = ViscosityFiniteDensity(delta,tau)
    
    H        = ViscosityFiniteDensityConstants();
    PowerSum = 0;
    SumErrk  = 0;
    
    for k = 1:6
        c       = 0 ;
        SumErrm = 0 ;
        for m = 1:7
            Part        = H(k,m) * (delta - 1).^(m-1);
            [c,SumErrm] = KahanSum(c,Part,SumErrm);
        end
        
        Part               = c .* (tau - 1).^(k-1);
        [PowerSum,SumErrk] = KahanSum(PowerSum,Part,SumErrk);
    end

    mu1 = exp(delta.*PowerSum);
    
end