function iNDt = InternalEnergyTripleRND(delta)

    %   Load triple line data
    [delLt,delGt] = TriplePointDensitiesR();
    [iLt,iGt]     = TriplePointInternalEnergiesND();
    
    %   Quality
    vLt = 1./delLt;
    x = (1./delta - vLt)./(1./delGt - vLt);
    
    %   Evaluate
    iNDt = (1 - x).*iLt + x.*iGt;

    
end