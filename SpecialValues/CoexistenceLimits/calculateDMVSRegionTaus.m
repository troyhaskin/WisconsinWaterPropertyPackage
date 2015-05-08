function [] = calculateDMVSRegionTaus
    
    %   Biggest tau
    [delLt,~] = TriplePointDensitiesR() ;
    tauHi     = TriplePointTau()        ;
    
    %   Middle tau
    [~,tauMd] = tausAtSaturableDeltas();
    
    %   Smallest tau
    tauLo = fzero(@(tau) saturationLiquidDensity(tau) - delLt,+2.30026);
    
    printLoHiRange([tauLo;tauMd;tauHi],'tau');
    
end


function delL = saturationLiquidDensity(tau)
    [~,delL,~] = SaturationStateGivenTauRRND(tau);
end