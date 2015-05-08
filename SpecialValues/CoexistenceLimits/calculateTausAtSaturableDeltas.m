function [] = calculateTausAtSaturableDeltas()
    
    
    %   tau associated with delLo from saturableDeltas()
    tauLo = TriplePointTau();
    
    %   tau associated with delHi from saturableDeltas()
    tauHi = fminsearch(@(tau) -saturationLiquidDensity(tau) ,2.3348,struct('TolX',1E-17,'TolFun',1E-17));
    
    printLoHiRange([tauLo,tauHi],'delta''s tau',false);
    
end

function delL = saturationLiquidDensity(tau)
    [~,delL,~] = SaturationStateGivenTauRRND(tau);
end