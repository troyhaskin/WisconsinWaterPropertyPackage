function [] = calculateSaturableDeltas()
    
    %   Smallest saturable delta
    [~,~,delGt] = SaturationStateGivenTauRRND(TriplePointTau());
    
    
    %   Largest saturable delta
    delL = @(tau)saturationLiquidDensity(tau);
    Tmax = fminsearch(@(tau) -delL(tau) ,CriticalTemperature()/277,struct('TolX',1E-16,'TolFun',1E-1));
    delLmax = delL(Tmax);
    
    printLoHiRange([delGt,delLmax],'delta');
    
end

function delL = saturationLiquidDensity(tau)
    [~,delL,~] = SaturationStateGivenTauRRND(tau);
end