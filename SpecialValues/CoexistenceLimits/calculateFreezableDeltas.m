function [] = calculateFreezableDeltas()
    
    %   Calculate the densities from the end of the known melting line
    [Pnd,tau] = MeltingLineRND([0,0,0,0,2]);
    [delLt,~] = TriplePointDensitiesR();
    delL      = fzero(@(del) PressureOneRND(del,tau(end)) - Pnd(end),delLt);
    
    %   Display values
    printLoHiRange([delLt,delL],'freezing delta');
end