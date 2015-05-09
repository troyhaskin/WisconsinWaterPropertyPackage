function [] = calculateINDsAtSaturableDeltas()
    
    [delGmin,delLmax] = saturableDeltas()       ;
    [tauGmin,tauLmax] = tausAtSaturableDeltas() ;
    
    iNDs = InternalEnergyOneRND([delGmin;delLmax],[tauGmin;tauLmax]);
    
    printLoHiRange(iNDs,'delta''s dimensionless internal energy',false);
    
end