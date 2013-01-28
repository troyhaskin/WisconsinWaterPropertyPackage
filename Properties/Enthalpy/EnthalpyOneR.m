function h = EnthalpyOneR(delta,tau)
    
    R  = SpecificGasConstant();
    Tc = CriticalTemperature();
    
    hStar = R * Tc;
    hND   = EnthalpyOneRND(delta,tau);
    
    h = hND * hStar;
    
end