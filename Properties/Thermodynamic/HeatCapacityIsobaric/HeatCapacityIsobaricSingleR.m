function cp = HeatCapacityIsobaricSingleR(delta,tau)
    % Single phase property with reduced state specification
    
    R = SpecificGasConstant(); %[J/kg/K]
    
    cpND = HeatCapacityIsobaricSingleRND(delta,tau);
    cp   = cpND * R;
    
end