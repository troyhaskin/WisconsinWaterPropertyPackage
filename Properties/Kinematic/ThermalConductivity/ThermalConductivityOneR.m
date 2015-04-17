function lam = ThermalConductivityOneR(delta,tau)
    
    lamRef = DimensioningThermalConductivity()         ; %[W/m-K]
    lamND  = ThermalConductivityOneRND(delta,tau); %[-]
    
    lam = lamND * lamRef;%[W/m-K]
end