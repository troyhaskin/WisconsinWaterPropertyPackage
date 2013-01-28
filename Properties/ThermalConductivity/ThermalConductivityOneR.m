function lam = ThermalConductivityOneR(delta,tau)
    
    lamRef = ReferenceThermalConductivity()         ; %[W/m-K]
    lamND  = ThermalConductivityOneRND(delta,tau); %[-]
    
    lam = lamND * lamRef;%[W/m-K]
end