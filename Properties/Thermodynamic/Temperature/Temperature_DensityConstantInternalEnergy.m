function T_rhoConI = Temperature_DensityConstantInternalEnergy(rho,T,PhaseCheck)
    
    i_rho = InternalEnergy_Density(rho,T,PhaseCheck);
    i_T   = InternalEnergy_Temperature(rho,T,PhaseCheck);
    
    T_rhoConI = - i_rho ./ i_T;
    
end