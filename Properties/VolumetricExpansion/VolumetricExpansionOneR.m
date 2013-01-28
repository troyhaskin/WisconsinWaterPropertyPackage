function beta = VolumetricExpansionOneR(delta,tau)
      
    rhoc  = CriticalDensity()               ;
    rho_T = Density_TemperatureR(delta,tau) ;
    
    beta = -1./(rhoc*delta) .* rho_T;
    
end