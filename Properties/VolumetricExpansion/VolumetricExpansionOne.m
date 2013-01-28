function beta = VolumetricExpansionOne(rho,T)
      
    [~,rhoc,Tc] = Nondimensionalizers();
    
    delta = rho  / rhoc ;
    tau   = Tc  ./ T    ;
    rho_T = Density_TemperatureR(delta,tau) ;
    
    beta = -1./rho .* rho_T;
    
end