function i_delta = InternalEnergyOne_delta(rho,T)
    
    [~,rhoc,Tc] = Nondimensionalizers();
    
    delta = rho  / rhoc ;
    tau   = Tc  ./ T    ;
    
    i_delta = InternalEnergyOneR_delta(delta,tau);
    
end