function IntE = InternalEnergyOne(rho,T)
    
    [~,rhoc,Tc] = Nondimensionalizers();
    
    delta = rho  / rhoc;
    tau   = Tc  ./ T   ;
    
    IntE = InternalEnergyOneR(delta,tau);
    
end