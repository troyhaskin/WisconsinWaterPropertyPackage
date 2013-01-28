function s = EntropyOne(rho,T)
    
    [~,rhoc,Tc] = Nondimensionalizers();
    
    delta = rho  / rhoc;
    tau   = Tc  ./ T;
    
    s = EntropyOneR(delta,tau);
    
end