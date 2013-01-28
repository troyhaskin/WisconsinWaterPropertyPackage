function s_tau = EntropyOne_tau(rho,T)
    
    [~,rhoc,Tc] = Nondimensionalizers();
    
    delta = rho  / rhoc;
    tau   = Tc  ./ T;
    
    s_tau = EntropyOneR_tau(delta,tau);
    
end