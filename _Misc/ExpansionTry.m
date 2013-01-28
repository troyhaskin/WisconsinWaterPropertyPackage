function rho_T = ExpansionTry(rho,T)
    
    [~,rhoc,Tc] = Nondimensionalizers();
    
    delta = rho  / rhoc ;
    tau   = Tc  ./ T    ;
    
    P_delta = PressureOneR_delta(delta,tau);    
    P_tau   = PressureOneR_tau  (delta,tau);
    
    rho_T = P_tau ./ P_delta .* (rhoc./Tc) .* tau.^2 * (-1./rho);
    
end