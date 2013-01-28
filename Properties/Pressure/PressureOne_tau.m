function P_tau = PressureOne_tau(rho,T)
    
    [~,rhoc,Tc] = Nondimensionalizers();
    
    % Reduced quantities
    delta = rho  / rhoc ;
    tau   = Tc  ./ T    ;
    
    % Make the call with reduced quantities
    P_tau = PressureOneR_tau(delta,tau) ;

end