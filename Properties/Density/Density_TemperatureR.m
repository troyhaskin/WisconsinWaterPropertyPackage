function rho_T = Density_TemperatureR(delta,tau)
    % Isobaric derivative of density w.r.t. temperature  
    
    del_t = delta_tau(delta,tau)    ;
    tau_T = tauRR_Temperature(tau)  ;
    rho_d = Density_delta()         ;
    
    rho_T = del_t .* tau_T .* rho_d ;
    
end