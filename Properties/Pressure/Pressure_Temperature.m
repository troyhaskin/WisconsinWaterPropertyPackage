function P_T = Pressure_Temperature(rho,T,PhaseCheck)
    
    if nargin < 3
        PhaseCheck = true;
    end
    
    P_tau = Pressure_tau(rho,T,PhaseCheck);
    tau_T = tau_Temperature(T);
    P_T   = P_tau .* tau_T;
    
end