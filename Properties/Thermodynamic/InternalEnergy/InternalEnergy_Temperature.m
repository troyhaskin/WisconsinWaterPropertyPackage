function i_T = InternalEnergy_Temperature(rho,T,PhaseCheck)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if nargin < 3
        PhaseCheck = true;
    end
       
    i_tau = InternalEnergy_tau(rho,T,PhaseCheck);
    tau_T = tau_Temperature(T);
    
    i_T = i_tau .* tau_T;
    i_T = RestoreShape(i_T,GreatestProduct(SizeRho,SizeT));
    
end