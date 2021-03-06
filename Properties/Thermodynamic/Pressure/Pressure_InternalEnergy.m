function P_i = Pressure_InternalEnergy(rho,T,PhaseCheck)

    [rho,SizeRho,T,SizeT] = Columnify(rho,T);

    if nargin < 3
        PhaseCheck = true;
    end

    OnePhiHandle = @(delta,tau,~) PressureOneR_tau(delta,tau) ;
                                 
    P_tau = GenericTwoPhiProperty(rho,T,OnePhiHandle,'quality',PhaseCheck);
    tau_T = tau_Temperature(T);
    i_T   = InternalEnergy_TemperatureN(rho,T,PhaseCheck);
    P_i   = P_tau .* tau_T ./ i_T;
    
    P_i = RestoreShape(P_i,GreatestProduct(SizeRho,SizeT));
end