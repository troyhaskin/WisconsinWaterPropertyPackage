function P_delta = Pressure_delta(rho,T,PhaseCheck)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3)
        PhaseCheck = true;
    end
    
    [~,rhoc,Tc] = Nondimensionalizers();
    
    delta = rho  / rhoc ;
    tau   = Tc  ./ T    ;
    
    P_delta = PressureR_delta(delta,tau,PhaseCheck);
    P_delta = RestoreShape(P_delta,GreatestProduct(SizeRho,SizeT));
end

