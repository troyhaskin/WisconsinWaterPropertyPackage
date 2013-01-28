function P_rho = Pressure_Density(rho,T,PhaseCheck)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if (nargin < 3)
        PhaseCheck = true;
    end
    
    % Get some Constants
    [~,rhoc,Tc] = Nondimensionalizers();
    
    % Reduce the variables
    delta = rho  / rhoc ;
    tau   = Tc  ./ T    ;
    
    % Get the derivative pieces
    P_delta   = PressureR_delta(delta,tau,PhaseCheck);
    delta_rho = delta_Density();
    
    % Chain rule
    P_rho = P_delta .* delta_rho;
    
    % Reshape and output
    P_rho = RestoreShape(P_rho,GreatestProduct(SizeRho,SizeT));
    
end

