function Psi = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhiOption,PhaseCheck)
    
    if nargin < 5
        PhaseCheck = true;
    end
    
    % Reduce
    tau   = CriticalTemperature() ./ T  ;
    delta = rho / CriticalDensity()     ;

    %	Call reduced core method
    Psi = GenericTwoPhiPropertyR(delta,tau,OnePhiHandle,TwoPhiOption,PhaseCheck);

end
