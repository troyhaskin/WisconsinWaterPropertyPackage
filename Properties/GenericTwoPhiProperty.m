function Psi = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhiOption,PhaseCheck,twoPhiState)
    
    if (nargin < 6)
        twoPhiState = [];
    end

    if (nargin < 5)
        PhaseCheck = true;
    end
    
    % Reduce
    tau   = CriticalTemperature() ./ T  ;
    delta = rho / CriticalDensity()     ;

    %	Call reduced core method
    Psi = GenericTwoPhiPropertyR(delta,tau,OnePhiHandle,TwoPhiOption,PhaseCheck,twoPhiState);

end
