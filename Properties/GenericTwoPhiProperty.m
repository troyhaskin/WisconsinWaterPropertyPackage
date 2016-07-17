function Psi = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhiOption,PhaseCheck,state)
    
    if (nargin < 6) || isempty(state)
        state = [];
    else
        state = state.ND_;
    end

    if (nargin < 5)
        PhaseCheck = true;
    end
    
    % Reduce
    tau   = CriticalTemperature() ./ T  ;
    delta = rho / CriticalDensity()     ;

    %	Call reduced core method
    Psi = GenericTwoPhiPropertyR(delta,tau,OnePhiHandle,TwoPhiOption,PhaseCheck,state);

end
