function Psi = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhiOption,PhaseCheck)
    
    if nargin < 5
        PhaseCheck = true;
    end
    
    % -----------------------------------------------
    % begin: Allocation and Setup
    [~,rhoc,Tc] = Nondimensionalizers();
    [rho,T]     = BalanceSizes(rho,T);
    
    % -----------------------------------------------
    % begin: Reduce the Quantities
    
    tau   = Tc ./ T   ;
    delta = rho / rhoc;
    
    % end: Reduce the Quantities

    %	Call reduced core method
    Psi = GenericTwoPhiPropertyR(delta,tau,OnePhiHandle,TwoPhiOption,PhaseCheck);

end
