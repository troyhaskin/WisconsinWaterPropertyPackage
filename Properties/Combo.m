function varargout = Combo(rho,T,PhaseCheck)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    [rho,T] = BalanceSizes(rho,T);
    
    if (nargin < 3)
        PhaseCheck = true;
    end
    
    OnePhiHandle    = @(delta,tau,~) [PressureOneR(delta,tau),PressureOneR(delta,tau)];
	TwoPhiHandle	= @(Psat,~,~,~,~,~) [Psat,Psat];
    
    P = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhiHandle,PhaseCheck);
    
    varargout = OutputMultipleProperties(P,GreatestProduct(SizeRho,SizeT));
    
end