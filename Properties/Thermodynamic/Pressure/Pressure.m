function P = Pressure(rho,T,PhaseCheck)

    [rho,SizeRho,T,SizeT] = Columnify(rho,T);

    if (nargin < 3)
        PhaseCheck = true;
    end

    OnePhiHandle    = @(delta,tau,Mask) PressureOneR(delta,tau);
	TwoPhiHandle	= @(Psat,~,~,~,~,TwoPhase) SmartMask(Psat,TwoPhase);

    P = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhiHandle,PhaseCheck);
    P = RestoreShape(P,GreatestProduct(SizeRho,SizeT));

end