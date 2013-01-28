function P_delta = PressureR_delta(delta,tau,PhaseCheck)
    
    [delta,SizeDelta,tau,SizeTau] = Columnify(delta,tau);
    
    if (nargin < 3)
        PhaseCheck = true;
    end
    
    OnePhiHandle    = @(delta,tau,~) PressureOneR_delta(delta,tau);
	TwoPhiHandle	= @(~,~,~,~,~,TwoPhase) SmartMask(0,TwoPhase);
    
    P_delta = GenericTwoPhiPropertyR(delta,tau,OnePhiHandle,TwoPhiHandle,PhaseCheck);
    P_delta = RestoreShape(P_delta,GreatestProduct(SizeDelta,SizeTau));

end

