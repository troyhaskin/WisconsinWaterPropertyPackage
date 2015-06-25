function P = Pressure(rho,T,PhaseCheck,twoPhiState)

    if (nargin < 4)
        twoPhiState = [];
    end
    if (nargin < 3) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);

    onePhiHandle = @(delta,tau,Mask) PressureOneR(delta,tau)            ;
	twoPhiHandle = @(Pnd,~,~,~,~,TwoPhase) Pnd*DimensioningPressure()   ;

    P = GenericTwoPhiProperty(rho,T,onePhiHandle,twoPhiHandle,PhaseCheck,twoPhiState);
    P = RestoreShape(P,GreatestProduct(SizeRho,SizeT));

end