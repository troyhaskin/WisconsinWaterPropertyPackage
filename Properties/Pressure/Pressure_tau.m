function P_tau = Pressure_tau(rho,T,PhaseCheck)
    
    [rho,SizeRho,T,SizeT] = Columnify(rho,T);
    
    if nargin < 3
        PhaseCheck = true;
    end
            
    OnePhiHandle = @(delta,tau,Mask) PressureOneRND_tau(delta,tau);
    TwoPhiHandle = @(Psat,delL,delG,~,tauSat,~) ClausiusClapeyronRRND(Psat,tauSat,delL,delG);
       
    Pnd_tau = GenericTwoPhiProperty(rho,T,OnePhiHandle,TwoPhiHandle,PhaseCheck);
    
	Pstar = DimensioningPressure();
    P_tau = Pnd_tau * Pstar;
    
    P_tau = RestoreShape(P_tau,GreatestProduct(SizeRho,SizeT));
    
end

% function P_tauSat = TwoPhasePressure_tauND(delL,delG,tauSat,Psat)
%     P_tauSat = ClausiusClapeyronRRND(Psat,tauSat,delL,delG);
% end





