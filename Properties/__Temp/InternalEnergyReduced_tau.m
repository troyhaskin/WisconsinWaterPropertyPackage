function didtau = InternalEnergyFromDeltaTau_tau(delta,tau,PhaseCheck)
    
    if nargin < 3
        PhaseCheck = true;
    end
    
    Tc        = CriticalTemperature()   ; %[K]
    R         = SpecificGasConstant()   ; %[J/kg/K]        
            
    OnePhiHandle = @(delta,tau,Mask) SinglePhase_tau(delta,tau,R*Tc);                                
    didtau = GenericTwoPhiPropertyFromDeltaTau(delta,tau,OnePhiHandle,'quality',PhaseCheck);
    
end

function i = SinglePhase_tau(delta,tau,RTc)
    i = (HelmholtzIdealGas_tt(delta,tau)    +...
         HelmholtzResidual_tt(delta,tau))   .* RTc;
end

% function dimaxdtau = MixtureDerivative(rhol,rhog,rho,T,TwoPhase)
%     x = QualityFromDensity(rho,rhol,rhog);
% end
