function P = SoundSpeed(rho,T,PhaseCheck)
    
    if nargin < 3
        PhaseCheck = true;
    end
    
    R     = SpecificGasConstant()   ; %[J/kg/K]
    
    RT              = R .* T;
    OnePhiHandle    = @(delta,tau,Mask) SinglePhase(delta,tau,Mask,RT);
%	TwoPhiHandle	= @(~,~,Psat,~,~,TwoPhase) SmartMask(Psat,TwoPhase);
    
    P = GenericTwoPhiProperty(rho,T,OnePhiHandle,'quality',PhaseCheck);
    
end

function c = SinglePhase(delta,tau,Mask,RT)
    
    Phi0_tt = HelmholtzIdealGas_tt(delta,tau);
    PhiR_d  = HelmholtzResidual_d (delta,tau);
    PhiR_dd = HelmholtzResidual_dd(delta,tau);
    PhiR_dt = HelmholtzResidual_dt(delta,tau);
    PhiR_tt = HelmholtzResidual_tt(delta,tau);

    Part = 1 + delta .* PhiR_r - delta .* tau .* PhiR_dt    ;
    Part = Part ./ (tau.^2 .* (Phi0_tt + PhiR_dd))          ;
    c = 1 + 2* delta.*PhiR_d + delta.^2 .* PhiR_dd - Part   ;   % dimensionless
    c = c .* RT(Mask)                                       ;   % dimensional
end