function P = SoundSpeed(rho,T,PhaseCheck)
    
    if nargin < 3
        PhaseCheck = true;
    end
    
    OnePhiHandle    = @(delta,tau,Mask) SoundSpeedOneR(delta,tau);
    P = GenericTwoPhiProperty(rho,T,OnePhiHandle,'void',PhaseCheck);
    
end
