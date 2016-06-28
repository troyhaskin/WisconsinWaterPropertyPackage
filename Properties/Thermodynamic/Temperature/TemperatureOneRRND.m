function tau = TemperatureOneRRND(delta,iND,tau0)
    
    %   Doesn't account for DMVS region
    if (nargin < 3) || isempty(tau0)
        [~,~,tau0] = EstimateDelLDelGTauFromDel(delta);
    end
    
    %   Iteration set-up
    tolerance = 1E-12   ;
    maxIter   = 100     ;
    updater   = @(tau,mask) update(tau,mask,delta,iND);
    tau       = NewtonUpdater(updater,tau0,tolerance,maxIter);

end

function [dtau,Norm] = update(tau,mask,delta,iNDgiven)
    
    PhiI_t           = HelmholtzIdealGas_t(delta(mask),tau)         ;
    PhiI_tt          = HelmholtzIdealGas_tt(delta(mask),tau)        ;
    [PhiR_t,PhiR_tt] = HelmholtzResidualCombo_t_tt(delta(mask),tau) ;
    iND              = PhiI_t  + PhiR_t                             ;
    iND_tau          = PhiI_tt + PhiR_tt                            ;
    
    residual  = iND - iNDgiven(mask);
    dresidual = iND_tau             ;
    
    dtau = residual ./ dresidual;
    Norm = abs(dtau)            ;
end
