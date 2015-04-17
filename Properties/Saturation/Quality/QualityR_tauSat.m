function x_tauSat = QualityR_tauSat(tau,delMix,Psat,delL,delG,il,ig,delL_tau,delG_tau)
    
    % If the saturation state was not specified, calculate it
    if (nargin < 3)
        [Psat,delL,delG] = SaturationStateGivenTausat(tau);
    end
    
    % If the internal energies were not specified, calculate them.
    if (nargin < 8)
        il = InternalEnergyOneR(delL,tau);
        ig = InternalEnergyOneR(delG,tau);
    end
       
    %   In the two-phase region, density and temperature are not independent
    %   variables; therefore, the densities have a derivative w.r.t. tau.
    if (nargin < 8)
        % Get the derivatives
        delL_tau = SatLDensityRND_tau(tau,Psat,delL,delG,il,ig);
        delG_tau = SatGDensityRND_tau(tau,Psat,delL,delG,il,ig);
%         delL_tau = delta_tau(delL,tau);
%         delG_tau = delta_tau(delG,tau);
    end
    
    
    Top1     = (delL - delMix) .* delL .*delG_tau;
    Top2     = (delMix - delG) .* delG .*delL_tau;
    x_tauSat = (Top1 + Top2) ./ (delMix .* (delG - delL).^2);
        
end


