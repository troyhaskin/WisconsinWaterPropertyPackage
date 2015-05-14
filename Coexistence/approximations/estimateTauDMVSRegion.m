function tauBotTop = estimateTauDMVSRegion(delta)
    
    %   Values at maximum saturation density
    [~,delMax]   = saturableDeltas();
    [~,tauMax,~] = dmvsTaus()       ;
    
    %   Top correlation
    c      = [+4.9690879530167342E-01;+1.6446327056450776E+00;-5.3583496627473801E+00];
    tauTop = c(2)*(delMax - delta).^c(1).*exp(c(3)*(delMax - delta))+tauMax;
    
    %   Bottom correlation
    c      = [+5.0023513364394190E-01;-1.7048255682253683E+00;+6.1326675555440637E+00];
    tauBot = c(2)*(delMax - delta).^c(1).*exp(c(3)*(delMax - delta))+tauMax;
    
    % Concatenate
    tauBotTop = [tauBot;tauTop];

end