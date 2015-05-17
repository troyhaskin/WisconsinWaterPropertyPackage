function tauBotTop = estimateTauDMVSRegion(delta)
    
    %   Values at maximum saturation density
    [~,delMax]   = saturableDeltas();
    [~,tauMax,~] = dmvsTaus()       ;
    
    %   Top correlation
    c      = [+4.9691120113719317E-01;+1.6446991534509878E+00;-5.4106872917548392E+00];
    tauTop = c(2)*(delMax - delta).^c(1).*exp(c(3)*(delMax - delta))+tauMax;
    
    %   Bottom correlation
    c      = [+5.0023547580458538E-01;-1.7048300739047608E+00;+6.1327419909294481E+00];
    tauBot = c(2)*(delMax - delta).^c(1).*exp(c(3)*(delMax - delta))+tauMax;
    
    % Concatenate
    tauBotTop = [tauBot;tauTop];

end