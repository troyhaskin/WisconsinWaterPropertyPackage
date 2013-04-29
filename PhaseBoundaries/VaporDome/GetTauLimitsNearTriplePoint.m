function [tauLow,tauHigh] = GetTauLimitsNearTriplePoint(delL)

    % As long as rhoL is less than the triple point density, delL will uniquely determine
    % a saturation state.  However, the saturation curve becomes multi-valued as a
    % function of density for densities equal to or greater than the triple point density.
    % This function returns the upper and lower limits of tau given a delL in this
    % multi-valued regime via a high order polynomial fit.

    % Store original size and columnify the input vector
    OriginalSize = size(delL)   ;
    delL         = delL(:)      ;
    
    % Load constant values and generate the computer mask
    [delLt,~ ]   = TriplePointDensitiesR()              ;
    delLmax      = MaximumSaturationDensityR()          ;
    ComputeMask  = (delLt <= delL) & (delL <= delLmax)  ;

    % Allocate output vectors
    tauLow  = delL * 0;
    tauHigh = delL * 0;

    if any(ComputeMask)

        % Pull the valid values
        delL = delL(ComputeMask);
        N    = length(delL)     ;

        % Get the tauGuess from solving the lower-order, quadratic fit
        [p,mu] = GuessCoefficients(); 
        [tauHighGuess,tauLowGuess] = QuadraticSolve(p(1),p(2),p(3) - delL,mu);


        % Get higher order fit parameters
        [p,mu]      = HighOrderCoefficients(); % mu is not needed
        
        % Augment delL compute for high/low residuals and create guess vector
        delL     = [    delL     ;     delL   ];
        tauGuess = [tauHighGuess ; tauLowGuess];

        % Solve the system using the Update() function below
        tauSol = NewtonUpdater(@(tau,Mask) Update(tau,Mask,p,mu,delL),tauGuess,1E-14,100,true);

        % Assign converged values
        tauHigh(ComputeMask) = tauSol( 1 : N )  ;
        tauLow (ComputeMask) = tauSol(N+1:end)  ;

    end

    % Reshape to original input's size
    tauHigh = reshape(tauHigh,OriginalSize);
    tauLow  = reshape(tauLow ,OriginalSize);

end

function [dtau,Rnorm] = Update(tau,Mask,p,mu,delL)

    % Calculate the polynomial and its derivative via Horner's Method
    [f,df] = HornersMethodDerivative(tau,p,mu);

    % Newton search direction and convergence metric
    dtau  = (f - delL(Mask))./df    ;
    Rnorm = abs(dtau)               ;

end




function [p,mu] = GuessCoefficients()
    
    % Polynomial coefficient array
    p  = [  -1.3762983879685923E-04 ,...
        +5.9290807879570011E-06 ,...
        +3.1053575036704597E+00 ];
    
    % Shift and scale coefficient array
    mu = [  +2.3342570297432994E+00 ,...
        +1.9850972137560687E-02 ];
end

function [p,mu] = HighOrderCoefficients()
    
    % Polynomial coefficient array
    p  = [  -5.1215324577905310E-14 ,...
        -1.7428143309626436E-12 ,...
        -5.5977119121727325E-11 ,...
        -1.5911064197115862E-09 ,...
        -6.0198659362902838E-08 ,...
        -1.0502615227318029E-06 ,...
        -1.3743471210755714E-04 ,...
        +7.8256398854154390E-06 ,...
        +3.1053574386590754E+00 ];
    
    % Shift and scale coefficient array
    mu = [  +2.3342570297432994E+00 ,...
        +1.9850972137560687E-02 ];
end


function [] = CoefficientGeneration(LowOrder,HighOrder)
    clc;
    % Constant loads
    Tt   = TriplePointTemperature()     ; % Triple point temperature
    Tt2  = 281.31428603704              ; % Temperature where delL == delLt
    Tc   = CriticalTemperature()        ; % Critical point temperature
%     rhoc = CriticalDensity()            ; % Critical density
    [delLt,~] = TriplePointDensitiesR() ; % Triple point density
    
    % Temperature vectors
    T   = linspace(Tt,Tt2,1E3);
    tau = Tc./T;
    
    % Calculate the saturation densities for the given tau
    [~,delL,~] = SaturationStateGivenTausat(tau);
    

    % Lower-order Fit
    Order    = LowOrder;
    [p,~,mu] = polyfit(tau,delL,Order);
    delLfit  = polyval(p,tau,[],mu);
    Error    = delLfit - delL;
    
    fprintf('Statistics for fit of order %G:\n',Order);
    fprintf('    L_1   Residual:          %+10.4E\n',norm(Error, 1 ));
    fprintf('    L_2   Residual:          %+10.4E\n',norm(Error, 2 ));
    fprintf('    L_inf Residual:          %+10.4E\n',norm(Error,Inf));
    fprintf('    Triple point Residuals:  %+10.4E\n',abs(delLfit(1)   - delLt));
    fprintf('                             %+10.4E\n',abs(delLfit(end) - delLt));
    fprintf('    Fit Parameters:\n');
    fprintf('        Shift and scale (mu):\n');
    fprintf('            %+23.16E\n',mu);
    fprintf('        Polynomial coefficients (p):\n');
    fprintf('            %+23.16E\n',p);
    fprintf('\n\n');
    
    % High-order Fit
    Order    = HighOrder;
    [p,~,mu] = polyfit(tau,delL,Order);
    delLfit  = polyval(p,tau,[],mu);
    Error    = delLfit - delL;
    
    fprintf('Statistics for fit of order %G:\n',Order);
    fprintf('    L_1   Residual:          %+10.4E\n',norm(Error, 1 ));
    fprintf('    L_2   Residual:          %+10.4E\n',norm(Error, 2 ));
    fprintf('    L_inf Residual:          %+10.4E\n',norm(Error,Inf));
    fprintf('    Triple point Residuals:  %+10.4E\n',abs(delLfit(1)   - delLt));
    fprintf('                             %+10.4E\n',abs(delLfit(end) - delLt));
    fprintf('    Fit Parameters:\n');
    fprintf('        Shift and scale (mu):\n');
    fprintf('            %+23.16E\n',mu);
    fprintf('        Polynomial coefficients (p):\n');
    fprintf('            %+23.16E\n',p);
end

