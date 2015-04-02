function [] = VerifyValues_HelmholtzIdeal(filename)
    
    if (nargin > 0) && not(isempty(filename)) && ischar(filename)
        outputFileName  = [CallerDirectory(),filename];
        fileID          = fopen(outputFileName,'w','native','UTF-8');
        
        if (fileID ~= -1)
            Print = @(Format,varargin) fprintf(fileID,Format,varargin{:});
        else
            fprintf('\n');
            warning('WWPP:VerifyTwoPhaseValues:UnableToOpenFile',...
                ['**** The requested file could not be created, ',...
                'so the verification is being output to the command window. ****']);
            fprintf('\n');
            Print = @(Format,varargin) fprintf(Format,varargin{:});
        end
        
    else
        Print = @(Format,varargin) fprintf(Format,varargin{:});
        fileID = -1;
    end
    
    
    
    %    Set-up state
    Tc    = CriticalTemperature();
    rhoc  = CriticalDensity()    ;
    
    Ttest   = 500       ; %[K]
    rhoTest = 838.025   ; %[kg/m^3]
    
    delta = rhoTest/  rhoc  ;
    tau   =   Tc   / Ttest  ;
    
    
    
    % Load ideal values
    FileNameSinglePhase = 'IAPWSValues_HelmholtzIdeal.csv';
    IAPWS = importdata(FileNameSinglePhase, ',', 1);
    IAPWS = IAPWS.data(:);
    
    %   Calculate package values
    WWPP = [...
        HelmholtzIdealGas(delta,tau)    ;
        HelmholtzIdealGas_d(delta,tau)  ;
        HelmholtzIdealGas_dd(delta,tau) ;
        HelmholtzIdealGas_t(delta,tau)  ;
        HelmholtzIdealGas_tt(delta,tau) ;
        HelmholtzIdealGas_dt(delta,tau)];
    WWPP = SigFigRound(WWPP,9);
    
    %   Difference
    absoluteDifference = abs(IAPWS - WWPP);
    
    
    %   Output comparison
    Print('\n');
    Print('                      Helmholtz Ideal Gas Values            \n');
    Print('===================================================================\n');
    Print('              State: {rho = %+8.3f, T = %+8.3f} \n',rhoTest,Ttest);
    Print('-------------------------------------------------------------------\n');
    Print('                  WWPP               IAPWS             Difference\n');
    Print('-------------------------------------------------------------------\n');
    Print('Phi0      %+17.8E   %+17.8E   %+17.8E\n',WWPP(1),IAPWS(1),absoluteDifference(1));
    Print('Phi0_d    %+17.8E   %+17.8E   %+17.8E\n',WWPP(2),IAPWS(2),absoluteDifference(2));
    Print('Phi0_dd   %+17.8E   %+17.8E   %+17.8E\n',WWPP(3),IAPWS(3),absoluteDifference(3));
    Print('Phi0_t    %+17.8E   %+17.8E   %+17.8E\n',WWPP(4),IAPWS(4),absoluteDifference(4));
    Print('Phi0_tt   %+17.8E   %+17.8E   %+17.8E\n',WWPP(5),IAPWS(5),absoluteDifference(5));
    Print('Phi0_dt   %+17.8E   %+17.8E   %+17.8E\n',WWPP(6),IAPWS(6),absoluteDifference(6));
    Print('===================================================================\n');
    Print('\n');
    
    
    if (fileID ~= -1)
        fclose(fileID);
    end
    
end