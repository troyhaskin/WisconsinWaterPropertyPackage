function [] = determineCoefficientsForZeroReferenceState()

    %   Set-Up
    taut          = TriplePointTau();
    [n_o,gamma_o] = Coefficients_HelmholtzIdealGas();
    [~,delL,~]    = SaturationStateGivenTauRRND(taut);
    
    
    %   Solve for n_o(2) first since it will influence the entropy
    iND    = @(n_o2) Local_HelmholtzIdealGas_t(n_o2,n_o,gamma_o,taut) + HelmholtzResidual_t(delL,taut);
    n_o(2) = fzero(iND,n_o(2));
    
    
    %   Solve for n_o(1)
    PhiR = HelmholtzResidual(delL,taut);
    obj = @(n_o1) (Local_HelmholtzIdealGas(n_o1,n_o,gamma_o,delL,taut) + PhiR);
    n_o(1) = fzero(obj,n_o(1));

    
    %   Check values
    iNDt = iND(n_o(2));
    sNDt = taut*iNDt - obj(n_o(1));
    fprintf('Triple Point Internal Energy: %+23.17E\n',iNDt);
    fprintf('Triple Point Entropy:         %+23.17E\n',sNDt);
    fprintf('New Ideal Gas Constants:\n');
    fprintf('     n_o(1):  %+23.16E\n',n_o(1));
    fprintf('     n_o(2):  %+23.16E\n',n_o(2));
    
    
    coefficientFile = which('Coefficients_HelmholtzIdealGas');
    if not(isempty(coefficientFile))
        text = fileread(coefficientFile);
        text = regexprep(text,'(n\_o\(1\).+?)[\-\+]?8\.320446483\d+E\+\d+',['$1',num2str(n_o(1),'%+23.16E')]);
        text = regexprep(text,'(n\_o\(2\).+?)[\-\+]?6\.683210527\d+E\+\d+',['$1',num2str(n_o(2),'%+23.16E')]);
        text = strrep(text,'%','%%');
        delete(coefficientFile);
        fileID = fopen(coefficientFile,'w');
        
        if (fileID ~= -1)
            fprintf(fileID,text);
            fclose(fileID);
        end
    end
    
end


function HelmDeriv = Local_HelmholtzIdealGas_t(n_o2,n_o,gamma_o,tau)

    HelmDeriv   = n_o2  ;
    SumErr      = 0     ;
    
    
    [HelmDeriv,SumErr]  = KahanSum(HelmDeriv,n_o(3)./tau,SumErr);
    
    for k = 4:8
        %   The default equation from the IAPWS-95 document is this:
        %       Part    = n_o(k) * gamma_o(k) * ((1-exp(-gamma_o(k).*tau)).^(-1) - 1);
        %   However, the equivalent reformulation below is used since it is less 
        %   computationally intensive.
        %
        Part    = n_o(k) * gamma_o(k) ./ (exp(gamma_o(k)*tau) - 1);
        [HelmDeriv,SumErr] = KahanSum(HelmDeriv,Part,SumErr);
    end
end


function Helm = Local_HelmholtzIdealGas(n_o1,n_o,gamma_o,delta,tau)
    
    Helm    = log(delta)    ;
    SumErr  = 0             ;
    
    [Helm,SumErr]   = KahanSum(Helm,n_o1             ,SumErr);
    [Helm,SumErr]   = KahanSum(Helm,n_o(2) * tau     ,SumErr);
    [Helm,SumErr]   = KahanSum(Helm,n_o(3) * log(tau),SumErr);
        
    for k = 4:8
        Part   = n_o(k) * log(1- exp(-gamma_o(k)*tau));
        [Helm,SumErr] = KahanSum(Helm,Part,SumErr);
    end

end